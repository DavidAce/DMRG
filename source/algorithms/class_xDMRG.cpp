//
// Created by david on 2018-02-09.
//


#include <simulation/nmspc_settings.h>
#include <state/class_state_finite.h>
#include <tools/nmspc_tools.h>
#include <tools/finite/opt.h>
#include "class_xDMRG.h"


class_xDMRG::class_xDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
        : class_algorithm_finite(std::move(h5ppFile_), "xDMRG",SimulationType::xDMRG, settings::xdmrg::num_sites) {
    log->trace("Constructing class_xDMRG");
    settings::xdmrg::min_sweeps = std::max(settings::xdmrg::min_sweeps, (size_t)(std::log2(chi_max())));
}




void class_xDMRG::run_preprocessing() {
    log->info("Running {} preprocessing", sim_name);
    class_algorithm_finite::run_preprocessing();
    t_pre.tic();
    sim_status.energy_dens_target = settings::xdmrg::energy_density_target;
    sim_status.energy_dens_window = settings::xdmrg::energy_density_window;
    find_energy_range();
    tools::finite::mps::internals::seed_state_unused = true;
    reset_to_random_state_in_energy_window(settings::model::initial_parity_sector,false, "Initializing");
//    reset_to_random_state_in_energy_window("random",false, "Initializing");
//    inflate_initial_state();

    t_pre.toc();
    log->info("Finished {} preprocessing", sim_name);
}

void class_xDMRG::run_simulation()    {
    log->info("Starting {} simulation", sim_name);
    while(true) {
        log->trace("Starting step {}, iteration {}, direction {}", sim_status.step, sim_status.iteration, state->get_direction());
        single_xDMRG_step();
        write_measurements();
        write_state();
        write_status();
        write_logs();
        check_convergence();
        try_projection();
        backup_best_state(*state); //Should come after check_convergence
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if (state->position_is_any_edge())
        {
            if (sim_status.iteration >= settings::xdmrg::max_sweeps) {stop_reason = StopReason::MAX_ITERS; break;}
            if (sim_status.simulation_has_succeeded)                 {stop_reason = StopReason::SUCCEEDED; break;}
            if (sim_status.simulation_has_to_stop)                   {stop_reason = StopReason::SATURATED; break;}
        }

        log->trace("Finished step {}, iteration {}, direction {}", sim_status.step, sim_status.iteration, state->get_direction());

        move_center_point();
        update_bond_dimension_limit();
        sim_status.iteration     = state->get_sweeps();
        sim_status.position      = state->get_position();
        sim_status.moves         = state->get_moves();
        sim_status.step++;
    }
    switch(stop_reason){
        case StopReason::MAX_ITERS : log->info("Finished {} simulation -- reason: MAX ITERS",sim_name) ;break;
        case StopReason::SUCCEEDED : log->info("Finished {} simulation -- reason: SUCCEEDED",sim_name) ;break;
        case StopReason::SATURATED : log->info("Finished {} simulation -- reason: SATURATED",sim_name) ;break;
        default: log->info("Finished {} simulation -- reason: NONE GIVEN",sim_name);
    }
}




void class_xDMRG::single_xDMRG_step()
{
    using namespace tools::finite;

    t_run.tic();
    log->trace("Starting single xDMRG step");
//  log->debug("Variance accurate check before xDMRG step: {:.16f}", std::log10(measure::accurate::energy_variance_per_site(*state)));

//    auto optMode  = sim_status.iteration  < 2  ?  opt::OptMode::OVERLAP : opt::OptMode::VARIANCE;
//    auto optMode  = sim_status.iteration  < 2  ?  opt::OptMode::OVERLAP : opt::OptMode::VARIANCE;

    auto optMode    = opt::OptMode(opt::MODE::OVERLAP);
    auto optSpace   = opt::OptSpace(opt::SPACE::DIRECT);
    auto optType    = opt::OptType(opt::TYPE::CPLX);
    optMode         = sim_status.iteration  >= 2 and measure::energy_variance_per_site(*state) < 1e-2 ? opt::MODE::VARIANCE   : optMode.option;
    optMode         = sim_status.iteration  >= 5 or  measure::energy_variance_per_site(*state) < 1e-2 ? opt::MODE::VARIANCE   : optMode.option;
//    optSpace        = state->get_chi_lim() <= 16                                                      ? opt::SPACE::SUBSPACE  : optSpace.option;
    optSpace        = optMode == opt::MODE::OVERLAP                                                   ? opt::SPACE::SUBSPACE  : optSpace.option;
    optSpace        = sim_status.simulation_has_got_stuck                                             ? opt::SPACE::SUBSPACE  : optSpace.option;
    optSpace        = state->size_2site()  > settings::precision::maxSizePartDiag                     ? opt::SPACE::DIRECT    : optSpace.option;
    optSpace        = sim_status.variance_mpo_has_converged                                           ? opt::SPACE::DIRECT    : optSpace.option;
    optType         = state->isReal()                                                                 ? opt::TYPE::REAL       : optType.option;
    long threshold = 0;
    switch(optSpace.option){
        case  opt::SPACE::SUBSPACE : threshold = settings::precision::maxSizePartDiag; break;
        case  opt::SPACE::DIRECT   : threshold = settings::precision::maxSizeDirect  ; break;
    }

    Eigen::Tensor<Scalar,3> theta;

    std::list<size_t> max_num_sites_list = {2,settings::precision::maxSitesMultiDmrg};
    while (max_num_sites_list.front() >=  max_num_sites_list.back() and not max_num_sites_list.size()==1) max_num_sites_list.pop_back();

//    if(optSpace.option == opt::SPACE::DIRECT)  max_num_sites_list = {settings::precision::maxSitesMultiDmrg};
//    std::list<size_t> max_num_sites_list = {2,settings::precision::maxSitesMultiDmrg};
//    if(sim_status.iteration == 0) max_num_sites_list = {settings::precision::maxSitesMultiDmrg}; //You can take many sites in the beginning

    for (auto & max_num_sites : max_num_sites_list){
        auto old_num_sites = state->active_sites.size();
        auto old_prob_size = state->active_problem_size();

        if (max_num_sites > 2){
            optSpace  = opt::OptSpace::DIRECT;
            threshold = settings::precision::maxSizeDirect;
        }

        state->activate_sites(threshold, max_num_sites);

        if( state->active_sites.size()   == old_num_sites and
            state->active_problem_size() == old_prob_size){
            if(optSpace == opt::OptSpace::SUBSPACE){
                log->debug("Changing to DIRECT optimization to activate more sites");
                optSpace = opt::OptSpace::DIRECT;
                threshold = settings::precision::maxSizeDirect;
//                max_num_sites_list = {settings::precision::maxSitesMultiDmrg};
                state->activate_sites(threshold, settings::precision::maxSitesMultiDmrg);
            }else{
                log->debug("Keeping old theta: Can't activate more sites");
                theta = state->get_multitheta();
                break;
            }
        }
        if(optSpace ==  opt::OptSpace::SUBSPACE and max_num_sites > 2) log->warn("About to do subspace with too many sites!");
        theta = opt::find_excited_state(*state, sim_status, optMode, optSpace,optType);
        if(optSpace ==  opt::OptSpace::DIRECT){
            double variance_direct   = measure::energy_variance_per_site(*state,theta);
            double variance_old      = measure::energy_variance_per_site(*state);
            if (variance_direct < (1.0-1e-3) * variance_old ) {
                log->debug("Keeping DIRECT optimized state");
                state->tag_active_sites_have_been_updated(true);
            }else{
                log->debug("DIRECT optimization did not improve enough. Keep trying.");
                state->tag_active_sites_have_been_updated(false);
            }

        }

        if(optSpace == opt::OptSpace::SUBSPACE){
            if(optMode == opt::OptMode::OVERLAP){
                log->debug("Keeping state with best overlap");
                state->tag_active_sites_have_been_updated(true);
            }else{
                // Check if you ended up with a better state
                double variance_new   = measure::energy_variance_per_site(*state,theta);
                double variance_old   = measure::energy_variance_per_site(*state);
                if (variance_new >= variance_old){
                    // State got worse.
                    log->debug("State got worse during SUBSPACE optimization");
                    if (sim_status.simulation_has_got_stuck){
                        //  Keep the bad state anyway (use this state as a perturbation to jump out of local minima)
                        log->debug("Keeping state anyway due to saturation");
                        state->tag_active_sites_have_been_updated(true);
                    }else{
                        // Check what DIRECT optimization has to offer
                        log->debug("Checking what DIRECT optimization can achieve");

                        auto theta_direct        = opt::find_excited_state(*state, sim_status, optMode, opt::OptSpace(opt::SPACE::DIRECT),optType);
                        double variance_direct   = measure::energy_variance_per_site(*state,theta_direct);
                        if (variance_direct < (1.0-1e-3) * variance_old ){
                            log->debug("Keeping DIRECT optimized state");
                            state->tag_active_sites_have_been_updated(true);
                            theta = theta_direct;
                        }else{
                            log->debug("DIRECT optimization did not improve enough. Try more sites");
                            state->tag_active_sites_have_been_updated(false);
                        }
                    }
                }else{
                    log->debug("State got better during SUBSPACE optimization");
                    state->tag_active_sites_have_been_updated(true);
                }
            }
        }


        if(state->active_sites_updated()){
            break;
        }
        if(& max_num_sites == &max_num_sites_list.back()){
            log->debug("Keeping last theta: Failed to find better theta and reached max number of sites");
            break;
        }
    }


//    if (optMode == opt::OptMode::OVERLAP){
//        sim_status.chi_temp = 4 * (1+sim_status.iteration);
//        log->debug("Forcing bond dimension down to {}",sim_status.chi_temp);
//    }


    log->debug("Variance check before truncate            : {:.16f}", std::log10(measure::energy_variance_per_site(*state,theta)));
//
    opt::truncate_theta(theta, *state);
//    tools::finite::mps::rebuild_environments(*state);
//    state->unset_measurements();
    log->debug("Variance check after truncate + move + env: {:.16f}", std::log10(measure::energy_variance_per_site(*state)));

    if(std::abs(tools::finite::measure::norm(*state) - 1.0) > settings::precision::maxNormError){
        tools::log->warn("Norm too large: {:.18f}",tools::finite::measure::norm(*state) );
        tools::finite::mps::normalize(*state);
        tools::finite::mps::rebuild_environments(*state);
    }
    if (state->position_is_the_left_edge()){
        tools::finite::mpo::reduce_mpo_energy(*state);
        log->debug("Variance check after reduce          : {:.16f}", std::log10(measure::energy_variance_per_site(*state)));
    }
    debug::check_integrity(*state);

    sim_status.energy_dens        = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);



    t_run.toc();
    sim_status.wall_time = t_tot.get_age();
    sim_status.simu_time = t_run.get_measured_time();

}


void class_xDMRG::check_convergence(){
    t_con.tic();
    if(state->position_is_any_edge()){
        check_convergence_variance(); // HEY! THIS IS NEW --- Trying state_backup for smoother slopes
        check_convergence_entg_entropy();
    }

//    double variance_new  = tools::finite::measure::energy_variance_per_site(*state);
//    double variance_best = tools::finite::measure::energy_variance_per_site(*state_backup);
//    if(variance_new < 0.1 * variance_best and sim_status.simulation_has_got_stuck){
//        log->debug("Simulation got unstuck!");
//        sim_status.simulation_has_got_stuck = false;
//        if(variance_new > settings::precision::varianceConvergenceThreshold){
//            sim_status.variance_mpo_has_saturated = false;
//            sim_status.variance_mpo_saturated_for = 0;
//        }
//    }


    sim_status.energy_dens = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);
    bool outside_of_window = std::abs(sim_status.energy_dens - sim_status.energy_dens_target)  > sim_status.energy_dens_window;
    if (sim_status.iteration > 2)
    {
        if (    outside_of_window
            and (sim_status.variance_mpo_has_saturated or
                 sim_status.variance_mpo_has_converged or
                 tools::finite::measure::energy_variance_per_site(*state) < 1e-4))
        {
            sim_status.energy_dens_window = std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5);
            std::string reason = fmt::format("saturated outside of energy window. Energy density: {}, Energy window: {} --> {}",
                    sim_status.energy_dens, sim_status.energy_dens_window, std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5) );
            reset_to_random_state_in_energy_window(settings::model::initial_parity_sector, false, reason);
        }
//        else
//        if( not     state->all_sites_updated()
//            and     sim_status.simulation_has_got_stuck
//            and     tools::finite::measure::energy_variance_per_site(*state) > 1e-4)
//        {
//            sim_status.energy_dens_window = std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5);
//            std::string reason = fmt::format("could not update all sites. Energy density: {}, Energy window: {} --> {}",
//                     sim_status.energy_dens, sim_status.energy_dens_window, std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5) );
//            reset_to_random_state_in_energy_window(settings::model::initial_parity_sector, false, reason);
//        }
    }

    sim_status.simulation_has_converged = sim_status.variance_mpo_has_converged and
                                          sim_status.entanglement_has_converged;

    sim_status.simulation_has_saturated = (sim_status.variance_mpo_saturated_for >= min_saturation_iters and
                                           sim_status.entanglement_saturated_for >= min_saturation_iters) or
                                          (sim_status.variance_mpo_saturated_for >= max_saturation_iters  or
                                           sim_status.entanglement_saturated_for >= max_saturation_iters)   ;


    sim_status.simulation_has_succeeded = sim_status.simulation_has_converged and
                                          sim_status.simulation_has_saturated;


    sim_status.simulation_has_got_stuck = sim_status.simulation_has_saturated and not
                                          sim_status.simulation_has_succeeded;


    if(state->position_is_any_edge()) {
        sim_status.simulation_has_stuck_for = sim_status.simulation_has_got_stuck ? sim_status.simulation_has_stuck_for + 1 : 0;
    }

        sim_status.simulation_has_to_stop = sim_status.simulation_has_stuck_for >= max_stuck_iters or
                                           (sim_status.simulation_has_stuck_for >= min_stuck_iters and
                                            sim_status.chi_lim_has_reached_chi_max);

    //                                        and (sim_status.variance_mpo_saturated_for >= max_saturation_iters and
    //                                             sim_status.entanglement_saturated_for >= max_saturation_iters);
    //    unstuck = unstuck == true  and  sim_status.simulation_has_got_stuck == false; // We were stuck, but no longer.
    //    if (unstuck) has_projected = false;


    log->debug("Simulation has converged: {}", sim_status.simulation_has_converged);
    log->debug("Simulation has saturated: {}", sim_status.simulation_has_saturated);
    log->debug("Simulation has succeeded: {}", sim_status.simulation_has_succeeded);
    log->debug("Simulation has got stuck: {}", sim_status.simulation_has_got_stuck);
    log->debug("Simulation has stuck for: {}", sim_status.simulation_has_stuck_for);
    log->debug("Simulation has to stop  : {}", sim_status.simulation_has_to_stop);


//    if (    sim_status.num_resets < settings::precision::maxResets
//            and tools::finite::measure::energy_variance_per_site(*state) > 1e-10)
//    {
//        std::string reason = fmt::format("simulation has saturated with bad precision",
//                                         sim_status.energy_dens, sim_status.energy_dens_window, sim_status.energy_dens_window);
//        reset_to_random_state_in_energy_window(settings::model::initial_parity_sector, false, reason);
//    }



    t_con.toc();
}

void class_xDMRG::try_projection(){
    if(settings::model::projection_trial_when_stuck and
        sim_status.simulation_has_got_stuck and
        state->position_is_any_edge())
    {
        log->debug("Trying projection to {}", settings::model::target_parity_sector);
        auto state_projected = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::model::target_parity_sector);
        double variance_projected = tools::finite::measure::energy_variance_per_site(state_projected);
        double variance_original  = tools::finite::measure::energy_variance_per_site(*state);
        if (variance_projected < variance_original){
            log->info("Projection succeeded, variance improved {} -> {}",
                      std::log10(variance_original), std::log10(variance_projected));
            *state = state_projected;
        }else{
            log->info("Projection did not improve variance {} -> {}",
                      std::log10(variance_original), std::log10(variance_projected));
        }
    }

}


void class_xDMRG::inflate_initial_state(){
    tools::log->trace("Inflating bond dimension");
    // Inflate by projecting randomly. Each projection doubles the bond dimension
    for (int i = 0; i < 4; i++){
        *state = tools::finite::ops::get_projection_to_closest_parity_sector(*state, "random" );
        log->debug("χ = {}"         , tools::finite::measure::bond_dimensions(*state));
    }
    *state = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::model::initial_parity_sector);
}


void class_xDMRG::reset_to_random_state_in_energy_window(const std::string &parity_sector,bool inflate, std::string reason ){
    log->info("Resetting to product state -- Reason: {}", reason);
//    sim_status.energy_dens_window = std::min(growth_factor*sim_status.energy_dens_window, 0.5);

    int counter = 0;
    bool outside_of_window = true;
    log->info("Energy initial (per site) = {:.16f} | density = {:.8f} | retries = {}", tools::finite::measure::energy_per_site(*state), sim_status.energy_dens,counter );


    while(outside_of_window){
        reset_to_random_state(parity_sector);
        if (inflate) inflate_initial_state();

        sim_status.energy_dens = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);
        outside_of_window = std::abs(sim_status.energy_dens - sim_status.energy_dens_target)  >= sim_status.energy_dens_window;
        counter++;
        if (counter % 10 == 0 and energy_window_growth_factor != 1.0) {
            log->info("Resetting to product state -- can't find state in energy window.  Increasing energy window: {} --> {}",
                    sim_status.energy_dens_window, std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5) );
            sim_status.energy_dens_window = std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5);
        }
        if(counter >= 2000) throw std::runtime_error("Failed to find initial state in energy window");
    }
    log->info("Energy initial (per site) = {:.16f} | density = {:.8f} | retries = {}", tools::finite::measure::energy_per_site(*state), sim_status.energy_dens,counter );
    clear_saturation_status();
    has_projected   = false;
    sim_status.num_resets++;
    sim_status.energy_ubound      = sim_status.energy_target + sim_status.energy_dens_window * (sim_status.energy_max-sim_status.energy_min);
    sim_status.energy_lbound      = sim_status.energy_target - sim_status.energy_dens_window * (sim_status.energy_max-sim_status.energy_min);

    log->info("Number of product state resets: {}",sim_status.num_resets );
}


void class_xDMRG::find_energy_range() {
    log->trace("Finding energy range");
    if (state->get_length() != num_sites()) throw std::runtime_error("find_energy_range: state lenght mismatch");
    size_t max_sweeps_during_f_range = 4;
    sim_status.iteration = state->reset_sweeps();
    sim_status.moves      = state->reset_moves();
    reset_to_random_state("random");
    update_bond_dimension_limit(16);

    // Find energy minimum
    while(true) {
        class_algorithm_finite::single_DMRG_step("SR");
        print_status_update();
        // It's important not to perform the last moves.
        // That last state would not get optimized
        if(state->position_is_any_edge()){
            if(sim_status.iteration >= max_sweeps_during_f_range
               or tools::finite::measure::energy_variance_per_site(*state) < 1e-8)
            {break;}
        }
        move_center_point();
        sim_status.iteration = state->get_sweeps();

    }
    sim_status.energy_min = tools::finite::measure::energy_per_site(*state);
    reset_to_random_state("random");
    // Find energy maximum
    while(true) {
        class_algorithm_finite::single_DMRG_step("LR");
        print_status_update();
        // It's important not to perform the last moves.
        // That last state would not get optimized
        if(state->position_is_any_edge()){
            if(sim_status.iteration >= max_sweeps_during_f_range
               or tools::finite::measure::energy_variance_per_site(*state) < 1e-8)
            {break;}
        }

        move_center_point();
        sim_status.iteration = state->get_sweeps();
    }
    sim_status.energy_max         = tools::finite::measure::energy_per_site(*state);
    sim_status.energy_target      = sim_status.energy_min    + sim_status.energy_dens_target  * (sim_status.energy_max - sim_status.energy_min);
    sim_status.energy_ubound      = sim_status.energy_target + sim_status.energy_dens_window  * (sim_status.energy_max - sim_status.energy_min);
    sim_status.energy_lbound      = sim_status.energy_target - sim_status.energy_dens_window  * (sim_status.energy_max - sim_status.energy_min);
    log->info("Energy minimum (per site) = {}", sim_status.energy_min);
    log->info("Energy maximum (per site) = {}", sim_status.energy_max);
    log->info("Energy target  (per site) = {}", sim_status.energy_target);
    log->info("Energy lbound  (per site) = {}", sim_status.energy_lbound);
    log->info("Energy ubound  (per site) = {}", sim_status.energy_ubound);

}




bool   class_xDMRG::sim_on()              {return settings::xdmrg::on;}
long   class_xDMRG::chi_max()             {return settings::xdmrg::chi_max;}
size_t class_xDMRG::num_sites()           {return settings::xdmrg::num_sites;}
size_t class_xDMRG::write_freq()          {return settings::xdmrg::write_freq;}
size_t class_xDMRG::print_freq()          {return settings::xdmrg::print_freq;}
bool   class_xDMRG::chi_grow()            {return settings::xdmrg::chi_grow;}
long   class_xDMRG::chi_init()            {return settings::xdmrg::chi_init;}
bool   class_xDMRG::store_wave_function() {return settings::xdmrg::store_wavefn;}



