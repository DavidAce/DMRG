//
// Created by david on 2018-02-09.
//

#include "class_xDMRG.h"
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_finite.h>
#include <tools/finite/mps.h>
#include <tools/finite/mpo.h>
#include <tools/finite/opt.h>
#include <tools/finite/ops.h>
#include <tools/finite/measure.h>
#include <tools/finite/debug.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

class_xDMRG::class_xDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
        : class_algorithm_finite(std::move(h5ppFile_), "xDMRG",SimulationType::xDMRG, settings::xdmrg::num_sites) {
    log->trace("Constructing class_xDMRG");
    settings::xdmrg::min_sweeps = std::max(settings::xdmrg::min_sweeps, (size_t)(std::log2(chi_max())));
}




void class_xDMRG::run_preprocessing() {
    log->info("Running {} preprocessing", sim_name);
    class_algorithm_finite::run_preprocessing();
    tools::common::profile::t_pre->tic();
    sim_status.energy_dens_target = settings::xdmrg::energy_density_target;
    sim_status.energy_dens_window = settings::xdmrg::energy_density_window;
    find_energy_range();
    state->set_chi_max(chi_max());
    sim_status.chi_max = chi_max();
    update_bond_dimension_limit(chi_init());
    if(settings::model::state_number >= 0){
        reset_to_initial_state();
    }else{
        reset_to_random_state_in_energy_window(settings::model::initial_parity_sector,false, "Initializing");
    }
    auto spin_components = tools::finite::measure::spin_components(*state);
    log->info("Initial spin components: {}", spin_components);

    tools::common::profile::t_pre->toc();
    log->info("Finished {} preprocessing", sim_name);
}



void class_xDMRG::run_simulation()    {
    log->info("Starting {} simulation", sim_name);
    while(true) {
        single_xDMRG_step();
        print_status_update();
        write_state();
        write_measurements();
        write_sim_status();
        write_profiling();
        copy_from_tmp();
        check_convergence();
        update_truncation_limit();     //Will update SVD threshold iff the state precision is being limited by truncation error
        update_bond_dimension_limit(); //Will update bond dimension if the state precision is being limited by bond dimension
        try_projection();
        try_chi_quench();
        try_damping();
        try_perturbation();
        // It's important not to perform the last move.
        // That last state would not get optimized
        if (state->position_is_any_edge() and not state->is_perturbed() and not state->is_damped())
        {
            if (sim_status.iteration >= settings::xdmrg::max_sweeps)    {stop_reason = StopReason::MAX_ITERS; break;}
            if (sim_status.simulation_has_succeeded)                    {stop_reason = StopReason::SUCCEEDED; break;}
            if (sim_status.simulation_has_to_stop)                      {stop_reason = StopReason::SATURATED; break;}
            if (sim_status.num_resets > settings::precision::max_resets){stop_reason = StopReason::MAX_RESET; break;}
        }
//        log->info("Truncation errors: {}", state->get_truncation_errors());
//        log->info("Bond dimensions  : {}", tools::finite::measure::bond_dimensions(*state));
        log->trace("Finished step {}, iteration {}, position {}, direction {}", sim_status.step, sim_status.iteration, state->get_position(), state->get_direction());
        move_center_point();

        sim_status.iteration     = state->get_sweeps();
        sim_status.position      = state->get_position();
        sim_status.moves         = state->get_moves();
        sim_status.step++;
    }

    switch(stop_reason){
        case StopReason::MAX_ITERS : log->info("Finished {} simulation -- reason: MAX ITERS",sim_name) ;break;
        case StopReason::SUCCEEDED : log->info("Finished {} simulation -- reason: SUCCEEDED",sim_name) ;break;
        case StopReason::SATURATED : log->info("Finished {} simulation -- reason: SATURATED",sim_name) ;break;
        case StopReason::MAX_RESET : log->info("Finished {} simulation -- reason: MAX RESET",sim_name) ;break;
        default: log->info("Finished {} simulation -- reason: NONE GIVEN",sim_name);
    }
}



void class_xDMRG::single_xDMRG_step()
{
    log->debug("Starting xDMRG step {} | iteration {} | position {} | direction {}", sim_status.step, sim_status.iteration,state->get_position(), state->get_direction());

    using namespace tools::finite;
    using namespace tools::finite::opt;
    tools::common::profile::t_sim->tic();

    // Set the fastest mode by default
    opt::OptMode optMode    = opt::OptMode::VARIANCE;
    opt::OptSpace optSpace  = opt::OptSpace::DIRECT;
    opt::OptType optType    = opt::OptType::CPLX;


    // Setup normal conditions
//    if(state->get_chi_lim() <= 32){
//        optMode  = OptMode::VARIANCE;
//        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
//    }
    if(state->get_chi_lim() <= 12){
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

//    if(state->get_chi_lim() <= 8){
//        optMode  = OptMode::VARIANCE;
//        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
//    }

    if(state->get_chi_lim() <= 6 or sim_status.iteration < 2){
        optMode  = OptMode::OVERLAP;
        optSpace = OptSpace::SUBSPACE_ONLY;
    }


    //Setup strong overrides to normal conditions, e.g., for experiments like chi quench

//    if(chi_quench_steps > 0){
//        optMode  = OptMode::VARIANCE;
//        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
//    }


    if(sim_status.variance_mpo_has_converged){
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::DIRECT;
    }

    if(state->isReal())    optType  = OptType::REAL;


    size_t threshold = 0;
    switch(optSpace){
        case  OptSpace::DIRECT   : threshold = settings::precision::max_size_direct  ; break;
        case  OptSpace::SUBSPACE_ONLY: threshold = settings::precision::max_size_part_diag; break;
        case  OptSpace::SUBSPACE_AND_DIRECT: threshold = settings::precision::max_size_part_diag; break;
    }

    std::list<size_t> max_num_sites_list;
    // Generate a list of maximum number of active sites to try
    if(chi_quench_steps > 0)
        max_num_sites_list = {settings::precision::max_sites_multidmrg};
    else if(sim_status.simulation_has_got_stuck)
        max_num_sites_list = {settings::precision::max_sites_multidmrg};
    else
        max_num_sites_list = {2};

    // Make sure not to use SUBSPACE if the 2-site problem size is huge
    // When this happens, we can't use SUBSPACE even with two sites,
    // so we might as well give it to DIRECT, which handles larger problems.
    if(state->size_2site()  > threshold) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::DIRECT;
    }


    max_num_sites_list.sort();
    max_num_sites_list.unique();
    max_num_sites_list.remove_if([](auto &elem){return elem > settings::precision::max_sites_multidmrg;});
    if(max_num_sites_list.empty()) throw std::runtime_error("No sites selected for multisite xDMRG");

    log->debug("Possible multisite step sizes: {}", max_num_sites_list);
    size_t theta_count = 0;
    double variance_old = measure::energy_variance_per_site(*state);
    std::map<double,std::tuple<Eigen::Tensor<Scalar,3>, std::list<size_t>, size_t>> results;
    for (auto & max_num_sites : max_num_sites_list){
        if(optMode == opt::OptMode::OVERLAP and optSpace == opt::OptSpace::DIRECT)
            throw std::logic_error("OVERLAP mode and DIRECT space are incompatible");

        auto old_num_sites = state->active_sites.size();
        auto old_prob_size = state->active_problem_size();
        state->activate_sites(threshold, max_num_sites);
//        if(chi_quench_steps == 0) state->activate_sites(threshold, max_num_sites);
//        else state->activate_truncated_sites(threshold,chi_lim_quench_ahead, max_num_sites);
//         Reduce bond dimensions for some sites ahead
//        if(chi_quench_steps > 0) tools::finite::mps::truncate_active_sites(*state, chi_lim_quench_ahead);

        //Check that we are not about to solve the same problem again
        if(not results.empty() and
           state->active_sites.size()   == old_num_sites and
           state->active_problem_size() == old_prob_size) {
            // If we reached this point we have exhausted the number of sites available
            log->debug("Can't activate more sites");
            break;
        }

        // If the problem is small enough we can use a safer OVERLAP step
//        if(state->active_problem_size() <= settings::precision::max_size_part_diag){
//            optMode  = OptMode::VARIANCE;
//            optSpace = OptSpace::SUBSPACE_ONLY;
//        }

//        Eigen::Tensor<Scalar,3> theta;
//        if(chi_quench_steps > 0)
//            theta = tools::finite::opt::internal::ham_sq_optimization(*state,optType, optMode,optSpace,"SR");
//        else
        auto theta = opt::find_excited_state(*state, sim_status, optMode, optSpace,optType);
        double variance_new = measure::energy_variance_per_site(*state,theta);
        results.insert({variance_new,{theta,state->active_sites,theta_count++}});

//        if(chi_quench_steps > 0){
//            double energy_new = measure::energy_per_site(*state,theta);
//            Eigen::Tensor<Scalar,3> theta_random(state->active_dimensions());
//            Eigen::Map<Eigen::VectorXcd> theta_random_map(theta_random.data(),theta_random.size());
//            Eigen::Map<Eigen::VectorXcd> theta_map(theta.data(),theta.size());
//            for (int trial = 0; trial < 10; trial++){
//                theta_random.setRandom();
//                theta_random_map = (theta_map + variance_new * theta_random_map.real()).normalized();
//                theta_random = opt::internal::ceres_direct_optimization(*state,theta_random, sim_status, optType, OptMode::VARIANCE, OptSpace::DIRECT);
//                double variance_random = measure::energy_variance_per_site(*state,theta_random);
//                double energy_random   = measure::energy_per_site(*state,theta_random);
//                double energy_diff = std::abs(energy_random - energy_new);
//                double overlap = std::real(theta_map.dot(theta_random_map));
//                if(energy_diff < std::sqrt(variance_new) and overlap > 0.9)
//                    results.insert({variance_random,{theta_random,state->active_sites,theta_count++}});
//            }
//        }



        // We can now decide if we are happy with the result or not.
        if (variance_new < variance_old) {
            log->debug("State improved during {} optimization",optSpace);
            break;
        }else{
            log->debug("State did not improve during {} optimization", optSpace);
            continue;
        }
    }
    int result_count = 0;
    for (auto & result : results){
        tools::log->info("Result {:3} candidate {:3} variance {:.16f}", result_count++, std::get<2>(result.second),std::log10(result.first));
    }
    //Check the contents of results.
//    auto[variance_new,theta] = std::make_pair(results.begin()->first,results.begin()->second);
    state->clear_cache();
    state->clear_measurements();
    auto variance_new = results.begin()->first;
    const auto &theta = std::get<0>(results.begin()->second);
    state->active_sites = std::get<1>(results.begin()->second);

    if(std::log10(variance_new) < std::log10(variance_old) - 1e-2)
        state->tag_active_sites_have_been_updated(true);

    //Truncate theta down to chi_lim
    size_t chi_lim = state->get_chi_lim();

    //Truncate even more if doing chi quench
//    if(chi_quench_steps > 0) chi_lim = chi_lim_quench_trail;

    //Do the truncation with SVD
    auto variance_before_svd = tools::finite::measure::energy_variance_per_site(*state,theta);
    log->trace("Variance check before truncate  : {:.16f}", std::log10(variance_before_svd));
    opt::truncate_theta(theta, *state,chi_lim);
    auto variance_after_svd = tools::finite::measure::energy_variance_per_site(*state);
    state->set_truncated_variance( (variance_after_svd - variance_before_svd) / variance_after_svd );
    log->trace("Variance check after truncate   : {:.16f}", std::log10(variance_after_svd));
    log->debug("Variance loss {:.16f}", state->get_truncated_variance() );

    //Normalize if unity was lost for some reason (numerical error buildup)
    if(std::abs(tools::finite::measure::norm(*state) - 1.0) > settings::precision::max_norm_error){
        tools::log->warn("Norm too large: {:.18f}",tools::finite::measure::norm(*state) );
        tools::finite::mps::normalize(*state);
        tools::finite::mps::rebuild_environments(*state);
    }
    if(settings::precision::use_reduced_energy and state->position_is_any_edge()){
        tools::finite::mpo::reduce_mpo_energy(*state);
    }

    debug::check_integrity(*state);
    sim_status.energy_dens        = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);



    tools::common::profile::t_sim->toc();
    sim_status.wall_time = tools::common::profile::t_tot->get_age();
    sim_status.simu_time = tools::common::profile::t_sim->get_measured_time();

}



void class_xDMRG::check_convergence(){
    tools::common::profile::t_con->tic();
    if(state->position_is_any_edge()){
        check_convergence_variance();
        check_convergence_entg_entropy();
    }

    sim_status.energy_dens = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);
    bool outside_of_window = std::abs(sim_status.energy_dens - sim_status.energy_dens_target)  > sim_status.energy_dens_window;
    if (sim_status.iteration > 2 and state->position_is_any_edge())
    {
        if (    outside_of_window
            and (sim_status.variance_mpo_has_saturated or
                 sim_status.variance_mpo_has_converged or
                 tools::finite::measure::energy_variance_per_site(*state) < 1e-4))
        {
            double old_energy_dens_window = sim_status.energy_dens_window;
            double new_energy_dens_window = std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5);
            std::string reason = fmt::format("saturated outside of energy window {} ± {}", sim_status.energy_dens_target,sim_status.energy_dens_window);
            log->info("Increasing energy window: {} --> {}",old_energy_dens_window, new_energy_dens_window);
            sim_status.energy_dens_window = new_energy_dens_window;
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

    sim_status.simulation_has_saturated = ((sim_status.variance_mpo_saturated_for >= min_saturation_iters and
                                           sim_status.entanglement_saturated_for >= min_saturation_iters) or
                                           (state->get_sweeps() > settings::xdmrg::min_sweeps and not state->any_sites_updated()));


    sim_status.simulation_has_succeeded = sim_status.simulation_has_converged and
                                          sim_status.simulation_has_saturated;


    sim_status.simulation_has_got_stuck = sim_status.simulation_has_saturated and not
                                          sim_status.simulation_has_succeeded;


    if(state->position_is_any_edge()) {
        sim_status.simulation_has_stuck_for = sim_status.simulation_has_got_stuck ? sim_status.simulation_has_stuck_for + 1 : 0;
    }







    sim_status.simulation_has_to_stop = sim_status.simulation_has_stuck_for >= max_stuck_iters;

    if(state->position_is_any_edge()) {
        log->debug("Simulation has converged: {}", sim_status.simulation_has_converged);
        log->debug("Simulation has saturated: {}", sim_status.simulation_has_saturated);
        log->debug("Simulation has succeeded: {}", sim_status.simulation_has_succeeded);
        log->debug("Simulation has got stuck: {}", sim_status.simulation_has_got_stuck);
        log->debug("Simulation has stuck for: {}", sim_status.simulation_has_stuck_for);
        log->debug("Simulation has to stop  : {}", sim_status.simulation_has_to_stop);
    }


//    if (    sim_status.num_resets < settings::precision::max_resets
//            and tools::finite::measure::energy_variance_per_site(*state) > 1e-10)
//    {
//        std::string reason = fmt::format("simulation has saturated with bad precision",
//                                         sim_status.energy_dens, sim_status.energy_dens_window, sim_status.energy_dens_window);
//        reset_to_random_state_in_energy_window(settings::model::initial_parity_sector, false, reason);
//    }



    tools::common::profile::t_con->toc();
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
    log->info("Searching for product state in normalized energy range: {} +- {} and sector {}", sim_status.energy_dens_target, sim_status.energy_dens_window,parity_sector);

    sim_status.num_resets++;
    if(sim_status.num_resets > settings::precision::max_resets){
        log->info("Not allowed more resets: num resets {} > max resets {}",sim_status.num_resets, settings::precision::max_resets);
        return;
    }

    int counter = 0;
    bool outside_of_window = true;


    while(true){
        reset_to_random_state(parity_sector);
        if (inflate) inflate_initial_state();
        sim_status.energy_dens = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);
        outside_of_window      = std::abs(sim_status.energy_dens - sim_status.energy_dens_target)  >= sim_status.energy_dens_window;
        if(not outside_of_window) break;
        counter++;
        if(counter >= 2000) throw std::runtime_error(fmt::format("Failed to find initial state in energy window after {}. retries: ", counter));
        if (counter % 10 == 0 and energy_window_growth_factor != 1.0) {
            double old_energy_dens_window = sim_status.energy_dens_window;
            double new_energy_dens_window = std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5);

            log->info("Can't find state in energy window.  Increasing energy window: {} --> {}",
                      old_energy_dens_window, new_energy_dens_window );
            sim_status.energy_dens_window = new_energy_dens_window;
        }
    }
    log->info("Energy initial (per site) = {:.16f} | density = {:.8f} | retries = {}", tools::finite::measure::energy_per_site(*state), sim_status.energy_dens,counter );
    clear_saturation_status();
    sim_status.energy_ubound      = sim_status.energy_target + sim_status.energy_dens_window * (sim_status.energy_max-sim_status.energy_min);
    sim_status.energy_lbound      = sim_status.energy_target - sim_status.energy_dens_window * (sim_status.energy_max-sim_status.energy_min);

    log->info("Number of product state resets: {}",sim_status.num_resets );
}


void class_xDMRG::find_energy_range() {
    log->trace("Finding energy range");
    if (state->get_length() != num_sites()) throw std::runtime_error("find_energy_range: state lenght mismatch");
    size_t max_sweeps_during_f_range = 4;
    sim_status.iteration  = state->reset_sweeps();
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
    state->reset_moves();
    state->reset_sweeps();
}




bool   class_xDMRG::sim_on()              {return settings::xdmrg::on;}
long   class_xDMRG::chi_max()             {return settings::xdmrg::chi_max;}
size_t class_xDMRG::num_sites()           {return settings::xdmrg::num_sites;}
size_t class_xDMRG::write_freq()          {return settings::xdmrg::write_freq;}
size_t class_xDMRG::print_freq()          {return settings::xdmrg::print_freq;}
bool   class_xDMRG::chi_grow()            {return settings::xdmrg::chi_grow;}
long   class_xDMRG::chi_init()            {return settings::xdmrg::chi_init;}
bool   class_xDMRG::store_wave_function() {return settings::xdmrg::store_wavefn;}



