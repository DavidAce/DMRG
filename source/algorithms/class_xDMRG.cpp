//
// Created by david on 2018-02-09.
//


#include <iomanip>
#include <simulation/nmspc_settings.h>
#include <state/class_mps_2site.h>
#include <state/class_finite_state.h>
#include <tools/nmspc_tools.h>
#include <tools/finite/opt.h>
#include <math/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <math/nmspc_math.h>
#include <h5pp/h5pp.h>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/ranges.h>
#include "class_xDMRG.h"



using namespace std;
using namespace Textra;

class_xDMRG::class_xDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
        : class_algorithm_finite(std::move(h5ppFile_), "xDMRG",SimulationType::xDMRG, settings::xdmrg::num_sites) {
    log->trace("Constructing class_xDMRG");
    settings::xdmrg::min_sweeps = std::max(settings::xdmrg::min_sweeps, (size_t)(std::log2(chi_max())));
}




void class_xDMRG::run_preprocessing() {
    log->info("Running {} preprocessing", sim_name);
    t_pre.tic();
    sim_status.energy_dens_target = settings::xdmrg::energy_density_target;
    sim_status.energy_dens_window = settings::xdmrg::energy_density_window;
    sim_status.chi_max = chi_max();
    state->set_chi_max(sim_status.chi_max);
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
        single_DMRG_step();
        write_measurements();
        write_state();
        write_status();
        write_logs();
        check_convergence();
        print_status_update();
        // It's important not to perform the last step.
        // That last state would not get optimized
        if (state->position_is_any_edge())
        {
            if (sim_status.iteration >= settings::xdmrg::max_sweeps) {stop_reason = StopReason::MAX_STEPS; break;}
            if (sim_status.simulation_has_succeeded)                 {stop_reason = StopReason::SUCCEEDED; break;}
            if (sim_status.simulation_has_to_stop)                   {stop_reason = StopReason::SATURATED; break;}
        }

        update_bond_dimension();
        log->trace("Finished step {}, iteration {}, direction {}", sim_status.step, sim_status.iteration, state->get_direction());
        sim_status.iteration     = state->get_sweeps();
        sim_status.position      = state->get_position();
        sim_status.moves         = state->get_moves();
        sim_status.step++;
    }
    switch(stop_reason){
        case StopReason::MAX_STEPS : log->info("Finished {} simulation -- reason: MAX_ITERS",sim_name) ;break;
        case StopReason::SUCCEEDED : log->info("Finished {} simulation -- reason: SUCCEEDED", sim_name) ;break;
        case StopReason::SATURATED : log->info("Finished {} simulation -- reason: SATURATED",sim_name) ;break;
        default: log->info("Finished {} simulation -- reason: NONE GIVEN",sim_name);
    }
}




void class_xDMRG::single_DMRG_step()
{
    using namespace tools::finite;

    t_run.tic();
    log->trace("Starting single xDMRG step");
//  log->debug("Variance accurate check before xDMRG step: {:.16f}", std::log10(measure::accurate::energy_variance_per_site(*state)));

//    auto optMode  = sim_status.iteration  < 2  ?  opt::OptMode::OVERLAP : opt::OptMode::VARIANCE;
//    auto optMode  = sim_status.iteration  < 2  ?  opt::OptMode::OVERLAP : opt::OptMode::VARIANCE;
    auto optMode    = measure::energy_variance_per_site(*state) > 1e-2  ?  opt::OptMode::OVERLAP : opt::OptMode::VARIANCE;

    auto optSpace = opt::OptSpace::SUBSPACE;
//    optSpace      = measure::energy_variance_per_site(*state) < settings::precision::varianceConvergenceThreshold         ?  opt::OptSpace::DIRECT  : optSpace;
    optSpace      = state->size_2site()  > settings::precision::maxSizePartDiag ? opt::OptSpace::DIRECT : optSpace;
//    optSpace      = sim_status.iteration >= settings::xdmrg::min_sweeps                                              ?  opt::OptSpace::DIRECT  : optSpace;
    auto optType  = state->isReal() ?  opt::OptType::REAL :  opt::OptType::CPLX;


    long threshold = 0;
    switch(optSpace){
        case  opt::OptSpace::SUBSPACE : threshold = settings::precision::maxSizePartDiag; break;
        case  opt::OptSpace::DIRECT   : threshold = settings::precision::maxSizeDirect  ; break;
    }

    debug::check_integrity(*state);
    Eigen::Tensor<Scalar,3> theta;
//    std::list<size_t> max_num_sites_list = math::range_list(2ul,settings::precision::maxSitesMultiDmrg,2ul);
    std::list<size_t> max_num_sites_list = {2,settings::precision::maxSitesMultiDmrg};
    if (max_num_sites_list.back() == settings::precision::maxSitesMultiDmrg) max_num_sites_list.pop_back();
    while(true){
        auto old_num_sites = state->active_sites.size();
        auto old_prob_size = state->active_problem_size();
        state->activate_sites(threshold, max_num_sites_list.front());


        if( state->active_sites.size()   == old_num_sites and
            state->active_problem_size() == old_prob_size){
            //Reached threshold
            if( optSpace == opt::OptSpace::SUBSPACE and
                old_prob_size > settings::precision::maxSizeFullDiag){
                //Switch to DIRECT
                optSpace  = opt::OptSpace::DIRECT;
                threshold = settings::precision::maxSizeDirect;
                log->debug("SUBSPACE threshold reached, switching to DIRECT mode");
                state->activate_sites(threshold, max_num_sites_list.front());
            }
            else{
                log->debug("Keeping last theta: Reached DIRECT threshold or already did Full Diag.");
                if(theta.size() == 0) throw std::logic_error("Theta is empty!");
                break;
            }
        }

        theta = opt::find_excited_state(*state, sim_status, optMode, optSpace,optType);
        max_num_sites_list.pop_front();
        if(state->active_sites_updated()){
            log->debug("Sites successfully updated");
            break;
        }
        if(max_num_sites_list.empty()){
            log->debug("Keeping last theta: failed to find better theta and maxSitesMultiDmrg reached");
            if(theta.size() == 0) throw std::logic_error("Theta is empty!");
            break;
        }
        if(state->get_direction() == 1  and state->get_position() - 1 + max_num_sites_list.front() >= state->get_length()){
            log->debug("Keeping last theta: can't activate more sites, reached the right edge");
            break;
        }
        if(state->get_direction() == -1 and state->get_position() + 2 - max_num_sites_list.front() < 0){
            log->debug("Keeping last theta: can't activate more sites, reached the left edge");
            break;
        }

    }


//    if (optMode == opt::OptMode::OVERLAP){
//        sim_status.chi_temp = 4 * (1+sim_status.iteration);
//        log->debug("Forcing bond dimension down to {}",sim_status.chi_temp);
//    }


    log->debug("Variance check before truncate       : {:.16f}", std::log10(measure::energy_variance_per_site(*state,theta)));

    opt::truncate_theta(theta, *state, sim_status.chi_temp, settings::precision::SVDThreshold);
    move_center_point();
    tools::finite::mps::rebuild_environments(*state);
    log->debug("Variance check after truncate + move : {:.16f}", std::log10(measure::energy_variance_per_site(*state)));

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

//    mps::rebuild_environments(*state);


    t_run.toc();
    sim_status.wall_time = t_tot.get_age();
    sim_status.simu_time = t_run.get_measured_time();

}


void class_xDMRG::check_convergence(){
    t_con.tic();
    if(state->position_is_any_edge()){
        check_convergence_variance();
        check_convergence_entg_entropy();
    }


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
        else
        if( not     state->all_sites_updated()
            and not sim_status.variance_mpo_has_converged
            and     sim_status.variance_mpo_has_saturated
            and     tools::finite::measure::energy_variance_per_site(*state) > 1e-4)
        {
            sim_status.energy_dens_window = std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5);
            std::string reason = fmt::format("could not update all sites. Energy density: {}, Energy window: {} --> {}",
                     sim_status.energy_dens, sim_status.energy_dens_window, std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5) );
            reset_to_random_state_in_energy_window(settings::model::initial_parity_sector, false, reason);
        }
    }

    sim_status.simulation_has_converged = sim_status.variance_mpo_has_converged and
                                          sim_status.entanglement_has_converged;
    sim_status.simulation_has_saturated = sim_status.variance_mpo_saturated_for >= min_saturation_iters and
                                          sim_status.entanglement_saturated_for >= min_saturation_iters;

    sim_status.simulation_has_succeeded = sim_status.simulation_has_converged and
                                          sim_status.simulation_has_saturated;



    log->debug("Simulation has converged: {}", sim_status.simulation_has_converged);
    log->debug("Simulation has saturated: {}", sim_status.simulation_has_saturated);
    log->debug("Simulation has succeeded: {}", sim_status.simulation_has_succeeded);

    if(        sim_status.bond_dimension_has_reached_max
       and     sim_status.simulation_has_saturated
       and not sim_status.simulation_has_converged)
    {
        if (        settings::model::project_when_saturated and
                not has_projected
            and not outside_of_window )
        {
            log->info("Projecting to {} due to saturation", settings::model::target_parity_sector);
            bool keep_bond_dimensions = true;
            *state = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::model::target_parity_sector,keep_bond_dimensions);
            has_projected = true;
        }
        else if (    sim_status.num_resets < settings::precision::maxResets
                 and tools::finite::measure::energy_variance_per_site(*state) > 1e-10)
        {
            std::string reason = fmt::format("simulation has saturated with bad precision",
                                             sim_status.energy_dens, sim_status.energy_dens_window, sim_status.energy_dens_window);
            reset_to_random_state_in_energy_window(settings::model::initial_parity_sector, false, reason);
        }

    }


    if (        settings::model::project_when_saturated
        and     state->position_is_any_edge()
        and     sim_status.variance_mpo_has_saturated
        and not sim_status.simulation_has_converged
        and not outside_of_window
        and not has_projected)
    {
        log->info("Projecting to {} due to saturation", settings::model::target_parity_sector);
        bool keep_bond_dimensions = true;
        *state = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::model::target_parity_sector,keep_bond_dimensions);
        has_projected = true;
    }




    sim_status.simulation_has_to_stop = sim_status.bond_dimension_has_reached_max
                                        and sim_status.simulation_has_saturated
                                        and (sim_status.variance_mpo_saturated_for >= max_saturation_iters or
                                             sim_status.entanglement_saturated_for >= max_saturation_iters);
    log->debug("Simulation has to stop: {}", sim_status.simulation_has_to_stop);

    t_con.toc();
}


void class_xDMRG::inflate_initial_state(){
    tools::log->trace("Inflating bond dimension");
    // Inflate by projecting randomly. Each projection doubles the bond dimension
    bool keep_bond_dimensions = false;
    for (int i = 0; i < 4; i++){
        *state = tools::finite::ops::get_projection_to_closest_parity_sector(*state, "random" ,keep_bond_dimensions);
        log->debug("χ = {}"         , tools::finite::measure::bond_dimensions(*state));
    }
    *state = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::model::initial_parity_sector ,keep_bond_dimensions);
}


void class_xDMRG::reset_to_random_state_in_energy_window(const std::string &parity_sector,bool inflate, std::string reason ){
    log->warn("Resetting to product state -- Reason: {}", reason);
//    sim_status.energy_dens_window = std::min(growth_factor*sim_status.energy_dens_window, 0.5);

    int counter = 0;
    bool outside_of_window = true;
    while(outside_of_window){
        reset_to_random_state(parity_sector);
        if (inflate) inflate_initial_state();

        sim_status.energy_dens = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);
        outside_of_window = std::abs(sim_status.energy_dens - sim_status.energy_dens_target)  >= sim_status.energy_dens_window;
        counter++;
        if (counter % 10 == 0) {
            log->info("Resetting to product state -- can't find state in energy window.  Increasing energy window: {} --> {}",
                    sim_status.energy_dens_window, std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5) );
            sim_status.energy_dens_window = std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5);
        }
    }
    log->info("Energy initial (per site) = {} | density = {} | retries = {}", tools::finite::measure::energy_per_site(*state), sim_status.energy_dens,counter );
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
    // Find energy minimum
    while(true) {
        class_algorithm_finite::single_DMRG_step("SR");
        print_status_update();
        // It's important not to perform the last moves.
        // That last state would not get optimized
        if(state->position_is_any_edge()){
            if(sim_status.iteration >= max_sweeps_during_f_range
               or state->measurements.energy_variance_per_site.value() < 1e-8)
            {break;}
        }
//        move_center_point();
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
               or state->measurements.energy_variance_per_site.value() < 1e-8)
            {break;}
        }

//        move_center_point();
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




//    int counterA = 0;
//    int counterB = 0;
//    tools::finite::mps::internals::seed_state_unused = true;
//    bool outside_of_window = true;
//    while(outside_of_window){
//        reset_to_random_state(settings::model::initial_parity_sector);
//        sim_status.energy_dens = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);
//        outside_of_window = std::abs(sim_status.energy_dens - sim_status.energy_dens_target)  >= sim_status.energy_dens_window;
//        counterA++;
//        counterB++;
//        if(counterA >= 100){
//            counterA = 0;
//            if(sim_status.energy_dens_window >= 0.5){break;}
//            sim_status.energy_dens_window = std::min(1.2*sim_status.energy_dens_window, 0.5);
//            sim_status.energy_ubound       = sim_status.energy_target +  sim_status.energy_dens_window*(sim_status.energy_max-sim_status.energy_min);
//            sim_status.energy_lbound       = sim_status.energy_target -  sim_status.energy_dens_window*(sim_status.energy_max-sim_status.energy_min);
//        }
//    }
//    log->info("Energy initial (per site) = {} | density = {} | retries = {}", tools::finite::measure::energy_per_site(*state), sim_status.energy_dens,counterB );
}




bool   class_xDMRG::sim_on()    {return settings::xdmrg::on;}
long   class_xDMRG::chi_max()   {return settings::xdmrg::chi_max;}
size_t class_xDMRG::num_sites() {return settings::xdmrg::num_sites;}
size_t class_xDMRG::write_freq(){return settings::xdmrg::write_freq;}
size_t class_xDMRG::print_freq(){return settings::xdmrg::print_freq;}
bool   class_xDMRG::chi_grow()  {return settings::xdmrg::chi_grow;}
bool   class_xDMRG::store_wave_function()  {return settings::fdmrg::store_wavefn;}



