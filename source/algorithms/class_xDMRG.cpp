//
// Created by david on 2018-02-09.
//


#include <iomanip>
#include <io/class_hdf5_table_buffer2.h>
#include <simulation/nmspc_settings.h>
#include <state/class_infinite_state.h>
#include <state/class_mps_2site.h>
#include <state/class_finite_state.h>
#include <state/tools/nmspc_tools.h>
#include <state/tools/finite/opt.h>
#include <math/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <general/nmspc_quantum_mechanics.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <general/nmspc_tensor_extra.h>
#include <h5pp/h5pp.h>
#include <spdlog/spdlog.h>

#include "class_xDMRG.h"



using namespace std;
using namespace Textra;

class_xDMRG::class_xDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
        : class_algorithm_finite(std::move(h5ppFile_), "xDMRG",SimulationType::xDMRG, settings::xdmrg::num_sites) {
    log->trace("Constructing class_xDMRG");
    settings::xdmrg::min_sweeps = std::max(settings::xdmrg::min_sweeps, 1+(size_t)(std::log2(chi_max())/2));
}







void class_xDMRG::run_preprocessing() {

    log->info("Starting {} preprocessing", sim_name);
    sim_status.energy_dens_target = settings::xdmrg::energy_density_target;
    sim_status.energy_dens_window = settings::xdmrg::energy_density_window;
    find_energy_range();
    log->info("Finished {} preprocessing", sim_name);
}

void class_xDMRG::run_simulation()    {
    log->info("Starting {} simulation", sim_name);
    while(true) {
        single_DMRG_step();
//        store_table_entry_progress();
        store_table_entry_site_state();
        store_profiling_totals();
        store_state_and_measurements_to_file();

        print_status_update();
        check_convergence();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if (state->position_is_any_edge())
        {
            if (sim_status.iteration >= settings::xdmrg::max_sweeps) {stop_reason = StopReason::MAX_STEPS; break;}
            if (sim_status.simulation_has_converged)                 {stop_reason = StopReason::CONVERGED; break;}
            if (sim_status.simulation_has_to_stop)                   {stop_reason = StopReason::SATURATED; break;}
        }

        update_bond_dimension();
        move_center_point();
        sim_status.iteration = state->get_sweeps();
        sim_status.position = state->get_position();
        log->trace("Finished step {}, iteration {}, direction {}",sim_status.step,sim_status.iteration,state->get_direction());
    }
    switch(stop_reason){
        case StopReason::MAX_STEPS : log->info("Finished {} simulation -- reason: MAX_ITERS",sim_name) ;break;
        case StopReason::CONVERGED : log->info("Finished {} simulation -- reason: CONVERGED",sim_name) ;break;
        case StopReason::SATURATED : log->info("Finished {} simulation -- reason: SATURATED",sim_name) ;break;
        default: log->info("Finished {} simulation -- reason: NONE GIVEN",sim_name);
    }
}




void class_xDMRG::single_DMRG_step()
/*!
 * \fn void single_DMRG_step()
 */
{
    using namespace tools::finite;

    t_sim.tic();
    log->trace("Starting single xDMRG step {}", sim_status.step);
    tools::log->debug("Variance accurate check before xDMRG step: {:.16f}", std::log10(measure::accurate::energy_variance_per_site(*state)));

    auto optMode  = sim_status.iteration  >= 2     ?   opt::OptMode::VARIANCE :  opt::OptMode::OVERLAP;
    auto optSpace = opt::OptSpace::SUBSPACE;
    optSpace =  tools::finite::measure::energy_variance_per_site(*state) <  1e-8 ?  opt::OptSpace::DIRECT  : optSpace;
    optSpace = sim_status.iteration              >= settings::xdmrg::min_sweeps  ?  opt::OptSpace::DIRECT  : optSpace;
    long threshold = 0;
    switch(optSpace){
        case  opt::OptSpace::SUBSPACE : threshold = settings::precision::MaxSizePartDiag; break;
        case  opt::OptSpace::DIRECT   : threshold = 2 * 2 * 64 * 128; break;
    }
    state->activate_sites(threshold);
    if (state->active_size() > settings::precision::MaxSizePartDiag) optSpace =  opt::OptSpace::DIRECT;

    auto optType = state->isReal() ?  opt::OptType::REAL :  opt::OptType::CPLX;

    Eigen::Tensor<Scalar,3> theta = opt::find_excited_state(*state, sim_status, optMode, optSpace,optType);
    opt::truncate_theta(theta, *state, sim_status.chi_temp, settings::precision::SVDThreshold);
    mps::rebuild_environments(*state);
    debug::check_integrity(*state);
    tools::log->debug("Variance accurate check after  xDMRG step: {:.16f}", std::log10(measure::accurate::energy_variance_per_site(*state)));

    state->unset_measurements();

    t_sim.toc();
    sim_status.wall_time = t_tot.get_age();
    sim_status.simu_time = t_sim.get_age();

}


void class_xDMRG::check_convergence(){

    t_sim.tic();
    t_con.tic();

    check_convergence_variance();
    check_convergence_entg_entropy();

    if (sim_status.iteration < settings::xdmrg::min_sweeps){
        clear_saturation_status();
    }


    sim_status.energy_dens = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);
    bool outside_of_window = std::abs(sim_status.energy_dens - sim_status.energy_dens_target)  > sim_status.energy_dens_window;
    if (outside_of_window
        and (   sim_status.iteration >= 2
                or tools::finite::measure::energy_variance_per_site(*state) < 1e-4
                or sim_status.variance_mpo_saturated_for > min_saturation_length
                or sim_status.variance_mpo_has_converged)
        )
    {
        double growth_factor = 1.10;
        log->info("Resetting to product state -- saturated outside of energy window. Energy density: {}, Energy window: {} --> {}",sim_status.energy_dens, sim_status.energy_dens_window, std::min(growth_factor*sim_status.energy_dens_window, 0.5) );
        sim_status.energy_dens_window = std::min(growth_factor*sim_status.energy_dens_window, 0.5);
        int counter = 0;
        while(outside_of_window){
            reset_to_random_state(settings::model::symmetry);
            sim_status.energy_dens = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);
            outside_of_window = std::abs(sim_status.energy_dens - sim_status.energy_dens_target)  >= sim_status.energy_dens_window;
            counter++;
            if (counter % 10 == 0) {
                log->info("Resetting to product state -- can't find state in energy window.  Increasing energy window: {} --> {}", sim_status.energy_dens_window, std::min(growth_factor*sim_status.energy_dens_window, 0.5) );
                sim_status.energy_dens_window = std::min(growth_factor*sim_status.energy_dens_window, 0.5);
            }
        }
        log->info("Energy initial (per site) = {} | density = {} | retries = {}", tools::finite::measure::energy_per_site(*state), sim_status.energy_dens,counter );
        clear_saturation_status();
        projected_during_saturation  = false;
        sim_status.energy_ubound      = sim_status.energy_target + sim_status.energy_dens_window * (sim_status.energy_max-sim_status.energy_min);
        sim_status.energy_lbound      = sim_status.energy_target - sim_status.energy_dens_window * (sim_status.energy_max-sim_status.energy_min);
    }


    if(    sim_status.variance_mpo_has_converged
       and sim_status.entanglement_has_converged)
    {
        log->debug("Simulation has converged");
        sim_status.simulation_has_converged = true;
    }

    if (    sim_status.variance_mpo_has_saturated
        and sim_status.entanglement_has_saturated
        and sim_status.bond_dimension_has_reached_max
        and sim_status.variance_mpo_saturated_for > max_saturation_length)
    {
        log->debug("Simulation has to stop");
        sim_status.simulation_has_to_stop = true;
    }




    if (state->position_is_any_edge()
        and sim_status.variance_mpo_has_saturated
        and not sim_status.simulation_has_converged
        and not outside_of_window
        and not projected_during_saturation)
    {
        log->info("Projecting to {} due to saturation", settings::model::symmetry);
        *state = tools::finite::ops::get_closest_parity_state(*state,settings::model::symmetry);
        projected_during_saturation = true;
    }



    t_con.toc();
    t_sim.toc();
}




void class_xDMRG::find_energy_range() {
    log->trace("Finding energy range");
    if (state->get_length() != num_sites()) throw std::runtime_error("find_energy_range: state lenght mismatch");
    size_t max_sweeps_during_f_range = 4;
    sim_status.iteration = state->reset_sweeps();

    // Find energy minimum
    while(true) {
        class_algorithm_finite::single_DMRG_step("SR");
        print_status_update();
        // It's important not to perform the last step.
        // That last state would not get optimized
        if(state->position_is_any_edge()){
            if(sim_status.iteration >= max_sweeps_during_f_range
               or state->measurements.energy_variance_per_site.value() < 1e-8)
            {break;}
        }
        move_center_point();
        sim_status.iteration = state->get_sweeps();

    }
    compute_observables();
    sim_status.energy_min = state->measurements.energy_per_site.value();

    reset_to_random_state("sx");
    // Find energy maximum
    while(true) {
        class_algorithm_finite::single_DMRG_step("LR");
        print_status_update();
        // It's important not to perform the last step.
        // That last state would not get optimized
        if(state->position_is_any_edge()){
            if(sim_status.iteration >= max_sweeps_during_f_range
               or state->measurements.energy_variance_per_site.value() < 1e-8)
            {break;}
        }

        move_center_point();
        sim_status.iteration = state->get_sweeps();
    }
    sim_status.energy_max         = tools::finite::measure::energy_per_site(*state);
    sim_status.energy_target      = sim_status.energy_min    + sim_status.energy_dens_target  * (sim_status.energy_max - sim_status.energy_min);
    sim_status.energy_ubound      = sim_status.energy_target + sim_status.energy_dens_window  * (sim_status.energy_max - sim_status.energy_min);
    sim_status.energy_lbound      = sim_status.energy_target - sim_status.energy_dens_window  * (sim_status.energy_max - sim_status.energy_min);
    sim_status.energy_dens        = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);
    log->info("Energy minimum (per site) = {}", sim_status.energy_min);
    log->info("Energy maximum (per site) = {}", sim_status.energy_max);
    log->info("Energy target  (per site) = {}", sim_status.energy_target);
    int counterA = 0;
    int counterB = 0;
    bool outside_of_window = std::abs(sim_status.energy_dens - sim_status.energy_dens_target)  >= sim_status.energy_dens_window;
    while(outside_of_window){
        reset_to_random_state("sx");
        sim_status.energy_dens = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);
        outside_of_window = std::abs(sim_status.energy_dens - sim_status.energy_dens_target)  >= sim_status.energy_dens_window;

        counterA++;
        counterB++;
        if(counterA >= 100){
            counterA = 0;
            if(sim_status.energy_dens_window >= 0.5){break;}
            sim_status.energy_dens_window = std::min(1.2*sim_status.energy_dens_window, 0.5);
            sim_status.energy_ubound       = sim_status.energy_target +  sim_status.energy_dens_window*(sim_status.energy_max-sim_status.energy_min);
            sim_status.energy_lbound       = sim_status.energy_target -  sim_status.energy_dens_window*(sim_status.energy_max-sim_status.energy_min);
        }
    }
    log->info("Energy initial (per site) = {} | density = {} | retries = {}", tools::finite::measure::energy_per_site(*state), sim_status.energy_dens,counterB );
}



bool   class_xDMRG::sim_on()    {return settings::xdmrg::on;}
long   class_xDMRG::chi_max()   {return settings::xdmrg::chi_max;}
size_t class_xDMRG::num_sites() {return settings::xdmrg::num_sites;}
size_t class_xDMRG::store_freq(){return settings::xdmrg::store_freq;}
size_t class_xDMRG::print_freq(){return settings::xdmrg::print_freq;}
bool   class_xDMRG::chi_grow()  {return settings::xdmrg::chi_grow;}
bool   class_xDMRG::store_wave_function()  {return settings::fdmrg::store_wavefn;}



