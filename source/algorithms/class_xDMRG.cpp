//
// Created by david on 2018-02-09.
//


#include <iomanip>
#include <io/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_state/class_superblock.h>
#include <mps_state/class_mps_2site.h>
#include <mps_state/class_finite_chain_state.h>
#include <mps_tools/nmspc_mps_tools.h>
#include <mps_tools/finite/opt.h>
#include <general/nmspc_math.h>
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
    sim_state.energy_dens_target = settings::xdmrg::energy_density_target;
    sim_state.energy_dens_window = settings::xdmrg::energy_density_window;
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

        check_convergence();
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if (state->position_is_any_edge())
        {
            if (sim_state.iteration >= settings::xdmrg::max_sweeps) {stop_reason = StopReason::MAX_STEPS; break;}
            if (sim_state.simulation_has_converged)                 {stop_reason = StopReason::CONVERGED; break;}
            if (sim_state.simulation_has_to_stop)                   {stop_reason = StopReason::SATURATED; break;}
        }

        update_bond_dimension();
        move_center_point();
        sim_state.iteration = state->get_sweeps();
        sim_state.position = state->get_position();
        log->trace("Finished step {}, iteration {}, direction {}",sim_state.step,sim_state.iteration,state->get_direction());
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
    mpstools::finite::debug::check_integrity(*state);


    t_sim.tic();
    t_opt.tic();
    log->trace("Starting single xDMRG step");
    using namespace  mpstools::finite::opt;
    using namespace  mpstools::finite::measure;

    // Table

    // Mode / Space |   FULL        PARTIAL     DIRECT
    // ---------------------------------------------------
    // OVERLAP      |   FO          FP          DV
    // VARIANCE     |   FV          FV          DV


    auto optMode  =  OptMode::OVERLAP;
    optMode  = sim_state.iteration              >= 2     ?  OptMode::VARIANCE : optMode;
    optMode  = energy_variance_per_site(*state) < 1e-6   ?  OptMode::VARIANCE : optMode;



    auto optSpace =  OptSpace::PARTIAL;
//    optSpace = sim_state.iteration              >= 2                            ? OptSpace::PARTIAL : optSpace;
//    optSpace = energy_variance_per_site(*state) <  1e-6                         ? OptSpace::DIRECT  : optSpace;
//    optSpace = sim_state.iteration              >= settings::xdmrg::min_sweeps  ? OptSpace::DIRECT  : optSpace;




    long threshold = 0;
    switch(optSpace){
        case OptSpace::FULL    : threshold = 2 * 2 * 16 * 16; break;
        case OptSpace::PARTIAL : threshold = 2 * 2 * 32 * 32; break;
        case OptSpace::DIRECT  : threshold = 2 * 2 * 64 * 64; break;
    }
    state->activate_sites(threshold);
    auto optType = state->isReal() ? OptType::REAL : OptType::CPLX;

    Eigen::Tensor<Scalar,3> theta;
    std::tie(theta, sim_state.energy_now) = mpstools::finite::opt::find_optimal_excited_state(*state,sim_state,optMode, optSpace,optType);
    sim_state.energy_dens = (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min);
    t_opt.toc();

    t_svd.tic();
    mpstools::finite::opt::truncate_theta(theta, *state, sim_state.chi_temp, settings::precision::SVDThreshold);
    t_svd.toc();

    state->unset_measurements();
    auto fullnorm = mpstools::finite::measure::norm(*state);
    if(std::abs(fullnorm - 1.0) > 1e-12) {
        throw std::runtime_error(fmt::format("Norm before rebuild of env too far from unity: {}",fullnorm));
    }
    mpstools::finite::mps::rebuild_environments(*state);
    mpstools::finite::debug::check_integrity(*state);
    state->unset_measurements();

    t_sim.toc();

    sim_state.wall_time = t_tot.get_age();
    sim_state.simu_time = t_sim.get_age();

}


void class_xDMRG::check_convergence(){

    t_sim.tic();
    t_con.tic();

    check_convergence_variance();
    check_convergence_entg_entropy();

    if (sim_state.iteration < settings::xdmrg::min_sweeps){
        clear_saturation_status();
    }


    bool outside_of_window = std::abs(sim_state.energy_dens - sim_state.energy_dens_target)  > sim_state.energy_dens_window;
    if (outside_of_window
        and (   sim_state.iteration >= 2
                or mpstools::finite::measure::energy_variance_per_site(*state) < 1e-4
                or sim_state.variance_mpo_saturated_for > min_saturation_length
                or sim_state.variance_mpo_has_converged)
        )
    {
        double growth_factor = 1.10;
        log->info("Resetting to product state -- saturated outside of energy window. Energy density: {}, Energy window: {} --> {}",sim_state.energy_dens, sim_state.energy_dens_window, std::min(growth_factor*sim_state.energy_dens_window, 0.5) );
        sim_state.energy_dens_window = std::min(growth_factor*sim_state.energy_dens_window, 0.5);
        int counter = 0;
        while(outside_of_window){
            reset_to_random_state(settings::model::symmetry);
            sim_state.energy_now  = mpstools::finite::measure::energy_per_site(*state);
            sim_state.energy_dens = (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min);
            outside_of_window = std::abs(sim_state.energy_dens - sim_state.energy_dens_target)  >= sim_state.energy_dens_window;
            counter++;
            if (counter % 10 == 0) {
                log->info("Resetting to product state -- can't find state in energy window.  Increasing energy window: {} --> {}", sim_state.energy_dens_window, std::min(growth_factor*sim_state.energy_dens_window, 0.5) );
                sim_state.energy_dens_window = std::min(growth_factor*sim_state.energy_dens_window, 0.5);
            }
        }
        log->info("Energy initial (per site) = {} | density = {} | retries = {}", sim_state.energy_now, sim_state.energy_dens,counter );
        clear_saturation_status();
        projected_during_saturation  = false;
        sim_state.energy_ubound      = sim_state.energy_target + sim_state.energy_dens_window * (sim_state.energy_max-sim_state.energy_min);
        sim_state.energy_lbound      = sim_state.energy_target - sim_state.energy_dens_window * (sim_state.energy_max-sim_state.energy_min);
    }





    if(    sim_state.variance_mpo_has_converged
       and sim_state.entanglement_has_converged)
    {
        log->debug("Simulation has converged");
        sim_state.simulation_has_converged = true;
    }

    if (    sim_state.variance_mpo_has_saturated
        and sim_state.entanglement_has_saturated
        and sim_state.bond_dimension_has_reached_max
        and sim_state.variance_mpo_saturated_for > max_saturation_length)
    {
        log->debug("Simulation has to stop");
        sim_state.simulation_has_to_stop = true;
    }




    if (state->position_is_any_edge()
        and sim_state.variance_mpo_has_saturated
        and not sim_state.simulation_has_converged
        and not outside_of_window
        and not projected_during_saturation)
    {
        log->info("Projecting to {} due to saturation", settings::model::symmetry);
        *state = mpstools::finite::ops::get_closest_parity_state(*state,settings::model::symmetry);
        projected_during_saturation = true;
    }



    t_con.toc();
    t_sim.toc();
}




void class_xDMRG::find_energy_range() {
    log->trace("Finding energy range");
    assert(state->get_length() == num_sites());
    size_t max_sweeps_during_f_range = 4;
    sim_state.iteration = state->reset_sweeps();

    // Find energy minimum
    while(true) {
        class_algorithm_finite::single_DMRG_step("SR");
        print_status_update();
        // It's important not to perform the last step.
        // That last state would not get optimized
        if(state->position_is_any_edge()){
            if(sim_state.iteration >= max_sweeps_during_f_range
               or state->measurements.energy_variance_per_site.value() < 1e-8)
            {break;}
        }
        move_center_point();
        sim_state.iteration = state->get_sweeps();

    }
    compute_observables();
    sim_state.energy_min = state->measurements.energy_per_site.value();

    reset_to_random_state("sx");
    // Find energy maximum
    while(true) {
        class_algorithm_finite::single_DMRG_step("LR");
        print_status_update();
        // It's important not to perform the last step.
        // That last state would not get optimized
        if(state->position_is_any_edge()){
            if(sim_state.iteration >= max_sweeps_during_f_range
               or state->measurements.energy_variance_per_site.value() < 1e-8)
            {break;}
        }

        move_center_point();
        sim_state.iteration = state->get_sweeps();
    }
    compute_observables();
    sim_state.energy_max         = state->measurements.energy_per_site.value();
    sim_state.energy_now         = state->measurements.energy_per_site.value();
    sim_state.energy_target      = sim_state.energy_min    + sim_state.energy_dens_target  * (sim_state.energy_max - sim_state.energy_min);
    sim_state.energy_ubound      = sim_state.energy_target + sim_state.energy_dens_window  * (sim_state.energy_max - sim_state.energy_min);
    sim_state.energy_lbound      = sim_state.energy_target - sim_state.energy_dens_window  * (sim_state.energy_max - sim_state.energy_min);
    sim_state.energy_dens        = (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min);
    log->info("Energy minimum (per site) = {}", sim_state.energy_min);
    log->info("Energy maximum (per site) = {}", sim_state.energy_max);
    log->info("Energy target  (per site) = {}", sim_state.energy_target);
    int counterA = 0;
    int counterB = 0;
    bool outside_of_window = std::abs(sim_state.energy_dens - sim_state.energy_dens_target)  >= sim_state.energy_dens_window;
    while(outside_of_window){
        reset_to_random_state("sx");
        sim_state.energy_now  = mpstools::finite::measure::energy_per_site(*state);
        sim_state.energy_dens = (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min);
        outside_of_window = std::abs(sim_state.energy_dens - sim_state.energy_dens_target)  >= sim_state.energy_dens_window;

        counterA++;
        counterB++;
        if(counterA >= 100){
            counterA = 0;
            if(sim_state.energy_dens_window >= 0.5){break;}
            sim_state.energy_dens_window = std::min(1.2*sim_state.energy_dens_window, 0.5);
            sim_state.energy_ubound       = sim_state.energy_target +  sim_state.energy_dens_window*(sim_state.energy_max-sim_state.energy_min);
            sim_state.energy_lbound       = sim_state.energy_target -  sim_state.energy_dens_window*(sim_state.energy_max-sim_state.energy_min);
        }
    }
    log->info("Energy initial (per site) = {} | density = {} | retries = {}", sim_state.energy_now, sim_state.energy_dens,counterB );
}









bool   class_xDMRG::sim_on()    {return settings::xdmrg::on;}
long   class_xDMRG::chi_max()   {return settings::xdmrg::chi_max;}
size_t class_xDMRG::num_sites() {return settings::xdmrg::num_sites;}
size_t class_xDMRG::store_freq(){return settings::xdmrg::store_freq;}
size_t class_xDMRG::print_freq(){return settings::xdmrg::print_freq;}
bool   class_xDMRG::chi_grow()  {return settings::xdmrg::chi_grow;}
bool   class_xDMRG::store_wave_function()  {return settings::fdmrg::store_wavefn;}



