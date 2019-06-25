//
// Created by david on 2018-01-31.
//


#include <iomanip>
#include <io/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_state/class_superblock.h>
#include <mps_state/class_finite_chain_state.h>
#include <mps_tools/nmspc_mps_tools.h>
#include <general/nmspc_math.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_random_numbers.h>
#include <h5pp/h5pp.h>
#include "class_fDMRG.h"


using namespace std;
using namespace Textra;

class_fDMRG::class_fDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
        : class_algorithm_finite(std::move(h5ppFile_),"fDMRG", SimulationType::fDMRG) {
    mpstools::finite::mpo::initialize_mpo(*state, settings::model::model_type, settings::fdmrg::num_sites);
    mpstools::finite::mps::initialize_mps(*state, settings::model::symmetry,   settings::fdmrg::num_sites);
    min_saturation_length = 0.5 * settings::fdmrg::num_sites;
    max_saturation_length = 1.0 * settings::fdmrg::num_sites;
    settings::fdmrg::min_sweeps = std::max(settings::fdmrg::min_sweeps, 1+(size_t)(std::log2(chi_max())/2));
}



void class_fDMRG::run_simulation(){
    log->info("Starting {} simulation", sim_name);
    while(true) {
        single_DMRG_step("SR");
        copy_superblock_to_chain();         //Needs to occurr after update_MPS...
        store_table_entry_progress();
        store_chain_entry_to_file();
        store_profiling_deltas();
        store_state_and_measurements_to_file();

        check_convergence();
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if (sim_state.iteration >= settings::fdmrg::min_sweeps and state->position_is_the_middle_any_direction())
        {
            if (sim_state.iteration >= settings::fdmrg::max_sweeps) {stop_reason = StopReason::MAX_STEPS; break;}
            if (sim_state.simulation_has_converged)                 {stop_reason = StopReason::CONVERGED; break;}
            if (sim_state.simulation_has_to_stop)                   {stop_reason = StopReason::SATURATED; break;}
        }
        update_bond_dimension(min_saturation_length);
        move_center_point();
        sim_state.iteration = state->get_sweeps();
        sim_state.position  = state->get_position();
        sim_state.step++;
        log->trace("Finished step {}, iteration {}",sim_state.step,sim_state.iteration);
    }
    switch(stop_reason){
        case StopReason::MAX_STEPS : log->info("Finished {} simulation -- reason: MAX_STEPS",sim_name) ;break;
        case StopReason::CONVERGED : log->info("Finished {} simulation -- reason: CONVERGED",sim_name) ;break;
        case StopReason::SATURATED : log->info("Finished {} simulation -- reason: SATURATED",sim_name) ;break;
        default: log->info("Finished {} simulation -- reason: NONE GIVEN",sim_name);
    }


}


void class_fDMRG::check_convergence(){
    t_sim.tic();
    t_con.tic();
    check_convergence_variance();
    if(state->position_is_the_middle_any_direction()){
        check_convergence_entg_entropy();
    }

    if (sim_state.iteration <= settings::xdmrg::min_sweeps){
        clear_saturation_status();
    }


    if(     sim_state.variance_mpo_has_converged
            and sim_state.entanglement_has_converged)
    {
        log->debug("Simulation has converged");
        sim_state.simulation_has_converged = true;
    }

    if (        sim_state.variance_mpo_has_saturated
                and sim_state.entanglement_has_saturated
                and sim_state.bond_dimension_has_reached_max
                and sim_state.variance_mpo_saturated_for > max_saturation_length)
    {
        log->debug("Simulation has to stop");
        sim_state.simulation_has_to_stop = true;
    }


    t_con.toc();
    t_sim.toc();

}



bool   class_fDMRG::sim_on()    {return settings::fdmrg::on;}
long   class_fDMRG::chi_max()   {return settings::fdmrg::chi_max;}
size_t class_fDMRG::num_sites() {return settings::fdmrg::num_sites;}
size_t class_fDMRG::store_freq(){return settings::fdmrg::store_freq;}
size_t class_fDMRG::print_freq(){return settings::fdmrg::print_freq;}
bool   class_fDMRG::chi_grow()  {return settings::fdmrg::chi_grow;}
bool   class_fDMRG::store_wave_function()  {return settings::fdmrg::store_wavefn;}


