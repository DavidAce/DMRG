//
// Created by david on 2018-01-31.
//


#include <iomanip>
#include <io/class_hdf5_log_buffer.h>
#include <simulation/nmspc_settings.h>
#include <state/class_infinite_state.h>
#include <state/class_finite_state.h>
#include <tools/nmspc_tools.h>
#include <math/nmspc_math.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_random_numbers.h>
#include <h5pp/h5pp.h>
#include "class_fDMRG.h"


using namespace std;
using namespace Textra;

class_fDMRG::class_fDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
        : class_algorithm_finite(std::move(h5ppFile_),"fDMRG", SimulationType::fDMRG, settings::fdmrg::num_sites) {
    log->trace("Constructing class_fDMRG");
    log_dmrg       = std::make_unique<class_hdf5_log<class_log_dmrg>>        (h5pp_file, sim_name + "/measurements", "simulation_progress", sim_name);
    settings::fdmrg::min_sweeps = std::max(settings::fdmrg::min_sweeps, (size_t)(std::log2(chi_max())));
}



void class_fDMRG::run_simulation(){
    log->info("Starting {} simulation", sim_name);
    while(true) {
        single_DMRG_step("SR");
        write_measurements();
        write_state();
        write_status();
        write_logs();
        check_convergence();
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if (sim_status.iteration >= settings::fdmrg::min_sweeps and state->position_is_the_middle_any_direction())
        {
            if (sim_status.iteration >= settings::fdmrg::max_sweeps) {stop_reason = StopReason::MAX_STEPS; break;}
            if (sim_status.simulation_has_converged)                 {stop_reason = StopReason::CONVERGED; break;}
            if (sim_status.simulation_has_to_stop)                   {stop_reason = StopReason::SATURATED; break;}
        }
        update_bond_dimension();
        sim_status.iteration = state->get_sweeps();
        sim_status.position  = state->get_position();
        sim_status.step++;
        log->trace("Finished step {}, iteration {}",sim_status.step,sim_status.iteration);
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

    if(state->position_is_any_edge()){
        check_convergence_variance();
        check_convergence_entg_entropy();
    }

    if (sim_status.iteration <= settings::xdmrg::min_sweeps){
        clear_saturation_status();
    }



    if(    sim_status.variance_mpo_has_converged
           and sim_status.entanglement_has_converged
//           and sim_status.variance_mpo_saturated_for >= min_saturation_iters
//           and sim_status.entanglement_saturated_for >= min_saturation_iters
            )
    {
        log->debug("Simulation has converged");
        sim_status.simulation_has_converged = true;
    }

    if (sim_status.bond_dimension_has_reached_max
        and (  sim_status.variance_mpo_saturated_for >= max_saturation_iters
               or sim_status.entanglement_saturated_for >= max_saturation_iters)
            )
    {
        log->debug("Simulation has to stop");
        sim_status.simulation_has_to_stop = true;
    }



    if (state->position_is_any_edge()
        and sim_status.variance_mpo_has_saturated
        and not sim_status.simulation_has_converged
        and not projected_during_saturation)
    {
        log->info("Projecting to {} due to saturation", settings::model::target_parity_sector);
        *state = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::model::target_parity_sector);
        projected_during_saturation = true;
    }


    t_con.toc();
    t_sim.toc();

}


void class_fDMRG::write_logs(bool force){
    if(not force){
        if (not settings::output::save_logs){return;}
        if (math::mod(sim_status.iteration, write_freq()) != 0) {return;}
        if (settings::output::storage_level < StorageLevel::NORMAL){return;}
    }
    log_sim_status->append_record(sim_status);
//    log_profiling->append_record();
//    log_dmrg->append_record();

}





bool   class_fDMRG::sim_on()    {return settings::fdmrg::on;}
long   class_fDMRG::chi_max()   {return settings::fdmrg::chi_max;}
size_t class_fDMRG::num_sites() {return settings::fdmrg::num_sites;}
size_t class_fDMRG::write_freq(){return settings::fdmrg::write_freq;}
size_t class_fDMRG::print_freq(){return settings::fdmrg::print_freq;}
bool   class_fDMRG::chi_grow()  {return settings::fdmrg::chi_grow;}
bool   class_fDMRG::store_wave_function()  {return settings::fdmrg::store_wavefn;}


