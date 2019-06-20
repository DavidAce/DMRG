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
        : class_algorithm_base(std::move(h5ppFile_),"fDMRG", SimulationType::fDMRG) {

    table_fdmrg       = std::make_unique<class_hdf5_table<class_table_dmrg>>        (h5pp_file, sim_name + "/measurements", "simulation_progress", sim_name);
    table_fdmrg_chain = std::make_unique<class_hdf5_table<class_table_finite_chain>>(h5pp_file, sim_name + "/measurements", "simulation_progress_full_chain", sim_name);
    mpstools::finite::chain::initialize_state(*state, settings::model::model_type, settings::model::symmetry, settings::fdmrg::num_sites);
    mpstools::finite::chain::copy_state_to_superblock(*state,*superblock);

    min_saturation_length = 0.5 * settings::fdmrg::num_sites;
    max_saturation_length = 1.0 * settings::fdmrg::num_sites;
    settings::fdmrg::min_sweeps = std::max(settings::fdmrg::min_sweeps, 1+(size_t)(std::log2(chi_max())/2));
}



void class_fDMRG::run()
/*!
 * \brief Dispatches xDMRG stages.
 * This function manages the stages of simulation differently depending on whether
 * the data already existed in hdf5 storage or not.
 *
 * There can be two main scenarios that split into cases:
 * 1) The hdf5 file existed already and contains
 *      a) nothing recognizeable (previous crash?)       -- run full simulation from scratch.
 *      b) a converged simulation but no MPS             -- run full simulation from scratch.
 *      c) a not-yet-converged MPS                       -- resume simulation, reset the number of sweeps first.
 *      d) a converged MPS                               -- not much to do... run postprocessing
 * 2) The hdf5 file did not exist                        -- run full simulation from scratch.

 *
 */
{
    if (!settings::fdmrg::on) { return; }
    t_tot.tic();

    if (h5pp_file->getCreateMode() == h5pp::CreateMode::OPEN){
        // This is case 1
        bool finOK_exists = h5pp_file->linkExists("common/finOK");
        bool simOK_exists = h5pp_file->linkExists(sim_name + "/simOK");
        bool mps_exists   = h5pp_file->linkExists(sim_name + "/state/mps");
        bool finOK = false;
        bool simOK = false;
        if(finOK_exists) finOK = h5pp_file->readDataset<bool>("common/finOK");
        if(simOK_exists) simOK = h5pp_file->readDataset<bool>(sim_name + "/simOK");



        if (not simOK or not finOK ){
            //Case 1 a -- run full simulation from scratch.
            log->trace("Case 1a");
            run_preprocessing();
            run_simulation();
        }else if(simOK and not mps_exists){
            // Case 1 b
            log->trace("Case 1b");
            run_preprocessing();
            run_simulation();
        }else if(simOK and mps_exists){
            // We can go ahead and load the state from hdf5
            log->info("Loading MPS from file");
            try{
                mpstools::finite::io::load_from_hdf5(*h5pp_file, *state, *superblock, sim_state, sim_name);
            }
            catch(std::exception &ex){
                log->error("Failed to load from hdf5: {}", ex.what());
                throw std::runtime_error("Failed to resume from file: " + std::string(ex.what()));
            }
            catch(...){log->error("Unknown error when trying to resume from file.");}

            bool convergence_was_reached  = h5pp_file->readDataset<bool>(sim_name + "/sim_state/simulation_has_converged");
            if(not convergence_was_reached){
                // Case 1 c -- resume simulation, reset the number of sweeps first.
                log->trace("Case 1c");
                settings::fdmrg::max_sweeps += state->get_sweeps();
                run_simulation();

            }else {
                // Case 1 d -- not much else to do.. redo postprocessing for good measure.
                log->trace("Case 1d");
            }
        }
    }else {
        // This is case 2
        log->trace("Case 2");
        run_preprocessing();
        run_simulation();
    }
    run_postprocessing();
    print_status_full();
    t_tot.toc();
    print_profiling();
}


void class_fDMRG::run_preprocessing() {
    log->info("Running {} preprocessing",sim_name);
    mpstools::finite::print::print_hamiltonians(*state);
    log->info("Finished {} preprocessing", sim_name);
}

void class_fDMRG::run_simulation(){
    log->info("Starting {} simulation", sim_name);
    while(true) {
        single_DMRG_step();
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
        enlarge_environment(state->get_direction());
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
void class_fDMRG::run_postprocessing(){
    log->info("Running {} postprocessing",sim_name);
    mpstools::finite::debug::check_integrity(*state,*superblock,sim_state);
    state->unset_measurements();
    state->do_all_measurements();
    mpstools::finite::io::write_all_measurements(*state, *h5pp_file, sim_name);
    mpstools::infinite::io::write_all_measurements(*superblock, *h5pp_file, sim_name);

    mpstools::finite::debug::print_parity_properties(*state);
    mpstools::finite::io::write_closest_parity_projection(*state, *h5pp_file, sim_name, settings::model::symmetry);

    //  Write the wavefunction (this is only defined for short enough chain ( L < 14 say)
    if(settings::fdmrg::store_wavefn){
        h5pp_file->writeDataset(mpstools::finite::measure::mps_wavefn(*state), sim_name + "/state/psi");
    }
    print_status_full();
    print_profiling();
    h5pp_file->writeDataset(true, sim_name + "/simOK");
    log->info("Finished {} postprocessing",sim_name);
}




void class_fDMRG::initialize_chain() {
    while(true){
        insert_superblock_to_chain();
        if (superblock->environment_size + 2ul < (unsigned long) settings::fdmrg::num_sites) {
            enlarge_environment();
            swap();
        } else {
            break;
        }
    }
}

void class_fDMRG::check_convergence(){
    t_sim.tic();
    t_con.tic();
    check_convergence_variance_mpo();
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



void class_fDMRG::store_state_and_measurements_to_file(bool force){
    if(not force){
        if (not settings::hdf5::save_progress){return;}
        if (Math::mod(sim_state.iteration, settings::fdmrg::store_freq) != 0) {return;}
        if (not state->position_is_any_edge()) {return;}
        if (settings::fdmrg::store_freq == 0){return;}
        if (settings::hdf5::storage_level <= StorageLevel::NONE){return;}
    }
    compute_observables(*state);
    log->trace("Storing all measurements to file");
    t_sto.tic();
    h5pp_file->writeDataset(false, sim_name + "/simOK");
    mpstools::finite::io::write_all_measurements(*state, *h5pp_file, sim_name);
    mpstools::finite::io::write_closest_parity_projection(*state, *h5pp_file, sim_name, settings::model::symmetry);
    //  Write the wavefunction (this is only defined for short enough chain ( L < 14 say)
    if(settings::xdmrg::store_wavefn){
        h5pp_file->writeDataset(mpstools::finite::measure::mps_wavefn(*state), sim_name + "/state/psi");
    }
    mpstools::finite::io::write_all_state(*state, *h5pp_file, sim_name);
    t_sto.toc();
    store_algorithm_state_to_file();
    h5pp_file->writeDataset(true, sim_name + "/simOK");
}


void class_fDMRG::store_table_entry_progress(bool force){
    if (not force){
        if (Math::mod(sim_state.iteration, settings::fdmrg::store_freq) != 0) {return;}
        if (not state->position_is_the_middle()) {return;}
        if (settings::fdmrg::store_freq == 0){return;}
    }

    compute_observables(*superblock);
    t_sto.tic();
    table_fdmrg->append_record(
            sim_state.iteration,
            state->get_length(),
            state->get_position(),
            superblock->measurements.bond_dimension.value(),
            settings::fdmrg::chi_max,
            superblock->measurements.energy_per_site_mpo.value(),
            superblock->measurements.energy_per_site_ham.value(),
            superblock->measurements.energy_per_site_mom.value(),
            sim_state.energy_min,
            sim_state.energy_max,
            sim_state.energy_target,
            superblock->measurements.energy_variance_per_site_mpo.value(),
            superblock->measurements.energy_variance_per_site_ham.value(),
            superblock->measurements.energy_variance_per_site_mom.value(),
            superblock->measurements.current_entanglement_entropy.value(),
            superblock->measurements.truncation_error.value(),
            t_tot.get_age());
    t_sto.toc();
}

void class_fDMRG::store_chain_entry_to_file(bool force){
    if (not force){
        if (Math::mod(sim_state.iteration, settings::fdmrg::store_freq) != 0) {return;}
        if (settings::fdmrg::store_freq == 0){return;}
        if (not (state->get_direction() == 1 or state->position_is_the_right_edge())){return;}
    }
   t_sto.tic();
    table_fdmrg_chain->append_record(
            sim_state.iteration,
            state->get_length(),
            state->get_position(),
            mpstools::common::measure::bond_dimension(*superblock),
            superblock->E_optimal / mpstools::common::measure::length(*superblock) * 2.0,
            mpstools::common::measure::current_entanglement_entropy(*superblock),
            mpstools::common::measure::truncation_error(*superblock)
        );
    t_sto.toc();
}



long   class_fDMRG::chi_max()   {return settings::fdmrg::chi_max;}
size_t class_fDMRG::num_sites() {return settings::fdmrg::num_sites;}
size_t class_fDMRG::store_freq(){return settings::fdmrg::store_freq;}
size_t class_fDMRG::print_freq(){return settings::fdmrg::print_freq;}
bool   class_fDMRG::chi_grow()  {return settings::fdmrg::chi_grow;}


void class_fDMRG::print_profiling(){
    if (settings::profiling::on) {
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_ste.print_time_w_percent(t_tot);
        t_prt.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        print_profiling_sim(t_sim);
        mpstools::common::profiling::obs::print_profiling(t_obs);
    }
}

void class_fDMRG::print_profiling_sim(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "Simulation breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_opt.print_time_w_percent(t_parent);
        t_svd.print_time_w_percent(t_parent);
        t_env.print_time_w_percent(t_parent);
        t_mps.print_time_w_percent(t_parent);
        t_con.print_time_w_percent(t_parent);
        std::cout << std::endl;
    }
}