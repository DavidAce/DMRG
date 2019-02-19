//
// Created by david on 2018-01-31.
//


#include <iomanip>
#include <io/class_hdf5_file.h>
#include <io/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_finite_chain_state.h>
#include <mps_routines/nmspc_mps_tools.h>
#include <general/nmspc_math.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_random_numbers.h>
#include "class_fDMRG.h"

using namespace std;
using namespace Textra;

class_fDMRG::class_fDMRG(std::shared_ptr<class_hdf5_file> hdf5_)
        : class_algorithm_base(std::move(hdf5_),"fDMRG", SimulationType::fDMRG) {
    set_logger(sim_name);

    table_fdmrg       = std::make_unique<class_hdf5_table<class_table_dmrg>>        (hdf5, sim_name + "/measurements", "simulation_progress");
    table_fdmrg_chain = std::make_unique<class_hdf5_table<class_table_finite_chain>>(hdf5, sim_name + "/measurements", "simulation_progress_full_chain");
    state             = std::make_shared<class_finite_chain_state>();
    superblock        = std::make_shared<class_superblock>(sim_type);
    MPS_Tools::Finite::Chain::initialize_state(*state,settings::model::model_type, settings::fdmrg::num_sites, settings::model::seed);
    MPS_Tools::Finite::Chain::copy_state_to_superblock(*state,*superblock);

    min_saturation_length = 1 * (int)(0.5 * settings::fdmrg::num_sites);
    max_saturation_length = 1 * (int)(1.0 * settings::fdmrg::num_sites);
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

    if (hdf5->fileMode == class_hdf5_file::FileMode::OPEN){
        // This is case 1
        bool fileOK;
        hdf5->read_dataset(fileOK, "common/fileOK");
        bool simOK = hdf5->link_exists(sim_name + "/simOK");
        bool mpsOK = hdf5->link_exists(sim_name + "/state/full/mps");
//        hdf5->print_contents_of_group(sim_name);

        if (not simOK){
            //Case 1 a -- run full simulation from scratch.
            spdlog::trace("Case 1a");
            run_preprocessing();
            run_simulation();
            run_postprocessing();
        }else if(simOK and not mpsOK){
            // Case 1 b
            spdlog::trace("Case 1b");
            run_preprocessing();
            run_simulation();
            run_postprocessing();
        }else if(simOK and mpsOK){
            // We can go ahead and load the state from hdf5
            spdlog::info("Loading MPS from file");
            try{
                MPS_Tools::Finite::Hdf5::load_from_hdf5(*state, *superblock, sim_state, *hdf5, sim_name);
            }
            catch(std::exception &ex){
                spdlog::error("Failed to load from hdf5: {}", ex.what());
                throw std::runtime_error("Failed to resume from file: " + std::string(ex.what()));
            }
            catch(...){spdlog::error("Unknown error when trying to resume from file.");}

            bool convergence_was_reached;
            hdf5->read_dataset(convergence_was_reached,sim_name + "/sim_state/simulation_has_converged");
            if(not convergence_was_reached){
                // Case 1 c -- resume simulation, reset the number of sweeps first.
                spdlog::trace("Case 1c");
                settings::fdmrg::max_sweeps += state->get_sweeps();
                run_simulation();
                run_postprocessing();

            }else {
                // Case 1 d -- not much else to do.. redo postprocessing for good measure.
                spdlog::trace("Case 1d");
                run_postprocessing();
            }
        }
    }else {
        // This is case 2
        spdlog::trace("Case 2");
        run_preprocessing();
        run_simulation();
        run_postprocessing();
    }
    t_tot.toc();
}


void class_fDMRG::run_preprocessing() {
    spdlog::info("Running {} preprocessing",sim_name);
    initialize_superblock(settings::model::initial_state);
    initialize_chain();
    set_random_fields_in_chain_mpo();
    MPS_Tools::Finite::Print::print_hamiltonians(*state);
    spdlog::info("Finished {} preprocessing", sim_name);
}

void class_fDMRG::run_simulation(){
    spdlog::info("Starting {} simulation", sim_name);
    while(true) {
        single_DMRG_step();
        copy_superblock_to_chain();         //Needs to occurr after update_MPS...
        store_progress_to_file();
        store_chain_entry_to_file();
        store_profiling_to_file_delta();
        store_state_to_file();

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
        spdlog::trace("Finished step {}, iteration {}",sim_state.step,sim_state.iteration);
    }
    switch(stop_reason){
        case StopReason::MAX_STEPS : spdlog::info("Finished {} simulation -- reason: MAX_STEPS",sim_name) ;break;
        case StopReason::CONVERGED : spdlog::info("Finished {} simulation -- reason: CONVERGED",sim_name) ;break;
        case StopReason::SATURATED : spdlog::info("Finished {} simulation -- reason: SATURATED",sim_name) ;break;
        default: spdlog::info("Finished {} simulation -- reason: NONE GIVEN",sim_name);
    }


}
void class_fDMRG::run_postprocessing(){
    spdlog::info("Running {} postprocessing",sim_name);
    MPS_Tools::Finite::Debug::check_integrity(*state,*superblock,sim_state);
    state->set_measured_false();
    state->do_all_measurements();
    MPS_Tools::Finite::Hdf5::write_all_measurements(*state,*hdf5,sim_name);
    MPS_Tools::Infinite::Hdf5::write_all_measurements(*superblock,*hdf5,sim_name);

    MPS_Tools::Finite::Ops::print_parity_properties(*state);
    MPS_Tools::Finite::Hdf5::write_all_parity_projections(*state,*superblock,*hdf5,sim_name);

    //  Write the wavefunction (this is only defined for short enough chain ( L < 14 say)
    if(settings::fdmrg::store_wavefn){
        hdf5->write_dataset(MPS_Tools::Finite::Measure::mps_wavefn(*state), sim_name + "/state/full/wavefunction");
    }
    print_status_full();
    print_profiling();
    spdlog::info("Finished {} postprocessing",sim_name);
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
    if (sim_state.iteration <= settings::fdmrg::min_sweeps){clear_saturation_status();}
    check_convergence_variance_mpo();
    if(sim_state.variance_mpo_has_converged)
    {
        spdlog::info("Simulation has converged");
        sim_state.simulation_has_converged = true;
    }

    else if (sim_state.variance_mpo_has_saturated and
             sim_state.bond_dimension_has_reached_max and
             sim_state.variance_mpo_saturated_for > max_saturation_length
            )
    {
        spdlog::info("Simulation has to stop");
        sim_state.simulation_has_to_stop = true;
    }
    t_con.toc();
    t_sim.toc();

}



void class_fDMRG::store_state_to_file(bool force){
    if(not force){
        if (Math::mod(sim_state.iteration, settings::fdmrg::store_freq) != 0) {return;}
        if (not state->position_is_the_middle_any_direction()) {return;}
        if (settings::fdmrg::store_freq == 0){return;}
    }
    spdlog::trace("Storing storing mps to file");
    t_sto.tic();

    MPS_Tools::Finite::Hdf5::write_all_state(*state, *hdf5, sim_name);
    if (settings::hdf5::resume_from_file){
        MPS_Tools::Finite::Hdf5::write_full_mps(*state,*hdf5,sim_name);
        MPS_Tools::Finite::Hdf5::write_full_mpo(*state,*hdf5,sim_name);
    }
    t_sto.toc();
    store_sim_to_file();
}


void class_fDMRG::store_progress_to_file(bool force){
    if (not force){
        if (Math::mod(sim_state.iteration, settings::fdmrg::store_freq) != 0) {return;}
        if (not state->position_is_the_middle()) {return;}
        if (settings::fdmrg::store_freq == 0){return;}
    }

    compute_observables();
    t_sto.tic();
    table_fdmrg->append_record(
            sim_state.iteration,
            state->get_length(),
            state->get_position(),
            superblock->measurements.bond_dimension,
            settings::fdmrg::chi_max,
            superblock->measurements.energy_per_site_mpo,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            superblock->measurements.energy_variance_per_site_mpo,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            superblock->measurements.current_entanglement_entropy,
            superblock->measurements.truncation_error,
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
            MPS_Tools::Common::Measure::bond_dimension(*superblock),
            superblock->E_optimal / MPS_Tools::Common::Measure::length(*superblock) * 2.0,
            MPS_Tools::Common::Measure::current_entanglement_entropy(*superblock),
            MPS_Tools::Common::Measure::truncation_error(*superblock)
        );
    t_sto.toc();
}



long   class_fDMRG::chi_max()   {return settings::fdmrg::chi_max;}
int    class_fDMRG::num_sites() {return settings::fdmrg::num_sites;}
int    class_fDMRG::store_freq(){return settings::fdmrg::store_freq;}
int    class_fDMRG::print_freq(){return settings::fdmrg::print_freq;}
bool   class_fDMRG::chi_grow()  {return settings::fdmrg::chi_grow;}


void class_fDMRG::print_profiling(){
    if (settings::profiling::on) {
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_ste.print_time_w_percent(t_tot);
        t_prt.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        std::cout << std::endl;
        print_profiling_sim(t_sim);
        superblock->print_profiling(t_obs);
        std::cout << std::endl;
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