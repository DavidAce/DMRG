//
// Created by david on 2019-06-24.
//

#include "class_algorithm_finite.h"
#include <h5pp/h5pp.h>
#include <mps_tools/nmspc_mps_tools.h>
#include <mps_state/class_finite_chain_state.h>
#include <general/nmspc_math.h>

class_algorithm_finite::class_algorithm_finite(std::shared_ptr<h5pp::File> h5ppFile_, std::string sim_name, SimulationType sim_type)
    : class_algorithm_base(std::move(h5ppFile_), sim_name,sim_type)

{
    log->trace("Constructing default state");
    state            = std::make_shared<class_finite_chain_state>();
    table_dmrg       = std::make_unique<class_hdf5_table<class_table_dmrg>>        (h5pp_file, sim_name + "/measurements", "simulation_progress", sim_name);
    table_dmrg_chain = std::make_unique<class_hdf5_table<class_table_finite_chain>>(h5pp_file, sim_name + "/measurements", "simulation_progress_full_chain", sim_name);


    min_saturation_length = 1;// * (int)(1.0 * settings::xdmrg::num_sites);
    max_saturation_length = 2;// * (int)(2.0 * settings::xdmrg::num_sites);

    settings::xdmrg::min_sweeps = std::max(settings::xdmrg::min_sweeps, 1+(size_t)(std::log2(chi_max())/2));


    class_xDMRG::class_xDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
    : class_algorithm_finite(std::move(h5ppFile_), "xDMRG",SimulationType::xDMRG) {

    }

}

void class_algorithm_finite::run()
/*!
 * \brief Dispatches finite DMRG stages.
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
    if (!sim_on()) { return; }
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


        if (not simOK or not finOK){
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
            log->trace("Loading MPS from file");
            try{
                mpstools::finite::io::load_from_hdf5(*h5pp_file, *state, sim_state, sim_name);
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
                settings::xdmrg::max_sweeps += state->get_sweeps();
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



void class_algorithm_finite::initialize_chain() {
    while(true){
        insert_superblock_to_chain();
        if (state->get_length()  < (unsigned long) num_sites()) {
            enlarge_environment();
            swap();
        } else {
            break;
        }
    }
}

void class_algorithm_finite::check_convergence_entg_entropy(double slope_threshold) {
    //Based on the the slope of entanglement middle_entanglement_entropy
    // This one is cheap to compute.
    log->debug("Checking convergence of entanglement");

    slope_threshold = std::isnan(slope_threshold) ? settings::precision::EntEntrSaturationThreshold  : slope_threshold;
    double entropysum = 0;
    for(auto &entropy : mpstools::finite::measure::entanglement_entropies(*state)){entropysum += entropy;}
    entropysum /= state->get_length();
    check_saturation_using_slope(BS_vec,
                                 S_vec,
                                 XS_vec,
                                 entropysum,
                                 sim_state.step,
                                 1,
                                 slope_threshold,
                                 S_slope,
                                 sim_state.entanglement_has_saturated);
    sim_state.entanglement_has_converged = sim_state.entanglement_has_saturated;
}


void class_algorithm_finite::run_postprocessing(){
    log->info("Running {} postprocessing",sim_name);
    mpstools::finite::debug::check_integrity(*state,sim_state);

//    mpstools::finite::debug::check_normalization_routine(*state);

    state->unset_measurements();
    log->info("Storing state and measurements to file");
    store_state_and_measurements_to_file(true);
    log->info("Finished {} postprocessing",sim_name);
}

void class_algorithm_finite::compute_observables(){
    log->trace("Starting all measurements on current superblock");
    t_sim.tic();
    t_obs.tic();
    state->do_all_measurements();
    t_obs.toc();
    t_sim.toc();
}


void class_algorithm_finite::store_table_entry_site_state(bool force){
    if (not force) {
        if (Math::mod(sim_state.iteration, store_freq()) != 0) { return; }
        if (store_freq() == 0) { return; }
//        if (not state->position_is_the_middle()) { return; }
        if (settings::hdf5::storage_level < StorageLevel::FULL){return;}
    }

    log->trace("Storing chain_entry to file");
    t_sto.tic();
    table_dmrg_chain->append_record(
            sim_state.iteration,
            state->get_length(),
            state->get_position(),
            mpstools::common::measure::bond_dimension(*state),
            sim_state.energy_now,
            mpstools::common::measure::current_entanglement_entropy(*superblock),
            mpstools::common::measure::truncation_error(*superblock)
    );
    t_sto.toc();
}



void class_algorithm_finite::print_profiling(){
    if (settings::profiling::on) {
        log->trace("Printing profiling information (tot)");
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

void class_algorithm_finite::print_profiling_sim(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        log->trace("Printing profiling information (sim)");
        std::cout << "\n Simulation breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_opt.print_time_w_percent(t_parent);
        t_eig.print_time_w_percent(t_parent);
        t_ham.print_time_w_percent(t_parent);
        t_svd.print_time_w_percent(t_parent);
        t_env.print_time_w_percent(t_parent);
        t_mps.print_time_w_percent(t_parent);
        t_con.print_time_w_percent(t_parent);
    }
}