//
// Created by david on 2019-06-24.
//

#include "class_algorithm_finite.h"
#include <h5pp/h5pp.h>
#include <math/nmspc_math.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_finite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite.h>

class_algorithm_finite::class_algorithm_finite(std::shared_ptr<h5pp::File> h5ppFile_, std::string sim_name, SimulationType sim_type)
    : class_algorithm_base(std::move(h5ppFile_), sim_name, sim_type)

{
    tools::log->trace("Constructing class_algorithm_finite");
    state = std::make_unique<class_state_finite>();
    state->set_chi_lim(2); // Can't call chi_init() <-- it's a pure virtual function
    if(state->hasNaN()) throw std::runtime_error("State has NAN's before initializing it");

    tools::finite::mps::initialize(*state, settings::model::model_type, settings::model::model_size, 0);
    tools::finite::mpo::initialize(*state, settings::model::model_type, settings::model::model_size, 0);
    tools::finite::mpo::randomize(*state);
    tools::finite::mps::random_product_state(*state, settings::strategy::initial_parity_sector, settings::input::bitfield,
                                             settings::strategy::use_pauli_eigvecs);
    tools::finite::debug::check_integrity(*state);
    S_mat.resize(state->get_length() + 1);
    X_mat.resize(state->get_length() + 1);

    tools::finite::print::print_hamiltonians(*state);
    tools::finite::print::print_state(*state);
}

/*!
 * \brief Dispatches finite DMRG stages.
 * This function manages the stages of simulation differently depending on whether
 * the data already existed in hdf5 storage or not.
 *
 * We start by searching the HDF5 group tree looking for
 *
 * These are the scenarios:
 * 1) The hdf5 file existed already and contains
 *      a) nothing recognizable (previous crash?)        -- delete the file -> run full simulation from scratch.
 *      b) a converged simulation but no MPS             -- run full simulation from scratch.
 *      c) a not-yet-converged MPS                       -- resume simulation, reset the number of sweeps first.
 *      d) a converged MPS                               -- not much to do... run postprocessing
 * 2) The hdf5 file did not exist                        -- run full simulation from scratch.

        D1: output file exists: check option in settings::output::create_mode:
            - case OPEN:        Load config and simulation state -> continue simulation
            - case RENAME:      Rename output file -> go to D2
            - case TRUNCATE:    Delete output file -> go to D2

        D2: output file does not exist: Check if coming from A2
            - Yes: load simulation state from given hdf5 file -> continue simulation
            - No: start new simulation

 */
void class_algorithm_finite::run() {
    tools::log->info("Starting {}", sim_name);
    tools::common::profile::t_tot->tic();
    if(settings::output::file_collision_policy == FileCollisionPolicy::RESUME and h5pp_file->linkExists("common/storage_level")) {
        // We may want to resume this simulation.
        // Resume can imply many things
        // 1) Resume a simulation which terminated prematurely
        // 2) Resume a previously successful simulation. This may be desireable if the config
        //    wants something that is not present in the file.
        //      a) A certain number of states
        //      b) A state inside of a particular energy window
        //      c) The ground or "roof" states
        // To guide the behavior, we check the setting ResumePolicy.

        try {
            auto state_prefix = tools::common::io::h5resume::find_resumable_state(*h5pp_file, sim_name);
            if(state_prefix.empty()) throw std::runtime_error("Could not resume: no valid resume candidates found");
            tools::log->info("Resuming state [{}]", state_prefix);
            tools::finite::io::h5resume::load_all(*h5pp_file, state_prefix, *state,sim_status);
            // Now we decide what to do
            // Let's consider case 1

//            if(sim_status.simulation_has_succeeded){
//
//            }
            exit(0);
            while(h5pp_file->linkExists(sim_name + state_name)) state_name = "state_" + std::to_string(sim_status.state_number++);
        } catch(std::exception &ex) {
            tools::log->info("Could not resume state from file [{}]: {}", h5pp_file->getFilePath(), ex.what());
            exit(0);
            run_preprocessing();
        }
    } else {
        run_preprocessing();
    }

    run_simulation();
    run_postprocessing();
    tools::common::profile::t_tot->toc();

    //
    //    if(not sim_on()) return;
    //    if(not settings::input::h5_load_filename.empty()){
    //        h5pp::File h5pp_load(settings::input::h5_load_filename,h5pp::AccessMode::READONLY,h5pp::CreateMode::OPEN);
    //        tools::finite::io::h5restore::load_all(h5pp_load, sim_name, sim_status, *state);
    //    }

    //
    //
    //    tools::log->info("Starting {}", sim_name);
    //    tools::common::profile::t_tot->tic();
    //    if(h5pp_file) {
    //        // This is case 1
    //        bool finOK_exists = h5pp_file->linkExists("common/finished_all");
    //        bool mps_exists   = h5pp_file->linkExists(sim_name + "/state/mps");
    //        bool finOK        = false;
    //        if(finOK_exists) finOK = h5pp_file->readDataset<bool>("common/finished_all");
    //
    //        if(not finOK) {
    //            // Case 1 a -- run full simulation from scratch.
    //            tools::log->trace("Case 1a");
    //            run_preprocessing();
    //            run_simulation();
    //        } else if(not mps_exists) {
    //            // Case 1 b
    //            tools::log->trace("Case 1b");
    //            run_preprocessing();
    //            run_simulation();
    //        } else if(mps_exists) {
    //            // We can go ahead and load the state from output
    //            tools::log->trace("Loading MPS from file");
    //            try {
    //                tools::finite::io::h5resume::load_all(*h5pp_file, sim_name, sim_status, *state);
    //            } catch(std::exception &ex) {
    //                tools::log->error("Failed to load from output: {}", ex.what());
    //                throw std::runtime_error("Failed to resume from file: " + std::string(ex.what()));
    //            } catch(...) {
    //                tools::log->error("Unknown error when trying to resume from file.");
    //            }
    //
    //            bool convergence_was_reached = h5pp_file->readDataset<bool>(sim_name + "/sim_status/simulation_has_converged");
    //            if(not convergence_was_reached) {
    //                // Case 1 c -- resume simulation, reset the number of sweeps first.
    //                tools::log->trace("Case 1c");
    //                settings::xdmrg::max_sweeps += state->get_iteration();
    //                run_simulation();
    //
    //            } else {
    //                // Case 1 d -- not much else to do.. redo postprocessing for good measure.
    //                tools::log->trace("Case 1d");
    //            }
    //        }
    //    } else {
    //        // This is case 2
    //        tools::log->trace("Case 2");
    //        run_preprocessing();
    //        run_simulation();
    //    }
    //    tools::common::profile::t_tot->toc();
    //    run_postprocessing();
}

void class_algorithm_finite::run_old()
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
    if(not sim_on()) {
        return;
    }
    tools::log->info("Starting {}", sim_name);
    tools::common::profile::t_tot->tic();
    if(h5pp_file) {
        // This is case 1
        bool finOK_exists = h5pp_file->linkExists("common/finished_all");
        bool mps_exists   = h5pp_file->linkExists(sim_name + "/state/mps");
        bool finOK        = false;
        if(finOK_exists) finOK = h5pp_file->readDataset<bool>("common/finished_all");

        if(not finOK) {
            // Case 1 a -- run full simulation from scratch.
            tools::log->trace("Case 1a");
            run_preprocessing();
            run_simulation();
        } else if(not mps_exists) {
            // Case 1 b
            tools::log->trace("Case 1b");
            run_preprocessing();
            run_simulation();
        } else if(mps_exists) {
            // We can go ahead and load the state from output
            tools::log->trace("Loading MPS from file");
            try {
                tools::finite::io::h5resume::load_all(*h5pp_file, sim_name, *state,sim_status);
            } catch(std::exception &ex) {
                tools::log->error("Failed to load from output: {}", ex.what());
                throw std::runtime_error("Failed to resume from file: " + std::string(ex.what()));
            } catch(...) {
                tools::log->error("Unknown error when trying to resume from file.");
            }

            bool convergence_was_reached = h5pp_file->readDataset<bool>(sim_name + "/sim_status/simulation_has_converged");
            if(not convergence_was_reached) {
                // Case 1 c -- resume simulation, reset the number of sweeps first.
                tools::log->trace("Case 1c");
                settings::xdmrg::max_sweeps += state->get_iteration();
                run_simulation();

            } else {
                // Case 1 d -- not much else to do.. redo postprocessing for good measure.
                tools::log->trace("Case 1d");
            }
        }
    } else {
        // This is case 2
        tools::log->trace("Case 2");
        run_preprocessing();
        run_simulation();
    }
    tools::common::profile::t_tot->toc();
    run_postprocessing();
}

void class_algorithm_finite::run_preprocessing() {
    tools::log->info("Running {} preprocessing (base)", sim_name);
    tools::common::profile::t_pre->tic();
    tools::finite::io::h5table::write_model(*h5pp_file, sim_name + "/model", settings::output::storage_level_results, *state);
    state->set_chi_max(chi_max());
    sim_status.chi_max = chi_max();
    update_bond_dimension_limit(chi_init());
    tools::common::profile::t_pre->toc();
    tools::log->info("Finished {} preprocessing (base)", sim_name);
}

void class_algorithm_finite::single_DMRG_step(const std::string & ritz) {
    /*!
     * \fn void single_DMRG_step(std::string ritz)
     */
    tools::log->trace("Starting single xDMRG step with ritz: [{}]", ritz);
    tools::common::profile::t_sim->tic();
    Eigen::Tensor<Scalar, 4> theta = tools::finite::opt::find_ground_state(*state, ritz);
    tools::finite::opt::truncate_theta(theta, *state);
    state->clear_measurements();
    tools::common::profile::t_sim->toc();
    sim_status.wall_time = tools::common::profile::t_tot->get_age();
    sim_status.simu_time = tools::common::profile::t_sim->get_measured_time();
}

void class_algorithm_finite::run_postprocessing() {
    tools::log->info("Running {} postprocessing", sim_name);
    tools::common::profile::t_pos->tic();
    sim_status.simulation_has_finished = true;
    tools::finite::debug::check_integrity(*state);
    state->clear_measurements();
    write_to_file(StorageReason::RESULTS);
    if(not has_projected) write_to_file(StorageReason::PROJ_STATE);
    print_status_full();
    tools::common::profile::t_pos->toc();
    tools::log->info("Finished {} postprocessing", sim_name);
}

void class_algorithm_finite::move_center_point(std::optional<size_t> num_moves) {
    if(not num_moves.has_value()) {
        if(state->active_sites.empty())
            num_moves = 1ul;
        else if(settings::precision::move_sites_multidmrg == "one")
            num_moves = 1ul;
        else if(settings::precision::move_sites_multidmrg == "mid")
            num_moves = std::max(1ul, (state->active_sites.size()) / 2);
        else if(settings::precision::move_sites_multidmrg == "max")
            num_moves = std::max(1ul, state->active_sites.size() - 2ul);
        else {
            throw std::logic_error("Specify how many sites multisite should move! Expected one of {one,mid,max}, got [" +
                                   settings::precision::move_sites_multidmrg + "]");
        }
    }

    tools::log->trace("Moving center point {} steps in direction {}", num_moves.value(), state->get_direction());
    state->clear_cache();
    try {
        for(size_t i = 0; i < num_moves.value(); i++) {
            if(state->position_is_any_edge()) {
                state->increment_iter();
                has_projected = false;
            }
            tools::finite::mps::move_center_point(*state);
            state->increment_step();
            if(chi_quench_steps > 0) chi_quench_steps--;
        }
        tools::finite::debug::check_integrity(*state);
        state->active_sites.clear();
    } catch(std::exception &e) {
        tools::finite::print::print_state(*state);
        throw std::runtime_error("Failed to move center point: " + std::string(e.what()));
    }
    sim_status.iter      = state->get_iteration();
    sim_status.step      = state->get_step();
    sim_status.position  = state->get_position();
    sim_status.direction = state->get_direction();


}

void class_algorithm_finite::update_truncation_limit() {
    if(not state->position_is_any_edge()) return;
    // Will update SVD threshold iff the energy variance is being limited by truncation error
    size_t bond_at_lim_count = state->num_bonds_at_limit();
    if(bond_at_lim_count > 0) return; // Return because we should rather increase bond dimension than lower the svd threshold
    tools::log->info("Truncated variances: {}", state->get_truncated_variances());
    //
    //    double truncation_threshold = settings::precision::svd_threshold;
    //    size_t trunc_bond_count  = state->num_sites_truncated(truncation_threshold);
    //    if(trunc_bond_count > 0) settings::precision::svd_threshold *= 0.5;
    //    tools::log->info("Lowered SVD threshold to {}",settings::precision::svd_threshold);
}

void class_algorithm_finite::update_bond_dimension_limit(std::optional<long> tmp_bond_limit) {
    if(tmp_bond_limit.has_value()) {
        state->set_chi_lim(tmp_bond_limit.value());
        sim_status.chi_lim = tmp_bond_limit.value();
        return;
    }
    try {
        long chi_lim_now = state->get_chi_lim();
        if(chi_lim_now < chi_init()) throw std::logic_error("Chi limit should be larger than chi init");
    } catch(std::exception &error) {
        // If we reached this stage, either
        // 1) chi_lim is not initialized yet
        // 2) chi_lim is initialized, but it is smaller than the init value found in settings
        // Either way, we should set chi_lim to be chi_init, unless chi_init is larger than tmp_bond_limit
        tools::log->info("Setting initial bond dimension limit: {}", chi_init());
        state->set_chi_lim(chi_init());
        sim_status.chi_lim = chi_init();
        return;
    }
    if(not state->position_is_any_edge()) return;
    sim_status.chi_lim_has_reached_chi_max = state->get_chi_lim() >= chi_max();
    if(not sim_status.chi_lim_has_reached_chi_max) {
        if(chi_grow()) {
            sim_status.chi_lim = state->get_chi_lim();

            // If we got here we want grow the bond dimension limit progressively during the simulation
            // Only increment the bond dimension if
            // * No experiments are on-going like perturbation or damping
            // * the simulation is stuck
            // * the state is limited by bond dimension
            if(state->is_damped()) {
                tools::log->info("State is undergoing disorder damping -- cannot increase bond dimension yet");
                return;
            }
            if(state->is_perturbed()) {
                tools::log->info("State is undergoing perturbation -- cannot increase bond dimension yet");
                return;
            }
            if(not sim_status.simulation_has_got_stuck and state->get_chi_lim() >= 16) {
                tools::log->info("State is not stuck yet. Kept current limit {}", state->get_chi_lim());
                return;
            }
            //            if(sim_status.simulation_has_stuck_for <= 1 and state->get_chi_lim() >= 16){
            //                tools::log->info("State has not been stuck for long enough. Kept current limit {}", state->get_chi_lim());
            //                return;
            //            }
            if(tools::log->level() <= spdlog::level::info) {
                double truncation_threshold = 2 * settings::precision::svd_threshold;
                size_t trunc_bond_count     = state->num_sites_truncated(truncation_threshold);
                size_t bond_at_lim_count    = state->num_bonds_at_limit();
                tools::log->info("Truncation threshold  : {:<.8e}", truncation_threshold);
                tools::log->info("Truncation errors     : {}", state->get_truncation_errors());
                tools::log->info("Bond dimensions       : {}", tools::finite::measure::bond_dimensions(*state));
                tools::log->info("Truncated bond count  : {} ", trunc_bond_count);
                tools::log->info("Bonds at limit  count : {} ", bond_at_lim_count);
                tools::log->info("Entanglement entropies: {} ", tools::finite::measure::entanglement_entropies(*state));
            }
            bool state_is_bond_limited = state->is_bond_limited(5 * settings::precision::svd_threshold);
            if(not state_is_bond_limited) {
                tools::log->info("State is not limited by its bond dimension. Kept current limit {}", state->get_chi_lim());
                return;
            }

            // Write current results before updating bond dimension
            write_to_file(StorageReason::CHI_UPDATE);
            //            write_measurements(true);
            //            write_sim_status(true);
            //            write_profiling(true);
            long chi_lim_new = std::min(state->get_chi_max(), state->get_chi_lim() * 2);
            if(state->get_chi_lim() < 16) chi_lim_new = std::min(state->get_chi_max(), state->get_chi_lim() + 1);
            tools::log->info("Updating bond dimension limit {} -> {}", state->get_chi_lim(), chi_lim_new);
            state->set_chi_lim(chi_lim_new);
            clear_saturation_status();
            sim_status.chi_lim_has_reached_chi_max = state->get_chi_lim() == chi_max();
            if(sim_status.chi_lim_has_reached_chi_max and has_projected) has_projected = false;
            tools::log->info("Projecting at site {} to direction {} after updating bond dimension to χ = {} ", state->get_position(),
                             settings::strategy::target_parity_sector, chi_lim_new);
            copy_from_tmp(StorageReason::CHI_UPDATE);
            if(not state->position_is_any_edge()) throw std::runtime_error("Update bond dimension: no longer at edge!");
            if(settings::strategy::project_on_chi_update) {
                *state        = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::strategy::target_parity_sector);
                has_projected = true;
            }
            write_to_file(StorageReason::PROJ_STATE);
            if(settings::strategy::randomize_on_chi_update and state->get_chi_lim() >= 32) reset_to_random_current_state();
        } else {
            // Here the settings specify to just set the limit to maximum chi directly
            tools::log->info("Setting bond dimension limit to maximum = {}", chi_max());
            state->set_chi_lim(chi_max());
        }
    } else {
        tools::log->debug("Chi limit has reached max: {} -> Kept current bond dimension limit {}", chi_max(), state->get_chi_lim());
    }
    sim_status.chi_lim = state->get_chi_lim();
    if(state->get_chi_lim() > state->get_chi_max())
        throw std::runtime_error(fmt::format("chi_lim is larger than chi_max! {} > {}", state->get_chi_lim(), state->get_chi_max()));
}

void class_algorithm_finite::reset_to_initial_state() {
    tools::log->trace("Resetting MPS to initial product state in parity sector: {}, state number {}", settings::strategy::initial_parity_sector,
                      settings::input::bitfield, settings::strategy::use_pauli_eigvecs);
    if(state->get_length() != settings::model::model_size) throw std::range_error("System size mismatch");
    // Initialize state
    tools::finite::mps::random_product_state(*state, settings::strategy::initial_parity_sector, settings::input::bitfield,
                                             settings::strategy::use_pauli_eigvecs);
    clear_saturation_status();
    state->lowest_recorded_variance = 1;
    sim_status.iter                 = state->reset_iter();
}

void class_algorithm_finite::reset_to_random_product_state(const std::string &parity_sector) {
    tools::log->trace("Resetting MPS to random product state in parity sector: {}", parity_sector);
    if(state->get_length() != settings::model::model_size) throw std::range_error("System size mismatch");
    // Randomize state
    tools::finite::mps::random_product_state(*state, parity_sector, -1, settings::strategy::use_pauli_eigvecs);
    clear_saturation_status();
    state->lowest_recorded_variance = 1;
    sim_status.iter                 = state->reset_iter();
    sim_status.step                 = state->reset_step();
    auto spin_components            = tools::finite::measure::spin_components(*state);
    tools::log->info("Successfully reset to product state with global spin components: {}", spin_components);
}

void class_algorithm_finite::reset_to_random_current_state(std::optional<double> chi_lim) {
    if(not state->position_is_any_edge()) return;
    tools::log->info("Resetting MPS by flipping random spins on current state");
    if(state->get_length() != settings::model::model_size) throw std::range_error("System size mismatch");
    tools::log->debug("Bond dimensions: {}", tools::finite::measure::bond_dimensions(*state));
    // Randomize state
    tools::log->info("Flipping random spins");
    //    tools::finite::mps::random_current_state(*state,"x");
    tools::finite::mps::random_current_state(*state, "x", "z");

    // Truncate even more on explicit request
    if(chi_lim) {
        size_t chi_lim_parsed = chi_lim.value() < 1 ? (size_t)(chi_lim.value() * (double) state->find_largest_chi()) : (size_t) chi_lim.value();
        tools::finite::mps::truncate_all_sites(*state, chi_lim_parsed);
    }
    tools::log->debug("Bond dimensions: {}", tools::finite::measure::bond_dimensions(*state));

    clear_saturation_status();
    state->lowest_recorded_variance = 1;
    sim_status.iter                 = state->reset_iter();
    auto spin_components            = tools::finite::measure::spin_components(*state);
    tools::log->info("Successfully reset to random state based on current state. New components: {}", spin_components);
    if(not state->position_is_any_edge()) throw std::runtime_error("Update bond dimension: no longer at edge!");
}

void class_algorithm_finite::try_projection() {
    if(not state->position_is_any_edge()) return;
    if(has_projected) return;
    if(settings::strategy::project_on_every_sweep or (settings::strategy::project_when_stuck and sim_status.simulation_has_got_stuck)) {
        tools::log->info("Trying projection to {}", settings::strategy::target_parity_sector);
        *state        = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::strategy::target_parity_sector);
        has_projected = true;
        write_to_file(StorageReason::PROJ_STATE);
    }
}

void class_algorithm_finite::try_bond_dimension_quench() {
    if(tools::finite::measure::energy_variance(*state) < 10 * settings::precision::variance_convergence_threshold) return;
    if(not settings::strategy::chi_quench_when_stuck) return;
    if(chi_quench_steps > 0) clear_saturation_status();
    if(not state->position_is_any_edge()) return;
    if(chi_quench_steps >= state->get_length()) {
        tools::log->info("Chi quench continues -- {} steps left", chi_quench_steps);
        tools::finite::mps::truncate_all_sites(*state, chi_lim_quench_ahead);
        return;
    }

    //    if(not sim_status.simulation_has_got_stuck) {
    //        tools::log->trace("Chi quench skipped: simulation not stuck");
    //        return;
    //    }
    if(sim_status.simulation_has_stuck_for <= 2) {
        tools::log->info("Chi quench skipped: simulation not been stuck for long enough");
        return;
    }
    if(num_chi_quenches >= max_chi_quenches) {
        tools::log->trace("Chi quench skipped: max number of chi quenches ({}) have been made already", num_chi_quenches);
        return;
    }
    double truncation_threshold = 5 * settings::precision::svd_threshold;
    size_t trunc_bond_count     = state->num_sites_truncated(truncation_threshold);
    size_t bond_at_lim_count    = state->num_bonds_at_limit();
    tools::log->trace("Truncation threshold : {:.4e}", std::pow(truncation_threshold, 2));
    tools::log->trace("Truncation errors    : {}", state->get_truncation_errors());
    tools::log->trace("Bond dimensions      : {}", tools::finite::measure::bond_dimensions(*state));
    tools::log->trace("Entanglement entr    : {}", tools::finite::measure::entanglement_entropies(*state));
    tools::log->trace("Truncated bond count : {} ", trunc_bond_count);
    tools::log->trace("Bonds at limit  count: {} ", bond_at_lim_count);
    if(state->is_bond_limited(truncation_threshold)) {
        tools::log->info("Chi quench skipped: state is bond limited - prefer updating bond dimension");
        return;
    }
    auto bond_dimensions    = tools::finite::measure::bond_dimensions(*state);
    auto max_bond_dimension = *max_element(std::begin(bond_dimensions), std::end(bond_dimensions));
    tools::log->info("Chi quench started");
    tools::finite::mps::truncate_all_sites(*state, max_bond_dimension / 2);
    tools::log->debug("Bond dimensions: {}", tools::finite::measure::bond_dimensions(*state));
    clear_saturation_status();
    chi_quench_steps = 1 * state->get_length();
    num_chi_quenches++;
}

void class_algorithm_finite::try_hamiltonian_perturbation() {
    if(not settings::strategy::perturb_when_stuck) return;
    if(not state->position_is_any_edge()) return;
    if(state->is_perturbed()) return;
    if(perturbation_steps > 0) return;
    if(not sim_status.simulation_has_got_stuck) {
        tools::log->info("Perturbation skipped: simulation not stuck");
        return;
    }
    if(num_perturbations >= max_perturbations) {
        tools::log->info("Perturbation skipped: max number of perturbation trials ({}) have been made already", num_perturbations);
        return;
    }
}

void class_algorithm_finite::try_disorder_damping() {
    if(not state->position_is_left_edge()) return;
    // If there are damping exponents to process, do so
    if(not damping_exponents.empty()) {
        tools::log->info("Setting damping exponent = {}", damping_exponents.back());
        state->damp_hamiltonian(damping_exponents.back(), 0);
        damping_exponents.pop_back();
        if(damping_exponents.empty() and state->is_damped()) throw std::logic_error("Damping trial ended but the state is still damped");
    }
    if(not settings::strategy::damping_when_stuck) return;
    if(state->is_damped()) return;
    if(not sim_status.simulation_has_got_stuck) {
        tools::log->info("Damping skipped: simulation not stuck");
        return;
    }
    if(num_dampings >= max_dampings) {
        tools::log->info("Damping skipped: max number of damping trials ({}) have been made already", num_dampings);
        return;
    }
    // damping_exponents = math::LogSpaced(6,0.0,0.25,0.001);
    // damping_exponents = {0.0,0.1};
    damping_exponents = {0.0, 0.01, 0.02, 0.04, 0.08, 0.16, 0.256};
    tools::log->info("Generating damping exponents = {}", damping_exponents);
    has_damped = true;
    clear_saturation_status();
    tools::log->info("Setting damping exponent = {}", damping_exponents.back());
    state->damp_hamiltonian(damping_exponents.back(), 0);
    damping_exponents.pop_back();
    if(damping_exponents.empty() and state->is_damped()) throw std::logic_error("Damping trial ended but the state is still damped");
}

void class_algorithm_finite::check_convergence_variance(double threshold, double slope_threshold) {
    // Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    if(not state->position_is_any_edge()) {
        return;
    }
    tools::log->debug("Checking convergence of variance mpo");
    threshold       = std::isnan(threshold) ? settings::precision::variance_convergence_threshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::variance_slope_threshold : slope_threshold;
    auto report     = check_saturation_using_slope(V_mpo_vec, X_mpo_vec, tools::finite::measure::energy_variance(*state), sim_status.iter, 1, slope_threshold);
    sim_status.variance_mpo_has_converged = tools::finite::measure::energy_variance(*state) < threshold;
    if(report.has_computed) {
        V_mpo_slopes.emplace_back(report.slope);
        auto last_nonconverged_ptr = std::find_if(V_mpo_vec.rbegin(), V_mpo_vec.rend(), [threshold](auto const &val) { return val > threshold; });
        auto last_nonsaturated_ptr =
            std::find_if(V_mpo_slopes.rbegin(), V_mpo_slopes.rend(), [slope_threshold](auto const &val) { return val > slope_threshold; });
        size_t converged_count                = (size_t) std::distance(V_mpo_vec.rbegin(), last_nonconverged_ptr);
        size_t saturated_count                = (size_t) std::distance(V_mpo_slopes.rbegin(), last_nonsaturated_ptr);
        sim_status.variance_mpo_has_saturated = report.slope < slope_threshold; // or saturated_count >= min_saturation_iters;
        sim_status.variance_mpo_saturated_for = std::max(converged_count, saturated_count);
        tools::log->debug("Variance slope details:");
        tools::log->debug(" -- relative slope    = {} %", report.slope);
        tools::log->debug(" -- tolerance         = {} %", slope_threshold);
        tools::log->debug(" -- last var average  = {} ", report.avgY);
        tools::log->debug(" -- check from        = {} ", report.check_from);
        tools::log->debug(" -- var history       = {} ", V_mpo_vec);
        tools::log->debug(" -- slope history     = {} ", V_mpo_slopes);
        tools::log->debug(" -- has saturated     = {} ", sim_status.variance_mpo_has_saturated);
        tools::log->debug(" -- has saturated for = {} ", sim_status.variance_mpo_saturated_for);
        tools::log->debug(" -- has converged     = {} ", sim_status.variance_mpo_has_converged);
        tools::log->debug(" -- has converged for = {} ", converged_count);
        if(V_mpo_vec.back() < threshold and sim_status.variance_mpo_saturated_for == 0) throw std::logic_error("Variance should have saturated");
        if(V_mpo_vec.back() < threshold and not sim_status.variance_mpo_has_converged) throw std::logic_error("Variance should have converged");
    }
}

void class_algorithm_finite::check_convergence_entg_entropy(double slope_threshold) {
    // Based on the the slope of entanglement entanglement_entropy_midchain
    // This one is cheap to compute.
    if(not state->position_is_any_edge()) {
        return;
    }
    tools::log->debug("Checking convergence of entanglement");

    slope_threshold                         = std::isnan(slope_threshold) ? settings::precision::entropy_slope_threshold : slope_threshold;
    auto                          entropies = tools::finite::measure::entanglement_entropies(*state);
    std::vector<SaturationReport> reports(entropies.size());

    for(size_t site = 0; site < entropies.size(); site++) {
        reports[site] = check_saturation_using_slope(S_mat[site], X_mat[site], entropies[site], sim_status.iter, 1, slope_threshold);
    }
    bool all_computed                     = std::all_of(reports.begin(), reports.end(), [](const SaturationReport r) { return r.has_computed; });
    sim_status.entanglement_has_saturated = false;
    if(all_computed) {
        // idx_max_slope is the index to the site with maximum slope
        auto idx_max_slope = (size_t) std::distance(
            reports.begin(),
            std::max_element(reports.begin(), reports.end(), [](const SaturationReport &r1, const SaturationReport &r2) { return r1.slope < r2.slope; }));
        // idx_max_slope is the index to the site with maximum slope
        //        size_t idx_min_satur = std::distance(reports.begin(),
        //                                             std::min_element(reports.begin(),reports.end(),
        //                                   [](const SaturationReport &r1, const SaturationReport &r2)
        //                                   {return r1.saturated_for < r2.saturated_for;}));

        S_slopes.push_back(reports[idx_max_slope].slope);
        auto last_nonsaturated_ptr = std::find_if(S_slopes.rbegin(), S_slopes.rend(), [slope_threshold](auto const &val) { return val > slope_threshold; });
        auto saturated_count       = (size_t) std::distance(S_slopes.rbegin(), last_nonsaturated_ptr);

        sim_status.entanglement_has_saturated = S_slopes.back() < slope_threshold;
        sim_status.entanglement_saturated_for = saturated_count;
        std::vector<double> all_avergs;
        std::vector<double> all_slopes;
        for(auto &r : reports) all_avergs.push_back(r.avgY);
        for(auto &r : reports) all_slopes.push_back(r.slope);
        tools::log->debug("Max slope of entanglement entropy at site {}: {:.8f} %", idx_max_slope, S_slopes.back());
        tools::log->debug("Entanglement slope details of worst slope:");
        tools::log->debug(" -- site              = {}", idx_max_slope);
        tools::log->debug(" -- relative slope    = {} %", reports[idx_max_slope].slope);
        tools::log->debug(" -- tolerance         = {} %", slope_threshold);
        tools::log->debug(" -- check from        = {} ", reports[idx_max_slope].check_from);
        tools::log->debug(" -- ent history       = {} ", S_mat[idx_max_slope]);
        tools::log->debug(" -- slope history     = {} ", S_slopes);
        tools::log->debug(" -- has saturated     = {} ", sim_status.entanglement_has_saturated);
        tools::log->debug(" -- has saturated for = {} (site {} )", sim_status.entanglement_saturated_for, saturated_count);
        tools::log->debug(" -- all averages      = {} ", all_avergs);
        tools::log->debug(" -- all slopes        = {} ", all_slopes);
        if(reports[idx_max_slope].slope > slope_threshold and sim_status.entanglement_has_saturated) throw std::logic_error("Not supposed to be saturated!!");
    }
    sim_status.entanglement_has_converged = sim_status.entanglement_has_saturated;
}

void class_algorithm_finite::clear_saturation_status() {
    tools::log->trace("Clearing saturation status");
    for(auto &mat : S_mat) {
        mat.clear();
    }
    for(auto &mat : X_mat) {
        mat.clear();
    }
    S_slopes.clear();

    V_mpo_vec.clear();
    X_mpo_vec.clear();
    V_mpo_slopes.clear();

    sim_status.entanglement_has_converged  = false;
    sim_status.entanglement_has_saturated  = false;
    sim_status.entanglement_saturated_for  = 0;
    sim_status.variance_mpo_has_converged  = false;
    sim_status.variance_mpo_has_saturated  = false;
    sim_status.variance_mpo_saturated_for  = 0;
    sim_status.chi_lim_has_reached_chi_max = false;
    sim_status.simulation_has_to_stop      = false;
    sim_status.simulation_has_got_stuck    = false;
    sim_status.simulation_has_converged    = false;
    sim_status.simulation_has_saturated    = false;
    sim_status.simulation_has_succeeded    = false;
    sim_status.simulation_has_stuck_for    = 0;
    has_projected                          = false;
    has_damped                             = false;
}

void class_algorithm_finite::write_to_file(StorageReason storage_reason) {
    if(storage_reason == StorageReason::PROJ_STATE and not has_projected) {
        auto state_projected = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::strategy::target_parity_sector);
        write_to_file(storage_reason, state_projected);
    } else {
        write_to_file(storage_reason, *state);
    }
}

void class_algorithm_finite::write_to_file(StorageReason storage_reason, const class_state_finite &state) {
    StorageLevel storage_level;
    std::string  state_prefix = sim_name + "/" + state_name;
    std::string  model_prefix = sim_name + "/model";
    switch(storage_reason) {
        case StorageReason::JOURNAL: {
            if(not state.position_is_any_edge()) return;
            if(math::mod(sim_status.iter, write_freq()) != 0) return;
            if(settings::output::storage_level_journal == StorageLevel::NONE) return;
            state_prefix  = state_prefix + "/journal";
            storage_level = settings::output::storage_level_journal;
            if(settings::output::journal_keep_only_last_iter)
                state_prefix.append("/iter_last");
            else
                state_prefix.append("/iter_" + std::to_string(sim_status.iter));
            break;
        }
        case StorageReason::RESULTS: {
            if(settings::output::storage_level_results == StorageLevel::NONE) return;
            storage_level = settings::output::storage_level_results;
            state_prefix  = state_prefix + "/results";
            break;
        }
        case StorageReason::CHI_UPDATE: {
            if(settings::output::storage_level_chi_update == StorageLevel::NONE) return;
            if(not chi_grow()) return;
            storage_level = settings::output::storage_level_chi_update;
            state_prefix  = state_prefix + "/results_chi_" + std::to_string(state.get_chi_lim());
            break;
        }
        case StorageReason::PROJ_STATE: {
            if(settings::output::storage_level_proj_state == StorageLevel::NONE) return;
            storage_level = settings::output::storage_level_proj_state;
            state_prefix  = state_prefix + "/projection";
            break;
        }
        case StorageReason::INIT_STATE: {
            if(settings::output::storage_level_init_state == StorageLevel::NONE) return;
            storage_level = settings::output::storage_level_init_state;
            state_prefix  = state_prefix + "/state_init";
            break;
        }
        case StorageReason::EMIN_STATE: {
            if(settings::output::storage_level_emin_state == StorageLevel::NONE) return;
            storage_level = settings::output::storage_level_emin_state;
            state_prefix  = sim_name + "/state_emin";
            break;
        }
        case StorageReason::EMAX_STATE: {
            if(settings::output::storage_level_emax_state == StorageLevel::NONE) return;
            storage_level = settings::output::storage_level_emax_state;
            state_prefix  = sim_name + "/state_emax";
            break;
        }
    }

    if(state_prefix.empty()) throw std::runtime_error("Prefix is empty");
    tools::finite::io::h5dset::write_mps(*h5pp_file, state_prefix, storage_level, state);
    tools::finite::io::h5dset::write_mpo(*h5pp_file, model_prefix, storage_level, state);
    tools::finite::io::h5dset::write_measurements(*h5pp_file, state_prefix, storage_level, state);
    tools::finite::io::h5table::write_measurements(*h5pp_file, state_prefix, storage_level, state,sim_status);
    tools::finite::io::h5table::write_sim_status(*h5pp_file, state_prefix, storage_level, sim_status);
    tools::finite::io::h5table::write_profiling(*h5pp_file, state_prefix, storage_level, sim_status);
    tools::finite::io::h5table::write_mem_usage(*h5pp_file, state_prefix, storage_level, sim_status);
    tools::common::io::h5attr::write_meta(*h5pp_file,sim_name, state_prefix, model_prefix, settings::model::model_type, storage_level, sim_status);
}

void class_algorithm_finite::copy_from_tmp(StorageReason storage_reason) {
    if(not h5pp_file) return;
    if(not settings::output::use_temp_dir) return;
    if(not state->position_is_any_edge()) return;
    if(settings::output::storage_level_results == StorageLevel::NONE and settings::output::storage_level_chi_update == StorageLevel::NONE and
       settings::output::storage_level_journal == StorageLevel::NONE and settings::output::storage_level_proj_state == StorageLevel::NONE)
        return;
    switch(storage_reason) {
        case StorageReason::JOURNAL:
            if(math::mod(sim_status.iter, settings::output::copy_from_temp_freq) != 0) return; // Check that we write according to the frequency given
        case StorageReason::RESULTS:
        case StorageReason::CHI_UPDATE:
        case StorageReason::PROJ_STATE:
        case StorageReason::INIT_STATE:
        case StorageReason::EMIN_STATE:
        case StorageReason::EMAX_STATE:break;
    }
    tools::common::io::h5tmp::copy_from_tmp(h5pp_file->getFilePath());
}

void class_algorithm_finite::print_status_update() {
    if(math::mod(sim_status.step, print_freq()) != 0) return;
    if(print_freq() == 0) return;

    using namespace tools::finite::measure;
    std::string report;
    report += fmt::format("{:<} ", sim_name);
    report += fmt::format("iter: {:<4} ", sim_status.iter);
    report += fmt::format("step: {:<5} ", sim_status.step);
    report += fmt::format("L: {} ", state->get_length());
    if(state->active_sites.empty())
        report += fmt::format("l: {:<2} ", state->get_position());
    else if(state->get_direction() > 0)
        report += fmt::format("l: [{:>2}-{:<2}] ", state->active_sites.front(), state->active_sites.back());
    else if(state->get_direction() < 0)
        report += fmt::format("l: [{:>2}-{:<2}] ", state->active_sites.back(), state->active_sites.front());
    report += fmt::format("E/L: {:<20.16f} ", energy_per_site(*state));
    if(sim_type == SimulationType::xDMRG) {
        report += fmt::format("ε: {:<6.4f} ", sim_status.energy_dens);
    }
    report += fmt::format("Sₑ(l): {:<10.8f} ", entanglement_entropy_current(*state));
    report += fmt::format("log₁₀ σ²(E)/L: {:<10.6f} [{:<10.6f}] ", std::log10(energy_variance_per_site(*state)),
                          std::log10(state->lowest_recorded_variance / static_cast<double>(state->get_length())));
    report += fmt::format("χmax: {:<3} χlim: {:<3} χ: {:<3} ", chi_max(), state->get_chi_lim(), bond_dimension_current(*state));
    report += fmt::format("log₁₀ trunc: {:<10.4f} ", std::log10(state->get_truncation_error(state->get_position())));
    report += fmt::format("stk: {:<1} ", sim_status.simulation_has_stuck_for);
    report += fmt::format("sat: [σ² {:<1} Sₑ {:<1}] ", sim_status.variance_mpo_saturated_for, sim_status.entanglement_saturated_for);
    report += fmt::format("con: {:<5} ", sim_status.simulation_has_converged);
    report += fmt::format("time: {:<8.2f}s ", tools::common::profile::t_tot->get_age());
    report += fmt::format("mem MB: [Rss {:<.1f} Peak {:<.1f} Vm {:<.1f}] ", tools::common::profile::mem_rss_in_mb(), tools::common::profile::mem_hwm_in_mb(),
                          tools::common::profile::mem_vm_in_mb());
    tools::log->info(report);
}

void class_algorithm_finite::print_status_full() {
    using namespace tools::finite::measure;
    tools::log->info("{:=^60}", "");
    tools::log->info("= {: ^56} =", "Final results [" + sim_name + "]");
    tools::log->info("{:=^60}", "");
    tools::log->info("--- Final results  --- {} ---", sim_name);
    tools::log->info("Sites                              = {}", state->get_length());
    tools::log->info("Iterations (sweeps)                = {}", sim_status.iter);
    tools::log->info("Steps                              = {}", sim_status.step);
    tools::log->info("Total time                         = {:<.1f} s = {:<.2f} min", tools::common::profile::t_tot->get_age(),
                     tools::common::profile::t_tot->get_age() / 60);
    tools::log->info("Energy per site E/L                = {:<.16f}", energy_per_site(*state));
    if(sim_type == SimulationType::xDMRG) {
        tools::log->info("Energy density (rescaled 0 to 1) ε = {:<6.4f}", energy_normalized(*state, sim_status));
    }
    tools::log->info("Variance per site log₁₀ σ²(E)/L    = {:<.16f}", std::log10(energy_variance_per_site(*state)));
    tools::log->info("Bond dimension maximum χmax        = {}", chi_max());
    tools::log->info("Bond dimensions χ                  = {}", bond_dimensions(*state));
    tools::log->info("Bond dimension  χ (mid)            = {}", bond_dimension_midchain(*state));
    tools::log->info("Entanglement entropies Sₑ          = {}", entanglement_entropies(*state));
    tools::log->info("Entanglement entropiy Sₑ (mid)     = {}", entanglement_entropy_midchain(*state));
    tools::log->info("Truncation Errors                  = {}", state->get_truncation_errors());
    tools::log->info("Simulation converged               = {:<}", sim_status.simulation_has_converged);
    tools::log->info("Simulation saturated               = {:<}", sim_status.simulation_has_saturated);
    tools::log->info("Simulation succeeded               = {:<}", sim_status.simulation_has_succeeded);
    tools::log->info("Simulation got stuck               = {:<}", sim_status.simulation_has_got_stuck);
    tools::log->info("σ² slope                           = {:<8.4f} %   Converged : {:<8}  Saturated: {:<8}", V_mpo_slopes.back(),
                     sim_status.variance_mpo_has_converged, sim_status.variance_mpo_has_saturated);
    tools::log->info("Sₑ slope                           = {:<8.4f} %   Converged : {:<8}  Saturated: {:<8}", S_slopes.back(),
                     sim_status.entanglement_has_converged, sim_status.entanglement_has_saturated);
}
