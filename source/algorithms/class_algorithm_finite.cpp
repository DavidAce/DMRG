//
// Created by david on 2019-06-24.
//

#include "class_algorithm_finite.h"
#include <h5pp/h5pp.h>
#include <math/nmspc_math.h>
#include <config/nmspc_settings.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/edges/class_edges_finite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite.h>

class_algorithm_finite::class_algorithm_finite(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type)
    : class_algorithm_base(std::move(h5ppFile_), algo_type), tensors(settings::model::model_type,settings::model::model_size,0) {
    tools::log->trace("Constructing class_algorithm_finite");
    tools::finite::mpo::randomize(*tensors.model);
    tools::finite::mps::random_product_state(*tensors.state, settings::strategy::initial_parity_sector, settings::input::bitfield,
                                             settings::strategy::use_pauli_eigvecs);
    tools::finite::debug::check_integrity(*tensors.state,*tensors.model,*tensors.edges);
    S_mat.resize(tensors.get_length() + 1);
    X_mat.resize(tensors.get_length() + 1);
    tools::finite::print::dimensions(*tensors.model);
    tools::finite::print::dimensions(tensors);
}


// We need to make a destructor manually for the enclosing class "class_model_finite"
// that encloses "class_model_base". Otherwise unique_ptr will forcibly inline its
// own default deleter.
// This allows us to forward declare the abstract base class "class_model_base"
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_algorithm_finite::~class_algorithm_finite() = default;




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
    tools::log->info("Starting {}", algo_name);
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
            auto state_prefix = tools::common::io::h5resume::find_resumable_state(*h5pp_file, algo_name);
            if(state_prefix.empty()) throw std::runtime_error("Could not resume: no valid resume candidates found");
            tools::log->info("Resuming state [{}]", state_prefix);
            tools::finite::io::h5resume::load_tensors(*h5pp_file, state_prefix, tensors, status);
            // Now we decide what to do
            // Let's consider case 1

            //            if(status.algorithm_has_succeeded){
            //
            //            }
            exit(0);
            while(h5pp_file->linkExists(algo_name + state_name)) state_name = "state_" + std::to_string(state_number++);
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
    //    if(not algo_on()) return;
    //    if(not settings::input::h5_load_filename.empty()){
    //        h5pp::File h5pp_load(settings::input::h5_load_filename,h5pp::AccessMode::READONLY,h5pp::CreateMode::OPEN);
    //        tools::finite::io::h5restore::load_tensors(h5pp_load, sim_name, status, *tensors.state);
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
    //                tools::finite::io::h5resume::load_tensors(*h5pp_file, sim_name, status, *tensors.state);
    //            } catch(std::exception &ex) {
    //                tools::log->error("Failed to load from output: {}", ex.what());
    //                throw std::runtime_error("Failed to resume from file: " + std::string(ex.what()));
    //            } catch(...) {
    //                tools::log->error("Unknown error when trying to resume from file.");
    //            }
    //
    //            bool convergence_was_reached = h5pp_file->readDataset<bool>(sim_name + "/status/simulation_has_converged");
    //            if(not convergence_was_reached) {
    //                // Case 1 c -- resume simulation, reset the number of sweeps first.
    //                tools::log->trace("Case 1c");
    //                settings::xdmrg::max_sweeps += tensors.state->get_iteration();
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
    if(not algo_on()) {
        return;
    }
    tools::log->info("Starting {}", algo_name);
    tools::common::profile::t_tot->tic();
    if(h5pp_file) {
        // This is case 1
        bool finOK_exists = h5pp_file->linkExists("common/finished_all");
        bool mps_exists   = h5pp_file->linkExists(algo_name + "/state/mps");
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
                tools::finite::io::h5resume::load_tensors(*h5pp_file, algo_name, tensors, status);
            } catch(std::exception &ex) {
                tools::log->error("Failed to load from output: {}", ex.what());
                throw std::runtime_error("Failed to resume from file: " + std::string(ex.what()));
            } catch(...) {
                tools::log->error("Unknown error when trying to resume from file.");
            }

            bool convergence_was_reached = h5pp_file->readDataset<bool>(algo_name + "/status/simulation_has_converged");
            if(not convergence_was_reached) {
                // Case 1 c -- resume simulation, reset the number of sweeps first.
                tools::log->trace("Case 1c");
                settings::xdmrg::max_iters += status.iter;
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
    tools::log->info("Running {} preprocessing (base)", algo_name);
    tools::common::profile::t_pre->tic();
    tools::finite::io::h5table::write_model(*h5pp_file, algo_name + "/model", settings::output::storage_level_model, *tensors.model);
    tools::finite::io::h5dset::write_model(*h5pp_file, algo_name + "/model", settings::output::storage_level_model, *tensors.model);
    tensors.state->set_chi_max(chi_max());
    status.chi_max = chi_max();
    update_bond_dimension_limit(chi_init());
    tools::common::profile::t_pre->toc();
    tools::log->info("Finished {} preprocessing (base)", algo_name);
}

void class_algorithm_finite::run_postprocessing() {
    tools::log->info("Running {} postprocessing", algo_name);
    tools::common::profile::t_pos->tic();
    status.algorithm_has_finished = true;
    tools::finite::debug::check_integrity(*tensors.state);
    tensors.clear_measurements();
    tensors.state->clear_measurements();
    write_to_file(StorageReason::FINISHED);
    if(not has_projected) write_to_file(StorageReason::PROJ_STATE);
    print_status_full();
    tools::common::profile::t_pos->toc();
    tools::log->info("Finished {} postprocessing", algo_name);
}

void class_algorithm_finite::move_center_point(std::optional<size_t> num_moves) {
    if(not num_moves.has_value()) {
        if(tensors.state->active_sites.empty())
            num_moves = 1ul;
        else if(settings::precision::move_multisite == "one")
            num_moves = 1ul;
        else if(settings::precision::move_multisite == "mid")
            num_moves = std::max(1ul, (tensors.state->active_sites.size()) / 2);
        else if(settings::precision::move_multisite == "max")
            num_moves = std::max(1ul, tensors.state->active_sites.size() - 2ul);
        else {
            throw std::logic_error("Specify how many sites multisite should move! Expected one of {one,mid,max}, got [" +
                                   settings::precision::move_multisite + "]");
        }
    }

    tools::log->trace("Moving center point {} steps in direction {}", num_moves.value(), tensors.state->get_direction());
    tensors.clear_cache();
    try {
        for(size_t i = 0; i < num_moves.value(); i++) {
            if(tensors.state->position_is_any_edge()) {
                tensors.state->increment_iter();
                has_projected = false;
            }
            tools::finite::mps::move_center_point(*tensors.state);
            tensors.state->increment_step();
            if(chi_quench_steps > 0) chi_quench_steps--;
        }
        tools::finite::debug::check_integrity(*tensors.state, *tensors.model, *tensors.edges);
        tensors.state->active_sites.clear();
    } catch(std::exception &e) {
        tools::finite::print::dimensions(tensors);
        throw std::runtime_error("Failed to move center point: " + std::string(e.what()));
    }
    status.iter      = tensors.state->get_iteration();
    status.step      = tensors.state->get_step();
    status.position  = tensors.state->get_position();
    status.direction = tensors.state->get_direction();
}

void class_algorithm_finite::update_truncation_limit() {
    if(not tensors.state->position_is_any_edge()) return;
    // Will update SVD threshold iff the energy variance is being limited by truncation error
    size_t bond_at_lim_count = tensors.state->num_bonds_at_limit();
    if(bond_at_lim_count > 0) return; // Return because we should rather increase bond dimension than lower the svd threshold
    tools::log->info("Truncated variances: {}", tensors.state->get_truncated_variances());
    //
    //    double truncation_threshold = settings::precision::svd_threshold;
    //    size_t trunc_bond_count  = tensors.state->num_sites_truncated(truncation_threshold);
    //    if(trunc_bond_count > 0) settings::precision::svd_threshold *= 0.5;
    //    tools::log->info("Lowered SVD threshold to {}",settings::precision::svd_threshold);
}

void class_algorithm_finite::update_bond_dimension_limit(std::optional<long> tmp_bond_limit) {
    if(tmp_bond_limit.has_value()) {
        tensors.state->set_chi_lim(tmp_bond_limit.value());
        status.chi_lim = tmp_bond_limit.value();
        return;
    }
    try {
        long chi_lim_now = tensors.state->get_chi_lim();
        if(chi_lim_now < chi_init()) throw std::logic_error("Chi limit should be larger than chi init");
    } catch(std::exception &error) {
        // If we reached this stage, either
        // 1) chi_lim is not initialized yet
        // 2) chi_lim is initialized, but it is smaller than the init value found in settings
        // Either way, we should set chi_lim to be chi_init, unless chi_init is larger than tmp_bond_limit
        tools::log->info("Setting initial bond dimension limit: {}", chi_init());
        tensors.state->set_chi_lim(chi_init());
        status.chi_lim = chi_init();
        return;
    }
    if(not tensors.state->position_is_any_edge()) return;
    status.chi_lim_has_reached_chi_max = tensors.state->get_chi_lim() >= chi_max();
    if(not status.chi_lim_has_reached_chi_max) {
        if(chi_grow()) {
            status.chi_lim = tensors.state->get_chi_lim();

            // If we got here we want grow the bond dimension limit progressively during the simulation
            // Only increment the bond dimension if
            // * No experiments are on-going like perturbation or damping
            // * the simulation is stuck
            // * the state is limited by bond dimension
            if(tensors.model->is_damped()) {
                tools::log->info("State is undergoing disorder damping -- cannot increase bond dimension yet");
                return;
            }
            if(tensors.model->is_perturbed()) {
                tools::log->info("State is undergoing perturbation -- cannot increase bond dimension yet");
                return;
            }
            if(not status.algorithm_has_got_stuck and tensors.state->get_chi_lim() >= 16) {
                tools::log->info("State is not stuck yet. Kept current limit {}", tensors.state->get_chi_lim());
                return;
            }
            //            if(status.algorithm_has_stuck_for <= 1 and tensors.state->get_chi_lim() >= 16){
            //                tools::log->info("State has not been stuck for long enough. Kept current limit {}", tensors.state->get_chi_lim());
            //                return;
            //            }
            if(tools::log->level() <= spdlog::level::info) {
                double truncation_threshold = 2 * settings::precision::svd_threshold;
                size_t trunc_bond_count     = tensors.state->num_sites_truncated(truncation_threshold);
                size_t bond_at_lim_count    = tensors.state->num_bonds_at_limit();
                tools::log->info("Truncation threshold  : {:<.8e}", truncation_threshold);
                tools::log->info("Truncation errors     : {}", tensors.state->get_truncation_errors());
                tools::log->info("Bond dimensions       : {}", tools::finite::measure::bond_dimensions(*tensors.state));
                tools::log->info("Truncated bond count  : {} ", trunc_bond_count);
                tools::log->info("Bonds at limit  count : {} ", bond_at_lim_count);
                tools::log->info("Entanglement entropies: {} ", tools::finite::measure::entanglement_entropies(*tensors.state));
            }
            bool state_is_bond_limited = tensors.state->is_bond_limited(5 * settings::precision::svd_threshold);
            if(not state_is_bond_limited) {
                tools::log->info("State is not limited by its bond dimension. Kept current limit {}", tensors.state->get_chi_lim());
                return;
            }

            // Write current results before updating bond dimension
            write_to_file(StorageReason::CHI_UPDATE);
            long chi_lim_new = std::min(tensors.state->get_chi_max(), tensors.state->get_chi_lim() * 2);
            if(tensors.state->get_chi_lim() < 16) chi_lim_new = std::min(tensors.state->get_chi_max(), tensors.state->get_chi_lim() + 1);
            tools::log->info("Updating bond dimension limit {} -> {}", tensors.state->get_chi_lim(), chi_lim_new);
            tensors.state->set_chi_lim(chi_lim_new);
            clear_saturation_status();
            status.chi_lim_has_reached_chi_max = tensors.state->get_chi_lim() == chi_max();
            if(status.chi_lim_has_reached_chi_max and has_projected) has_projected = false;
            tools::log->info("Projecting at site {} to direction {} after updating bond dimension to χ = {} ", tensors.state->get_position(),
                             settings::strategy::target_parity_sector, chi_lim_new);
            copy_from_tmp(StorageReason::CHI_UPDATE);
            if(not tensors.state->position_is_any_edge()) throw std::runtime_error("Update bond dimension: no longer at edge!");
            if(settings::strategy::project_on_chi_update) {
                *tensors.state        = tools::finite::ops::get_projection_to_closest_parity_sector(*tensors.state, settings::strategy::target_parity_sector);
                has_projected = true;
            }
            write_to_file(StorageReason::PROJ_STATE);
            if(settings::strategy::randomize_on_chi_update and tensors.state->get_chi_lim() >= 32) reset_to_random_current_state();
        } else {
            // Here the settings specify to just set the limit to maximum chi directly
            tools::log->info("Setting bond dimension limit to maximum = {}", chi_max());
            tensors.state->set_chi_lim(chi_max());
        }
    } else {
        tools::log->debug("Chi limit has reached max: {} -> Kept current bond dimension limit {}", chi_max(), tensors.state->get_chi_lim());
    }
    status.chi_lim = tensors.state->get_chi_lim();
    if(tensors.state->get_chi_lim() > tensors.state->get_chi_max())
        throw std::runtime_error(fmt::format("chi_lim is larger than chi_max! {} > {}", tensors.state->get_chi_lim(), tensors.state->get_chi_max()));
}

void class_algorithm_finite::reset_to_initial_state() {
    tools::log->trace("Resetting MPS to initial product state in parity sector: {}, state number {}", settings::strategy::initial_parity_sector,
                      settings::input::bitfield, settings::strategy::use_pauli_eigvecs);
    if(tensors.state->get_length() != settings::model::model_size) throw std::range_error("System size mismatch");
    // Initialize state
    tools::finite::mps::random_product_state(*tensors.state, settings::strategy::initial_parity_sector, settings::input::bitfield,
                                             settings::strategy::use_pauli_eigvecs);
    clear_saturation_status();
    tensors.state->lowest_recorded_variance = 1;
    status.iter                 = tensors.state->reset_iter();
}

void class_algorithm_finite::reset_to_random_product_state(const std::string &parity_sector) {
    tools::log->trace("Resetting MPS to random product state in parity sector: {}", parity_sector);
    // Randomize state
    tools::finite::mps::random_product_state(*tensors.state, parity_sector, -1, settings::strategy::use_pauli_eigvecs);
    clear_saturation_status();
    tensors.state->lowest_recorded_variance = 1;
    status.iter                 = tensors.state->reset_iter();
    status.step                 = tensors.state->reset_step();
    auto spin_components        = tools::finite::measure::spin_components(*tensors.state);
    tools::log->info("Successfully reset to product state with global spin components: {}", spin_components);
}

void class_algorithm_finite::reset_to_random_current_state(std::optional<double> chi_lim) {
    if(not tensors.position_is_any_edge()) return;
    tools::log->info("Resetting MPS by flipping random spins on current state");
    if(tensors.get_length() != settings::model::model_size) throw std::range_error("System size mismatch");
    tools::log->debug("Bond dimensions: {}", tools::finite::measure::bond_dimensions(*tensors.state));
    // Randomize state
    tools::log->info("Flipping random spins");
    //    tools::finite::mps::random_current_state(*tensors.state,"x");
    tools::finite::mps::random_current_state(*tensors.state, "x", "z");

    // Truncate even more on explicit request
    if(chi_lim) {
        size_t chi_lim_parsed = chi_lim.value() < 1 ? (size_t)(chi_lim.value() * (double) tensors.state->find_largest_chi()) : (size_t) chi_lim.value();
        tools::finite::mps::truncate_all_sites(*tensors.state, chi_lim_parsed);
    }
    tools::log->debug("Bond dimensions: {}", tools::finite::measure::bond_dimensions(*tensors.state));

    clear_saturation_status();
    tensors.state->lowest_recorded_variance = 1;
    status.iter                 = tensors.state->reset_iter();
    auto spin_components            = tools::finite::measure::spin_components(*tensors.state);
    tools::log->info("Successfully reset to random state based on current state. New components: {}", spin_components);
    if(not tensors.state->position_is_any_edge()) throw std::runtime_error("Update bond dimension: no longer at edge!");
}

void class_algorithm_finite::try_projection(class_state_finite & state) {
    if(not state.position_is_any_edge()) return;
    if(has_projected) return;
    if(settings::strategy::project_on_every_sweep or (settings::strategy::project_when_stuck and status.algorithm_has_got_stuck)) {
        tools::log->info("Trying projection to {}", settings::strategy::target_parity_sector);
        *tensors.state        = tools::finite::ops::get_projection_to_closest_parity_sector(*tensors.state, settings::strategy::target_parity_sector);
        has_projected = true;
        write_to_file(StorageReason::PROJ_STATE);
    }
}

void class_algorithm_finite::try_bond_dimension_quench(class_state_finite & state) {
    if(tools::finite::measure::energy_variance(*tensors.state, *tensors.model, *tensors.edges) < 10 * settings::precision::variance_convergence_threshold) return;
    if(not settings::strategy::chi_quench_when_stuck) return;
    if(chi_quench_steps > 0) clear_saturation_status();
    if(not state.position_is_any_edge()) return;
    if(chi_quench_steps >= state.get_length()) {
        tools::log->info("Chi quench continues -- {} steps left", chi_quench_steps);
        tools::finite::mps::truncate_all_sites(*tensors.state, chi_lim_quench_ahead);
        return;
    }

    //    if(not status.algorithm_has_got_stuck) {
    //        tools::log->trace("Chi quench skipped: simulation not stuck");
    //        return;
    //    }
    if(status.algorithm_has_stuck_for <= 2) {
        tools::log->info("Chi quench skipped: simulation not been stuck for long enough");
        return;
    }
    if(num_chi_quenches >= max_chi_quenches) {
        tools::log->trace("Chi quench skipped: max number of chi quenches ({}) have been made already", num_chi_quenches);
        return;
    }
    double truncation_threshold = 5 * settings::precision::svd_threshold;
    size_t trunc_bond_count     = state.num_sites_truncated(truncation_threshold);
    size_t bond_at_lim_count    = state.num_bonds_at_limit();
    tools::log->trace("Truncation threshold : {:.4e}", std::pow(truncation_threshold, 2));
    tools::log->trace("Truncation errors    : {}", state.get_truncation_errors());
    tools::log->trace("Bond dimensions      : {}", tools::finite::measure::bond_dimensions(*tensors.state));
    tools::log->trace("Entanglement entr    : {}", tools::finite::measure::entanglement_entropies(*tensors.state));
    tools::log->trace("Truncated bond count : {} ", trunc_bond_count);
    tools::log->trace("Bonds at limit  count: {} ", bond_at_lim_count);
    if(state.is_bond_limited(truncation_threshold)) {
        tools::log->info("Chi quench skipped: state is bond limited - prefer updating bond dimension");
        return;
    }
    auto bond_dimensions    = tools::finite::measure::bond_dimensions(*tensors.state);
    auto max_bond_dimension = *max_element(std::begin(bond_dimensions), std::end(bond_dimensions));
    tools::log->info("Chi quench started");
    tools::finite::mps::truncate_all_sites(*tensors.state, max_bond_dimension / 2);
    tools::log->debug("Bond dimensions: {}", tools::finite::measure::bond_dimensions(*tensors.state));
    clear_saturation_status();
    chi_quench_steps = 1 * state.get_length();
    num_chi_quenches++;
}

void class_algorithm_finite::try_hamiltonian_perturbation(class_state_finite & state) {
    if(not settings::strategy::perturb_when_stuck) return;
    if(not state.position_is_any_edge()) return;
    if(tensors.model->is_perturbed()) return;
    if(perturbation_steps > 0) return;
    if(not status.algorithm_has_got_stuck) {
        tools::log->info("Perturbation skipped: simulation not stuck");
        return;
    }
    if(num_perturbations >= max_perturbations) {
        tools::log->info("Perturbation skipped: max number of perturbation trials ({}) have been made already", num_perturbations);
        return;
    }
}

void class_algorithm_finite::try_disorder_damping(class_model_finite & model) {
//    if(not state.position_is_left_edge()) return;
    // If there are damping exponents to process, do so
    if(not damping_exponents.empty()) {
        tools::log->info("Setting damping exponent = {}", damping_exponents.back());
        model.damp_hamiltonian(damping_exponents.back(), 0);
        damping_exponents.pop_back();
        if(damping_exponents.empty() and model.is_damped()) throw std::logic_error("Damping trial ended but the state is still damped");
    }
    if(not settings::strategy::damping_when_stuck) return;
    if(model.is_damped()) return;
    if(not status.algorithm_has_got_stuck) {
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
    model.damp_hamiltonian(damping_exponents.back(), 0);
    damping_exponents.pop_back();
    if(damping_exponents.empty() and model.is_damped()) throw std::logic_error("Damping trial ended but the state is still damped");
}

void class_algorithm_finite::check_convergence_variance(double threshold, double slope_threshold) {
    // Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
//    if(not tensors.state->position_is_any_edge()) {
//        return;
//    }
    tools::log->debug("Checking convergence of variance mpo");
    threshold       = std::isnan(threshold) ? settings::precision::variance_convergence_threshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::variance_slope_threshold : slope_threshold;
    auto report     = check_saturation_using_slope(V_mpo_vec, X_mpo_vec, tools::finite::measure::energy_variance(tensors), status.iter, 1, slope_threshold);
    status.variance_mpo_has_converged = tools::finite::measure::energy_variance(*tensors.state, *tensors.model, *tensors.edges) < threshold;
    if(report.has_computed) {
        V_mpo_slopes.emplace_back(report.slope);
        auto last_nonconverged_ptr = std::find_if(V_mpo_vec.rbegin(), V_mpo_vec.rend(), [threshold](auto const &val) { return val > threshold; });
        auto last_nonsaturated_ptr =
            std::find_if(V_mpo_slopes.rbegin(), V_mpo_slopes.rend(), [slope_threshold](auto const &val) { return val > slope_threshold; });
        size_t converged_count                = (size_t) std::distance(V_mpo_vec.rbegin(), last_nonconverged_ptr);
        size_t saturated_count                = (size_t) std::distance(V_mpo_slopes.rbegin(), last_nonsaturated_ptr);
        status.variance_mpo_has_saturated = report.slope < slope_threshold; // or saturated_count >= min_saturation_iters;
        status.variance_mpo_saturated_for = std::max(converged_count, saturated_count);
        tools::log->debug("Variance slope details:");
        tools::log->debug(" -- relative slope    = {} %", report.slope);
        tools::log->debug(" -- tolerance         = {} %", slope_threshold);
        tools::log->debug(" -- last var average  = {} ", report.avgY);
        tools::log->debug(" -- check from        = {} ", report.check_from);
        tools::log->debug(" -- var history       = {} ", V_mpo_vec);
        tools::log->debug(" -- slope history     = {} ", V_mpo_slopes);
        tools::log->debug(" -- has saturated     = {} ", status.variance_mpo_has_saturated);
        tools::log->debug(" -- has saturated for = {} ", status.variance_mpo_saturated_for);
        tools::log->debug(" -- has converged     = {} ", status.variance_mpo_has_converged);
        tools::log->debug(" -- has converged for = {} ", converged_count);
        if(V_mpo_vec.back() < threshold and status.variance_mpo_saturated_for == 0) throw std::logic_error("Variance should have saturated");
        if(V_mpo_vec.back() < threshold and not status.variance_mpo_has_converged) throw std::logic_error("Variance should have converged");
    }
}

void class_algorithm_finite::check_convergence_entg_entropy(double slope_threshold) {
    // Based on the the slope of entanglement entanglement_entropy_midchain
    // This one is cheap to compute.
    if(not tensors.state->position_is_any_edge()) {
        return;
    }
    tools::log->debug("Checking convergence of entanglement");

    slope_threshold                         = std::isnan(slope_threshold) ? settings::precision::entropy_slope_threshold : slope_threshold;
    auto                          entropies = tools::finite::measure::entanglement_entropies(*tensors.state);
    std::vector<SaturationReport> reports(entropies.size());

    for(size_t site = 0; site < entropies.size(); site++) {
        reports[site] = check_saturation_using_slope(S_mat[site], X_mat[site], entropies[site], status.iter, 1, slope_threshold);
    }
    bool all_computed                     = std::all_of(reports.begin(), reports.end(), [](const SaturationReport r) { return r.has_computed; });
    status.entanglement_has_saturated = false;
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

        status.entanglement_has_saturated = S_slopes.back() < slope_threshold;
        status.entanglement_saturated_for = saturated_count;
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
        tools::log->debug(" -- has saturated     = {} ", status.entanglement_has_saturated);
        tools::log->debug(" -- has saturated for = {} (site {} )", status.entanglement_saturated_for, saturated_count);
        tools::log->debug(" -- all averages      = {} ", all_avergs);
        tools::log->debug(" -- all slopes        = {} ", all_slopes);
        if(reports[idx_max_slope].slope > slope_threshold and status.entanglement_has_saturated) throw std::logic_error("Not supposed to be saturated!!");
    }
    status.entanglement_has_converged = status.entanglement_has_saturated;
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

    status.entanglement_has_converged  = false;
    status.entanglement_has_saturated  = false;
    status.entanglement_saturated_for  = 0;
    status.variance_mpo_has_converged  = false;
    status.variance_mpo_has_saturated  = false;
    status.variance_mpo_saturated_for  = 0;
    status.chi_lim_has_reached_chi_max = false;
    status.algorithm_has_to_stop      = false;
    status.algorithm_has_got_stuck    = false;
    status.algorithm_has_converged    = false;
    status.algorithm_has_saturated    = false;
    status.algorithm_has_succeeded    = false;
    status.algorithm_has_stuck_for    = 0;
    has_projected                          = false;
    has_damped                             = false;
}

void class_algorithm_finite::write_to_file(StorageReason storage_reason) {
    if(storage_reason == StorageReason::PROJ_STATE and not has_projected) {
        auto state_projected = tools::finite::ops::get_projection_to_closest_parity_sector(*tensors.state, settings::strategy::target_parity_sector);
        write_to_file(storage_reason, state_projected);
    } else {
        write_to_file(storage_reason, *tensors.state);
    }
}

void class_algorithm_finite::write_to_file(StorageReason storage_reason, const class_state_finite &state) {
    StorageLevel      storage_level;
    const std::string table_prefix = algo_name + "/" + state_name;
    std::string       state_prefix = algo_name + "/" + state_name;  // May get modified
    std::string       model_prefix = algo_name + "/model";
    switch(storage_reason) {
        case StorageReason::FINISHED: {
            if(status.algorithm_has_succeeded)
                storage_level = settings::output::storage_level_good_state;
            else
                storage_level = settings::output::storage_level_fail_state;
            state_prefix.append("/finished");
            break;
        }
        case StorageReason::CHECKPOINT: {
            if(not state.position_is_any_edge()) return;
            if(math::mod(status.iter, settings::output::checkpoint_frequency) != 0) return;
            state_prefix.append("/checkpoint");
            storage_level = settings::output::storage_level_checkpoint;
            if(settings::output::checkpoint_keep_newest_only)
                state_prefix.append("/iter_last");
            else
                state_prefix.append("/iter_" + std::to_string(status.iter));
            break;
        }
        case StorageReason::CHI_UPDATE: {
            if(not chi_grow()) return;
            storage_level = settings::output::storage_level_checkpoint;
            state_prefix.append("/checkpoint");
            state_prefix.append("/chi_" + std::to_string(state.get_chi_lim()));
            break;
        }
        case StorageReason::PROJ_STATE: {
            storage_level = settings::output::storage_level_proj_state;
            state_prefix.append("/projection");
            break;
        }
        case StorageReason::INIT_STATE: {
            storage_level = settings::output::storage_level_init_state;
            state_prefix.append("/state_init");
            break;
        }
        case StorageReason::EMIN_STATE: {
            storage_level = settings::output::storage_level_emin_state;
            state_prefix  = algo_name + "/state_emin";
            break;
        }
        case StorageReason::EMAX_STATE: {
            storage_level = settings::output::storage_level_emax_state;
            state_prefix  = algo_name + "/state_emax";
            break;
        }
    }
    if(storage_level == StorageLevel::NONE) return;
    if(state_prefix.empty()) throw std::runtime_error("State prefix is empty");
    tools::finite::io::h5dset::write_state(*h5pp_file, state_prefix, storage_level, state);
    tools::finite::io::h5dset::write_ent(*h5pp_file, state_prefix, storage_level, state);
    tools::common::io::h5attr::write_meta(*h5pp_file, algo_name, state_prefix, model_prefix, settings::model::model_type, storage_level, status);


    // The main results have now been written. Next we append data to tables
    // Some storage reasons should not do this however. Like projection.
    // Also we can avoid repeated entries by only allowing fresh step numbers.
    if(storage_reason == StorageReason::PROJ_STATE) return;
    static size_t last_step_written = 0;
    if(status.step == last_step_written and last_step_written > 0) return;


//    tools::finite::io::h5table::write_measurements(*h5pp_file, table_prefix, storage_level, state, status);
    tools::finite::io::h5table::write_measurements(*h5pp_file, table_prefix, storage_level, tensors, status);
    tools::finite::io::h5table::write_sim_status(*h5pp_file, table_prefix, storage_level, status);
    tools::finite::io::h5table::write_profiling(*h5pp_file, table_prefix, storage_level, status);
    tools::finite::io::h5table::write_mem_usage(*h5pp_file, table_prefix, storage_level, status);
    last_step_written = status.step;
}

void class_algorithm_finite::copy_from_tmp(StorageReason storage_reason) {
    if(not h5pp_file) return;
    if(not settings::output::use_temp_dir) return;
    if(not tensors.state->position_is_any_edge()) return;
    switch(storage_reason) {
        case StorageReason::CHECKPOINT:
            if(math::mod(status.iter, settings::output::copy_from_temp_freq) != 0) return; // Check that we write according to the frequency given
        case StorageReason::FINISHED:
        case StorageReason::CHI_UPDATE:
        case StorageReason::PROJ_STATE:
        case StorageReason::INIT_STATE:
        case StorageReason::EMIN_STATE:
        case StorageReason::EMAX_STATE: break;
    }
    tools::common::io::h5tmp::copy_from_tmp(h5pp_file->getFilePath());
}

void class_algorithm_finite::print_status_update() {
    if(math::mod(status.step, print_freq()) != 0) return;
    if(print_freq() == 0) return;

    std::string report;
    report += fmt::format("{:<} ", algo_name);
    report += fmt::format("iter: {:<4} ", status.iter);
    report += fmt::format("step: {:<5} ", status.step);
    report += fmt::format("L: {} ", tensors.state->get_length());
    if(tensors.state->active_sites.empty())
        report += fmt::format("l: {:<2} ", tensors.state->get_position());
    else if(tensors.state->get_direction() > 0)
        report += fmt::format("l: [{:>2}-{:<2}] ", tensors.state->active_sites.front(), tensors.state->active_sites.back());
    else if(tensors.state->get_direction() < 0)
        report += fmt::format("l: [{:>2}-{:<2}] ", tensors.state->active_sites.back(), tensors.state->active_sites.front());
    report += fmt::format("E/L: {:<20.16f} ", tools::finite::measure::energy_per_site(*tensors.state, *tensors.model, *tensors.edges));
    if(algo_type == AlgorithmType::xDMRG) {
        report += fmt::format("ε: {:<6.4f} ", status.energy_dens);
    }
    report += fmt::format("Sₑ(l): {:<10.8f} ", tools::finite::measure::entanglement_entropy_current(*tensors.state));
    report += fmt::format("log₁₀ σ²(E)/L: {:<10.6f} [{:<10.6f}] ", std::log10(tools::finite::measure::energy_variance_per_site(*tensors.state, *tensors.model, *tensors.edges)),
                          std::log10(tensors.state->lowest_recorded_variance / static_cast<double>(tensors.state->get_length())));
    report += fmt::format("χmax: {:<3} χlim: {:<3} χ: {:<3} ", chi_max(), tensors.state->get_chi_lim(), tools::finite::measure::bond_dimension_current(*tensors.state));
    report += fmt::format("log₁₀ trunc: {:<10.4f} ", std::log10(tensors.state->get_truncation_error(tensors.state->get_position())));
    report += fmt::format("stk: {:<1} ", status.algorithm_has_stuck_for);
    report += fmt::format("sat: [σ² {:<1} Sₑ {:<1}] ", status.variance_mpo_saturated_for, status.entanglement_saturated_for);
    report += fmt::format("con: {:<5} ", status.algorithm_has_converged);
    report += fmt::format("time: {:<8.2f}s ", tools::common::profile::t_tot->get_age());
    report += fmt::format("mem MB: [Rss {:<.1f} Peak {:<.1f} Vm {:<.1f}] ", tools::common::profile::mem_rss_in_mb(), tools::common::profile::mem_hwm_in_mb(),
                          tools::common::profile::mem_vm_in_mb());
    tools::log->info(report);
}

void class_algorithm_finite::print_status_full() {
    tools::log->info("{:=^60}", "");
    tools::log->info("= {: ^56} =", "Final results [" + algo_name + "]");
    tools::log->info("{:=^60}", "");
    tools::log->info("--- Final results  --- {} ---", algo_name);
    tools::log->info("Sites                              = {}", tensors.state->get_length());
    tools::log->info("Iterations (sweeps)                = {}", status.iter);
    tools::log->info("Steps                              = {}", status.step);
    tools::log->info("Total time                         = {:<.1f} s = {:<.2f} min", tools::common::profile::t_tot->get_age(),
                     tools::common::profile::t_tot->get_age() / 60);
    tools::log->info("Energy per site E/L                = {:<.16f}", tools::finite::measure::energy_per_site(*tensors.state, *tensors.model, *tensors.edges));
    if(algo_type == AlgorithmType::xDMRG) {
        tools::log->info("Energy density (rescaled 0 to 1) ε = {:<6.4f}", tools::finite::measure::energy_normalized(*tensors.state, *tensors.model, *tensors.edges, status.energy_min, status.energy_max));
    }
    tools::log->info("Variance per site log₁₀ σ²(E)/L    = {:<.16f}", std::log10(tools::finite::measure::energy_variance_per_site(*tensors.state, *tensors.model, *tensors.edges)));
    tools::log->info("Bond dimension maximum χmax        = {}", chi_max());
    tools::log->info("Bond dimensions χ                  = {}", tools::finite::measure::bond_dimensions(*tensors.state));
    tools::log->info("Bond dimension  χ (mid)            = {}", tools::finite::measure::bond_dimension_midchain(*tensors.state));
    tools::log->info("Entanglement entropies Sₑ          = {}", tools::finite::measure::entanglement_entropies(*tensors.state));
    tools::log->info("Entanglement entropiy Sₑ (mid)     = {}", tools::finite::measure::entanglement_entropy_midchain(*tensors.state));
    tools::log->info("Truncation Errors                  = {}", tensors.state->get_truncation_errors());
    tools::log->info("Simulation converged               = {:<}", status.algorithm_has_converged);
    tools::log->info("Simulation saturated               = {:<}", status.algorithm_has_saturated);
    tools::log->info("Simulation succeeded               = {:<}", status.algorithm_has_succeeded);
    tools::log->info("Simulation got stuck               = {:<}", status.algorithm_has_got_stuck);
    tools::log->info("σ² slope                           = {:<8.4f} %   Converged : {:<8}  Saturated: {:<8}", V_mpo_slopes.back(),
                     status.variance_mpo_has_converged, status.variance_mpo_has_saturated);
    tools::log->info("Sₑ slope                           = {:<8.4f} %   Converged : {:<8}  Saturated: {:<8}", S_slopes.back(),
                     status.entanglement_has_converged, status.entanglement_has_saturated);
}
