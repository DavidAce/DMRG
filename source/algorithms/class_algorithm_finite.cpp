//
// Created by david on 2019-06-24.
//

#include "class_algorithm_finite.h"
#include <config/nmspc_settings.h>
#include <h5pp/h5pp.h>
#include <math/num.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/debug.h>
#include <tools/finite/env.h>
#include <tools/finite/io.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <tools/finite/print.h>

class_algorithm_finite::class_algorithm_finite(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type)
    : class_algorithm_base(std::move(h5ppFile_), algo_type) {
    tools::log->trace("Constructing class_algorithm_finite");
    tensors.initialize(settings::model::model_type, settings::model::model_size, 0);
    S_mat.resize(tensors.get_length() + 1);
    X_mat.resize(tensors.get_length() + 1);
    //    tools::finite::print::dimensions(*tensors.model);
    //    tools::finite::print::dimensions(tensors);
}

// We need to make a destructor manually for the enclosing class "class_model_finite"
// that encloses "class_model_base". Otherwise unique_ptr will forcibly inline its
// own default deleter.
// This allows us to forward declare the abstract base class "class_model_base"
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
// class_algorithm_finite::~class_algorithm_finite() = default;

void class_algorithm_finite::run()
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
{
    tools::log->info("Starting {}", algo_name);
    tools::common::profile::t_tot->tic();

    // We may want to resume this simulation.
    if(settings::output::file_collision_policy == FileCollisionPolicy::RESUME and h5pp_file->linkExists("common/storage_level")) {
        try {
            resume();
        } catch(std::exception &ex) {
            tools::log->info("Could not resume state from file [{}]: {}", h5pp_file->getFilePath(), ex.what());
            exit(1);
            //            run_default_task_list();
        }
    } else {
        run_default_task_list();
    }
    tools::common::profile::t_tot->toc();
}

void class_algorithm_finite::run_preprocessing() {
    tools::log->info("Running default preprocessing for {}", algo_name);
    tools::common::profile::prof[algo_type]["t_pre"]->tic();
    status.clear();
    randomize_model();
    init_bond_dimension_limits();
    randomize_state(ResetReason::INIT, settings::strategy::initial_state, settings::strategy::target_sector, settings::input::bitfield,
                    settings::strategy::use_eigenspinors);
    write_to_file(StorageReason::MODEL);
    tools::common::profile::prof[algo_type]["t_pre"]->toc();
    tools::log->info("Finished default preprocessing for {}", algo_name);
}

void class_algorithm_finite::run_postprocessing() {
    tools::log->info("Running default postprocessing for {}", algo_name);
    tools::common::profile::prof[algo_type]["t_pos"]->tic();
    write_to_file(StorageReason::CHECKPOINT, CopyPolicy::OFF);
    write_to_file(StorageReason::PROJ_STATE, CopyPolicy::OFF);
    write_to_file(StorageReason::FINISHED, CopyPolicy::FORCE);
    print_status_full();
    tools::common::profile::prof[algo_type]["t_pos"]->toc();
    tools::log->info("Finished default postprocessing for {}", algo_name);
}

void class_algorithm_finite::move_center_point(std::optional<size_t> num_moves) {
    if(not num_moves.has_value()) {
        if(tensors.active_sites.empty()) num_moves = 1ul;
        else if(settings::strategy::multisite_move == MultisiteMove::ONE)
            num_moves = 1ul;
        else if(settings::strategy::multisite_move == MultisiteMove::MID)
            num_moves = std::max<size_t>(1, (tensors.active_sites.size()) / 2);
        else if(settings::strategy::multisite_move == MultisiteMove::MAX) {
            size_t offset = 1ul;
            if(tensors.state->get_direction() == 1 and tensors.active_sites.back() == tensors.get_length() - 1) offset = 2ul;
            if(tensors.state->get_direction() == -1 and tensors.active_sites.front() == 0) offset = 2ul;
            num_moves = std::max<size_t>(1, tensors.active_sites.size() - offset);

        } else
            throw std::logic_error("Could not determine how many sites to move");
    }

    tools::log->debug("Moving center point {} steps in direction {}", num_moves.value(), tensors.state->get_direction());
    tensors.clear_cache();
    tensors.clear_measurements();
    try {
        for(size_t i = 0; i < num_moves.value(); i++) {
            if(tensors.position_is_any_edge()) {
                tensors.state->increment_iter();
                has_projected = false;
            }
            tensors.move_center_point(status.chi_lim);
            tensors.state->increment_step();
            if(chi_quench_steps > 0) chi_quench_steps--;
        }
        tools::finite::debug::check_integrity(*tensors.state, *tensors.model, *tensors.edges);
        tensors.active_sites.clear();
        tensors.state->active_sites.clear();
        tensors.model->active_sites.clear();
        tensors.edges->active_sites.clear();
    } catch(std::exception &e) {
        tools::finite::print::dimensions(tensors);
        throw std::runtime_error("Failed to move center point: " + std::string(e.what()));
    }
    status.iter      = tensors.state->get_iteration();
    status.step      = tensors.state->get_step();
    status.position  = tensors.state->get_position();
    status.direction = tensors.state->get_direction();
    has_projected    = false;
}

void class_algorithm_finite::reduce_mpo_energy() {
    // Reduce mpo energy to avoid catastrophic cancellation
    // Note that this operation makes the Hamiltonian nearly singular,
    // which is tough for Lanczos/Arnoldi iterations to handle
    if(settings::precision::use_reduced_energy and tensors.position_is_any_edge()) {
        tensors.reduce_mpo_energy();
    }
}

void class_algorithm_finite::update_bond_dimension_limit([[maybe_unused]] std::optional<long> tmp_bond_limit) {
    if(not tensors.position_is_any_edge()) return;
    status.chi_lim_has_reached_chi_max = status.chi_lim >= status.chi_lim_max;
    if(not cfg_chi_lim_grow()) {
        status.chi_lim = status.chi_lim_max;
        return;
    }
    if(status.chi_lim_has_reached_chi_max) return;

    // If we got here we want increase the bond dimension limit progressively during the simulation
    // Only increment the bond dimension if the following are all true
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
    if(not status.algorithm_has_got_stuck and status.chi_lim >= 16) {
        tools::log->info("State is not stuck yet. Kept current limit {}", status.chi_lim);
        return;
    }

    if(tools::log->level() < spdlog::level::info) {
        double truncation_threshold = 2 * settings::precision::svd_threshold;
        size_t trunc_bond_count     = tensors.state->num_sites_truncated(truncation_threshold);
        size_t bond_at_lim_count    = tensors.state->num_bonds_reached_chi(status.chi_lim);
        tools::log->info("Truncation threshold  : {:<.8e}", truncation_threshold);
        tools::log->info("Truncation errors     : {}", tensors.state->get_truncation_errors());
        tools::log->info("Bond dimensions       : {}", tools::finite::measure::bond_dimensions(*tensors.state));
        tools::log->info("Truncated bond count  : {} ", trunc_bond_count);
        tools::log->info("Bonds at limit  count : {} ", bond_at_lim_count);
        tools::log->info("Entanglement entropies: {} ", tools::finite::measure::entanglement_entropies(*tensors.state));
    }
    bool state_is_bond_limited = tensors.state->is_bond_limited(status.chi_lim, 5 * settings::precision::svd_threshold);
    if(not state_is_bond_limited) {
        tools::log->info("State is not limited by its bond dimension. Kept current limit {}", status.chi_lim);
        return;
    }

    // Write current results before updating bond dimension
    write_to_file(StorageReason::CHI_UPDATE);
    if(settings::strategy::randomize_on_chi_update and status.chi_lim >= 32)
        randomize_state(ResetReason::CHI_UPDATE, StateType::RANDOMIZE_PREVIOUS_STATE, std::nullopt, status.chi_lim);

    tools::log->info("Updating bond dimension limit {} -> {}", status.chi_lim, status.chi_lim * 2);
    status.chi_lim *= 2;
    clear_convergence_status();
    status.chi_lim_has_reached_chi_max = status.chi_lim == status.chi_lim_max;

    // Last sanity check before leaving here
    if(status.chi_lim > status.chi_lim_max)
        throw std::runtime_error(fmt::format("chi_lim is larger than cfg_chi_lim_max! {} > {}", status.chi_lim, status.chi_lim_max));
}

void class_algorithm_finite::randomize_model() {
    tools::log->info("Randomizing model");
    tensors.randomize_model();
    clear_convergence_status();
}

void class_algorithm_finite::randomize_state(ResetReason reason, StateType state_type, std::optional<std::string> sector, std::optional<long> chi_lim,
                                             std::optional<bool> use_eigenspinors, std::optional<long> bitfield, std::optional<double> svd_threshold) {
    tools::log->info("Randomizing state [{}] to type [{}] | Reason {} ...", state_name, enum2str(state_type), enum2str(reason));
    if(reason == ResetReason::SATURATED) {
        if(status.num_resets >= settings::strategy::max_resets)
            return tools::log->warn("Skipped reset: num resets {} >= max resets {}", status.num_resets, settings::strategy::max_resets);
        else
            status.num_resets++; // Only increment if doing it for saturation reasons
    }
    tools::common::profile::prof[algo_type]["t_rnd"]->tic();
    if(not sector) sector = settings::strategy::target_sector;
    if(not chi_lim) {
        if(state_type == StateType::RANDOMIZE_PREVIOUS_STATE)
            chi_lim = static_cast<long>(std::pow(2, std::floor(std::log2(tensors.state->find_largest_chi())))); // Nearest power of two from below
        else
            chi_lim = cfg_chi_lim_init();
    }
    if(not use_eigenspinors) use_eigenspinors = settings::strategy::use_eigenspinors;
    if(not bitfield) bitfield = settings::input::bitfield;
    if(not svd_threshold and state_type == StateType::RANDOMIZE_PREVIOUS_STATE) svd_threshold = 1e-4;

    tensors.randomize_state(state_type, sector.value(), chi_lim.value(), use_eigenspinors.value(), bitfield, svd_threshold);
    while(tensors.state->get_position() != 0 and tensors.state->get_direction() != 1) move_center_point();
    clear_convergence_status();
    status.reset();
    status.iter      = tensors.state->reset_iter();
    status.step      = tensors.state->reset_step();
    status.position  = tensors.state->get_position();
    status.direction = tensors.state->get_direction();
    if(cfg_chi_lim_grow()) status.chi_lim = chi_lim.value();
    if(reason == ResetReason::NEW_STATE) excited_state_number++;

    if(tensors.state->find_largest_chi() > chi_lim.value())
        throw std::runtime_error(
            fmt::format("Faulty truncation after randomize. Max found chi is {}, but chi limit is {}", tensors.state->find_largest_chi(), chi_lim.value()));

    tools::log->info("Randomizing state [{}] to type [{}] | Reason {} ... OK!", state_name, enum2str(state_type), enum2str(reason));
    tools::log->info("Spin components {}", tools::finite::measure::spin_components(*tensors.state));
    tools::log->info("Bond dimensions {}", tools::finite::measure::bond_dimensions(*tensors.state));
    tools::common::profile::prof[algo_type]["t_rnd"]->toc();
}

void class_algorithm_finite::try_projection() {
    if(not tensors.position_is_any_edge()) return;
    if(has_projected) return;
    if(settings::strategy::project_on_every_sweep or (settings::strategy::project_when_stuck and status.algorithm_has_got_stuck)) {
        tools::log->info("Trying projection to {}", settings::strategy::target_sector);
        tensors.project_to_nearest_sector(settings::strategy::target_sector);
        has_projected = true;
        write_to_file(StorageReason::PROJ_STATE);
    }
}

void class_algorithm_finite::try_bond_dimension_quench() {
    if(not settings::strategy::chi_quench_when_stuck) return;
    if(chi_quench_steps > 0) clear_convergence_status();
    if(not tensors.position_is_any_edge()) return;
    if(chi_quench_steps >= tensors.get_length()) {
        tools::log->info("Chi quench continues -- {} steps left", chi_quench_steps);
        tools::finite::mps::truncate_all_sites(*tensors.state, chi_lim_quench_ahead);
        return;
    }
    if(status.algorithm_has_stuck_for <= 2) {
        tools::log->info("Chi quench skipped: simulation not been stuck for long enough");
        return;
    }
    if(num_chi_quenches >= max_chi_quenches) {
        tools::log->trace("Chi quench skipped: max number of chi quenches ({}) have been made already", num_chi_quenches);
        return;
    }
    if(tools::finite::measure::energy_variance(tensors) < 10 * settings::precision::variance_convergence_threshold) return;
    double truncation_threshold = 5 * settings::precision::svd_threshold;
    size_t trunc_bond_count     = tensors.state->num_sites_truncated(truncation_threshold);
    size_t bond_at_lim_count    = tensors.state->num_bonds_reached_chi(status.chi_lim);
    tools::log->trace("Truncation threshold : {:.4e}", std::pow(truncation_threshold, 2));
    tools::log->trace("Truncation errors    : {}", tensors.state->get_truncation_errors());
    tools::log->trace("Bond dimensions      : {}", tools::finite::measure::bond_dimensions(*tensors.state));
    tools::log->trace("Entanglement entr    : {}", tools::finite::measure::entanglement_entropies(*tensors.state));
    tools::log->trace("Truncated bond count : {} ", trunc_bond_count);
    tools::log->trace("Bonds at limit  count: {} ", bond_at_lim_count);
    if(tensors.state->is_bond_limited(status.chi_lim, truncation_threshold)) {
        tools::log->info("Chi quench skipped: state is bond limited - prefer updating bond dimension");
        return;
    }
    auto bond_dimensions    = tools::finite::measure::bond_dimensions(*tensors.state);
    auto max_bond_dimension = *max_element(std::begin(bond_dimensions), std::end(bond_dimensions));
    tools::log->info("Chi quench started");
    tools::finite::mps::truncate_all_sites(*tensors.state, max_bond_dimension / 2);
    tools::log->debug("Bond dimensions: {}", tools::finite::measure::bond_dimensions(*tensors.state));
    clear_convergence_status();
    chi_quench_steps = 1 * tensors.get_length();
    num_chi_quenches++;
}

void class_algorithm_finite::try_hamiltonian_perturbation() {
    if(not settings::strategy::perturb_when_stuck) return;
    if(perturbation_steps-- > 0) return;
    if(not status.algorithm_has_got_stuck) {
        tools::log->info("Perturbation skipped: simulation not stuck");
        return;
    }
    if(num_perturbations >= max_perturbations) {
        tools::log->info("Perturbation skipped: max number of perturbation trials ({}) have been made already", num_perturbations);
        return;
    }
    if(tensors.model->is_perturbed()) {
        tensors.perturb_model_params(0, 0, PerturbMode::UNIFORM_RANDOM_PERCENTAGE);
    } else {
        tensors.perturb_model_params(1e-2, 1e-2, PerturbMode::UNIFORM_RANDOM_PERCENTAGE);
        perturbation_steps = 2 * (tensors.get_length() - 1);
    }
}

void class_algorithm_finite::try_disorder_damping() {
    // if(not state.position_is_left_edge()) return;
    // If there are damping exponents to process, do so
    if(not damping_exponents.empty()) {
        tools::log->info("Setting damping exponent = {}", damping_exponents.back());
        tensors.damp_model_disorder(damping_exponents.back(), 0);
        damping_exponents.pop_back();
        if(damping_exponents.empty() and tensors.model->is_damped()) throw std::logic_error("Damping trial ended but the state is still damped");
    }
    if(not settings::strategy::damping_when_stuck) return;
    if(tensors.model->is_damped()) return;
    if(not status.algorithm_has_got_stuck) {
        tools::log->info("Damping skipped: simulation not stuck");
        return;
    }
    if(num_dampings >= max_dampings) {
        tools::log->info("Damping skipped: max number of damping trials ({}) have been made already", num_dampings);
        return;
    }
    // damping_exponents = num::LogSpaced(6,0.0,0.25,0.001);
    // damping_exponents = {0.0,0.1};
    damping_exponents = {0.0, 0.01, 0.02, 0.04, 0.08, 0.16, 0.256};
    tools::log->info("Generating damping exponents = {}", damping_exponents);
    has_damped = true;
    clear_convergence_status();
    tools::log->info("Setting damping exponent = {}", damping_exponents.back());
    tensors.damp_model_disorder(damping_exponents.back(), 0);
    damping_exponents.pop_back();
    if(damping_exponents.empty() and tensors.model->is_damped()) throw std::logic_error("Damping trial ended but the state is still damped");
}

void class_algorithm_finite::check_convergence_variance(std::optional<double> threshold, std::optional<double> slope_threshold) {
    // Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    //    if(not tensors.state->position_is_any_edge()) {
    //        return;
    //    }
    tools::log->trace("Checking convergence of variance mpo");
    if(not threshold) threshold = settings::precision::variance_convergence_threshold;
    if(not slope_threshold) slope_threshold = settings::precision::variance_slope_threshold;
    auto report = check_saturation_using_slope(V_mpo_vec, X_mpo_vec, tools::finite::measure::energy_variance(tensors), status.iter, 1, slope_threshold.value());
    status.variance_mpo_has_converged = tools::finite::measure::energy_variance(tensors) < threshold;
    if(report.has_computed) {
        V_mpo_slopes.emplace_back(report.slope);
        auto last_nonconverged_ptr = std::find_if(V_mpo_vec.rbegin(), V_mpo_vec.rend(), [threshold](auto const &val) { return val > threshold; });
        auto last_nonsaturated_ptr =
            std::find_if(V_mpo_slopes.rbegin(), V_mpo_slopes.rend(), [slope_threshold](auto const &val) { return val > slope_threshold; });
        auto converged_count              = static_cast<size_t>(std::distance(V_mpo_vec.rbegin(), last_nonconverged_ptr));
        auto saturated_count              = static_cast<size_t>(std::distance(V_mpo_slopes.rbegin(), last_nonsaturated_ptr));
        status.variance_mpo_has_saturated = report.slope < slope_threshold; // or saturated_count >= min_saturation_iters;
        status.variance_mpo_saturated_for = std::max(converged_count, saturated_count);
        std::vector<double> V_mpo_vec_log10;
        for(const auto &v : V_mpo_vec) V_mpo_vec_log10.emplace_back(std::log10(v));
        if(tools::log->level() == spdlog::level::debug)
            tools::log->debug("Energy variance convergence:  rel. slope {:.3f} | slope tolerance {:.2f}", report.slope, slope_threshold.value());
        else if(tools::log->level() == spdlog::level::trace) {
            tools::log->trace("Energy variance slope details:");
            tools::log->trace(" -- relative slope     = {:.3f} %", report.slope);
            tools::log->trace(" -- slope tolerance    = {:.3f} %", slope_threshold.value());
            tools::log->trace(" -- avgerage from iter = {} ", report.check_from);
            tools::log->trace(" -- log10 var average  = {:.3f} ", std::log10(report.avgY));
            tools::log->trace(" -- log10 var history  = {:.3f} ", fmt::join(V_mpo_vec_log10, ", "));
            tools::log->trace(" -- slope history      = {:.3f} ", fmt::join(V_mpo_slopes, ", "));
            tools::log->trace(" -- has saturated      = {} for {} iters ", status.variance_mpo_has_saturated, status.variance_mpo_saturated_for);
            tools::log->trace(" -- has converged      = {} for {} iters ", status.variance_mpo_has_converged, converged_count);
        }
        if(V_mpo_vec.back() < threshold and status.variance_mpo_saturated_for == 0) throw std::logic_error("Variance should have saturated");
        if(V_mpo_vec.back() < threshold and not status.variance_mpo_has_converged) throw std::logic_error("Variance should have converged");
    }
}

void class_algorithm_finite::check_convergence_entg_entropy(std::optional<double> slope_threshold) {
    // Based on the the slope of entanglement entanglement_entropy_midchain
    // This one is cheap to compute.
    if(not tensors.state->position_is_any_edge()) { return; }
    tools::log->trace("Checking convergence of entanglement");
    if(not slope_threshold) slope_threshold = settings::precision::entropy_slope_threshold;
    auto                          entropies = tools::finite::measure::entanglement_entropies(*tensors.state);
    std::vector<SaturationReport> reports(entropies.size());

    for(size_t site = 0; site < entropies.size(); site++)
        reports[site] = check_saturation_using_slope(S_mat[site], X_mat[site], entropies[site], status.iter, 1, slope_threshold.value());

    bool all_computed                 = std::all_of(reports.begin(), reports.end(), [](const SaturationReport r) { return r.has_computed; });
    status.entanglement_has_saturated = false;
    if(all_computed) {
        // idx_max_slope is the index to the site with maximum slope
        auto idx_max_slope = (size_t) std::distance(
            reports.begin(),
            std::max_element(reports.begin(), reports.end(), [](const SaturationReport &r1, const SaturationReport &r2) { return r1.slope < r2.slope; }));

        S_slopes.push_back(reports[idx_max_slope].slope);
        auto last_nonsaturated_ptr = std::find_if(S_slopes.rbegin(), S_slopes.rend(), [slope_threshold](auto const &val) { return val > slope_threshold; });
        auto saturated_count       = static_cast<size_t>(std::distance(S_slopes.rbegin(), last_nonsaturated_ptr));

        status.entanglement_has_saturated = S_slopes.back() < slope_threshold;
        status.entanglement_saturated_for = saturated_count;
        std::vector<double> all_avergs;
        std::vector<double> all_slopes;
        for(auto &r : reports) all_avergs.push_back(r.avgY);
        for(auto &r : reports) all_slopes.push_back(r.slope);
        if(tools::log->level() == spdlog::level::debug)
            tools::log->debug("Entanglement ent. convergence at site {}: rel. slope {:.8f} | slope tolerance {:.2f}", idx_max_slope,
                              reports[idx_max_slope].slope, slope_threshold.value());
        else if(tools::log->level() == spdlog::level::trace) {
            tools::log->trace("Entanglement slope details:");
            tools::log->trace(" -- site              = {}", idx_max_slope);
            tools::log->trace(" -- relative slope    = {:.8f} %", reports[idx_max_slope].slope);
            tools::log->trace(" -- slope tolerance   = {:.8f} %", slope_threshold.value());
            tools::log->trace(" -- average from iter = {} ", reports[idx_max_slope].check_from);
            tools::log->trace(" -- ent history       = {:.6f} ", fmt::join(S_mat[idx_max_slope], ", "));
            tools::log->trace(" -- slope history     = {:.3f} ", fmt::join(S_slopes, ", "));
            tools::log->trace(" -- has saturated     = {} for {} iters (site {})", status.entanglement_has_saturated, status.entanglement_saturated_for,
                              idx_max_slope);
            tools::log->trace(" -- all averages      = {:.6f} ", fmt::join(all_avergs, ", "));
            tools::log->trace(" -- all slopes        = {:.3f} ", fmt::join(all_slopes, ", "));
        }
        if(reports[idx_max_slope].slope > slope_threshold and status.entanglement_has_saturated) throw std::logic_error("Not supposed to be saturated!!");
    }
    status.entanglement_has_converged = status.entanglement_has_saturated;
}

void class_algorithm_finite::clear_convergence_status() {
    tools::log->trace("Clearing convergence status");
    for(auto &mat : S_mat) { mat.clear(); }
    for(auto &mat : X_mat) { mat.clear(); }
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
    status.algorithm_has_to_stop       = false;
    status.algorithm_has_got_stuck     = false;
    status.algorithm_has_converged     = false;
    status.algorithm_has_saturated     = false;
    status.algorithm_has_succeeded     = false;
    status.algorithm_has_stuck_for     = 0;
    has_projected                      = false;
    has_damped                         = false;
}

void class_algorithm_finite::write_to_file(StorageReason storage_reason, std::optional<CopyPolicy> copy_policy) {
    write_to_file(storage_reason, *tensors.state, copy_policy);
}

void class_algorithm_finite::write_to_file(StorageReason storage_reason, const class_state_finite &state, std::optional<CopyPolicy> copy_policy,
                                           bool is_projection, const std::string &given_prefix) {
    // Setup this save
    StorageLevel             storage_level;
    std::string              state_prefix = algo_name + '/' + state_name; // May get modified
    std::string              model_prefix = algo_name + "/model";
    std::vector<std::string> table_prefxs = {algo_name + '/' + state_name + "/tables"}; // Common tables
    if(not given_prefix.empty()) state_prefix = given_prefix;
    switch(storage_reason) {
        case StorageReason::FINISHED: {
            if(status.algorithm_has_succeeded) storage_level = settings::output::storage_level_good_state;
            else
                storage_level = settings::output::storage_level_fail_state;
            state_prefix += "/finished";
            table_prefxs.emplace_back(state_prefix); // Appends to its own table as well as the common ones
            break;
        }
        case StorageReason::CHECKPOINT: {
            if(not state.position_is_any_edge()) return;
            if(num::mod(status.iter, settings::output::checkpoint_frequency) != 0) return;
            state_prefix += "/checkpoint";
            storage_level = settings::output::storage_level_checkpoint;
            if(settings::output::checkpoint_keep_newest_only) state_prefix += "/iter_last";
            else
                state_prefix += fmt::format("/iter_{}", status.iter);
            table_prefxs.emplace_back(state_prefix); // Appends to its own table as well as the common ones
            break;
        }
        case StorageReason::CHI_UPDATE: {
            if(not cfg_chi_lim_grow()) return;
            // If we have updated chi we may want to write a projection too
            storage_level = settings::output::storage_level_checkpoint;
            state_prefix += "/checkpoint";
            state_prefix += fmt::format("/chi_{}", status.chi_lim);
            table_prefxs = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::PROJ_STATE: {
            storage_level = settings::output::storage_level_proj_state;
            if(not is_projection and storage_level != StorageLevel::NONE) {
                auto state_projected = tools::finite::ops::get_projection_to_nearest_sector(*tensors.state, settings::strategy::target_sector);
                return write_to_file(storage_reason, state_projected, copy_policy, true, state_prefix);
            }
            state_prefix += "/projection";
            table_prefxs = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::INIT_STATE: {
            storage_level = settings::output::storage_level_init_state;
            state_prefix += "/state_init";
            table_prefxs = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::EMIN_STATE: {
            storage_level = settings::output::storage_level_emin_state;
            state_prefix  = algo_name + "/state_emin";
            table_prefxs  = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::EMAX_STATE: {
            storage_level = settings::output::storage_level_emax_state;
            state_prefix  = algo_name + "/state_emax";
            table_prefxs  = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::MODEL: {
            storage_level = settings::output::storage_level_model;
            tools::finite::io::h5table::save_model(*h5pp_file, model_prefix + "/hamiltonian", storage_level, *tensors.model);
            tools::finite::io::h5dset::save_model(*h5pp_file, model_prefix + "/mpo", storage_level, *tensors.model);
            return copy_from_tmp(storage_reason, CopyPolicy::TRY);
        }
    }
    if(storage_level == StorageLevel::NONE) return;
    if(state_prefix.empty()) throw std::runtime_error("State prefix is empty");
    tools::log->info("Writing to file: Reason [{}] | Level [{}] | hdf5 prefix [{}]", enum2str(storage_reason), enum2str(storage_level), state_prefix);
    // Start saving tensors and metadata
    tools::finite::io::h5dset::save_state(*h5pp_file, state_prefix, storage_level, state, status);
    tools::finite::io::h5dset::save_entgm(*h5pp_file, state_prefix, storage_level, state, status);
    tools::common::io::h5attr::save_meta(*h5pp_file, storage_level, storage_reason, settings::model::model_type, settings::model::model_size, algo_type,
                                         state_name, state_prefix, model_prefix, status);

    // The main results have now been written. Next we append data to tables
    for(const auto &table_prefix : table_prefxs) {
        tools::finite::io::h5table::save_measurements(*h5pp_file, table_prefix + "/measurements", storage_level, tensors, status);
        tools::finite::io::h5table::save_sim_status(*h5pp_file, table_prefix + "/status", storage_level, status);
        tools::finite::io::h5table::save_profiling(*h5pp_file, table_prefix + "/profiling", storage_level, status);
        tools::finite::io::h5table::save_mem_usage(*h5pp_file, table_prefix + "/mem_usage", storage_level, status);
    }
    // Copy from temporary location to destination depending on given policy
    copy_from_tmp(storage_reason, copy_policy);
}

void class_algorithm_finite::print_status_update() {
    if(num::mod(status.step, cfg_print_freq()) != 0) return;
    if(cfg_print_freq() == 0) return;

    std::string report;
    //    report += fmt::format("{:<} ", algo_name);
    report += fmt::format("{:<} ", state_name);
    report += fmt::format("iter: {:<4} ", status.iter);
    report += fmt::format("step: {:<5} ", status.step);
    report += fmt::format("L: {} ", tensors.state->get_length());
    if(tensors.state->active_sites.empty()) report += fmt::format("l: {:<2} ", tensors.state->get_position());
    else if(tensors.state->get_direction() > 0)
        report += fmt::format("l: [{:>2}-{:<2}] ", tensors.state->active_sites.front(), tensors.state->active_sites.back());
    else if(tensors.state->get_direction() < 0)
        report += fmt::format("l: [{:>2}-{:<2}] ", tensors.state->active_sites.back(), tensors.state->active_sites.front());
    report += fmt::format("E/L: {:<20.16f} ", tools::finite::measure::energy_per_site(tensors));
    if(algo_type == AlgorithmType::xDMRG) { report += fmt::format("ε: {:<6.4f} ", status.energy_dens); }
    report += fmt::format("Sₑ(l): {:<10.8f} ", tools::finite::measure::entanglement_entropy_current(*tensors.state));
    report += fmt::format("log₁₀ σ²(E)/L: {:<10.6f} [{:<10.6f}] ", std::log10(tools::finite::measure::energy_variance_per_site(tensors)),
                          std::log10(status.lowest_recorded_variance_per_site));
    report +=
        fmt::format("χmax: {:<3} χlim: {:<3} χ: {:<3} ", cfg_chi_lim_max(), status.chi_lim, tools::finite::measure::bond_dimension_current(*tensors.state));
    report += fmt::format("log₁₀ trunc: {:<10.4f} ", std::log10(tensors.state->get_truncation_error(tensors.state->get_position())));
    report += fmt::format("stk: {:<1} ", status.algorithm_has_stuck_for);
    report += fmt::format("sat: [σ² {:<1} Sₑ {:<1}] ", status.variance_mpo_saturated_for, status.entanglement_saturated_for);
    report += fmt::format("con: {:<5} ", status.algorithm_has_converged);
    report += fmt::format("time:{:>8.2f}s ", tools::common::profile::t_tot->get_measured_time());
    report += fmt::format("mem: [rss {:<.1f} peak {:<.1f} vm {:<.1f}] MB ", tools::common::profile::mem_rss_in_mb(), tools::common::profile::mem_hwm_in_mb(),
                          tools::common::profile::mem_vm_in_mb());
    tools::log->info(report);
}

void class_algorithm_finite::print_status_full() {
    tensors.do_all_measurements();
    tools::log->info("{:=^60}", "");
    tools::log->info("= {: ^56} =", "Final results [" + algo_name + "][" + state_name + "]");
    tools::log->info("{:=^60}", "");
    tools::log->info("Stop reason                        = {}", enum2str(stop_reason));
    tools::log->info("Sites                              = {}", tensors.state->get_length());
    tools::log->info("Position                           = {}", tensors.state->get_position());
    tools::log->info("Direction                          = {}", tensors.state->get_direction());
    tools::log->info("Iterations (full chain sweeps)     = {}", status.iter);
    tools::log->info("Steps (moves along the chain)      = {}", status.step);
    tools::log->info("Total time                         = {:<.1f} s = {:<.2f} min", tools::common::profile::t_tot->get_measured_time(),
                     tools::common::profile::t_tot->get_measured_time() / 60);
    tools::log->info("Energy per site E/L                = {:<.16f}", tools::finite::measure::energy_per_site(tensors));
    if(algo_type == AlgorithmType::xDMRG) {
        tools::log->info("Energy density (rescaled 0 to 1) ε = {:<6.4f}",
                         tools::finite::measure::energy_normalized(tensors, status.energy_min_per_site, status.energy_max_per_site));
    }
    tools::log->info("Variance per site log₁₀ σ²(E)/L    = {:<.16f}", std::log10(tools::finite::measure::energy_variance_per_site(tensors)));
    tools::log->info("Bond dimension maximum χmax        = {}", cfg_chi_lim_max());
    tools::log->info("Bond dimensions χ                  = {}", tools::finite::measure::bond_dimensions(*tensors.state));
    tools::log->info("Bond dimension  χ (mid)            = {}", tools::finite::measure::bond_dimension_midchain(*tensors.state));
    tools::log->info("Entanglement entropies Sₑ          = {:.6f}", fmt::join(tools::finite::measure::entanglement_entropies(*tensors.state), ", "));
    tools::log->info("Entanglement entropiy Sₑ (mid)     = {:.5f}", tools::finite::measure::entanglement_entropy_midchain(*tensors.state), ", ");
    tools::log->info("Truncation Errors                  = {:.3e}", fmt::join(tensors.state->get_truncation_errors(), ", "));
    tools::log->info("Algorithm has converged            = {:<}", status.algorithm_has_converged);
    tools::log->info("Algorithm has saturated            = {:<}", status.algorithm_has_saturated);
    tools::log->info("Algorithm has succeeded            = {:<}", status.algorithm_has_succeeded);
    tools::log->info("Algorithm has got stuck            = {:<}", status.algorithm_has_got_stuck);
    tools::log->info("σ²                                 = Converged : {:<8}  Saturated: {:<8}", status.variance_mpo_has_converged,
                     status.variance_mpo_has_saturated);
    tools::log->info("Sₑ                                 = Converged : {:<8}  Saturated: {:<8}", status.entanglement_has_converged,
                     status.entanglement_has_saturated);
    tools::log->info("{:=^60}", "");
}
