//
// Created by david on 2019-06-24.
//

#include "class_algorithm_finite.h"
#include <config/nmspc_settings.h>
#include <general/nmspc_exceptions.h>
#include <general/nmspc_iter.h>
#include <h5pp/h5pp.h>
#include <math/num.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
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
    tensors.state->set_algorithm(algo_type);
    entropy_iter.resize(tensors.get_length() + 1);

    max_stuck_iters      = settings::precision::max_stuck_iters;
    min_saturation_iters = settings::precision::min_saturation_iters;
    max_saturation_iters = settings::precision::max_saturation_iters;
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
    auto t_tot = tools::common::profile::t_tot->tic_token();

    // We may want to resume this simulation.
    if(settings::output::file_collision_policy == FileCollisionPolicy::RESUME and h5pp_file->linkExists("common/storage_level")) {
        try {
            resume();
        } catch(const except::state_error &ex) {
            throw except::resume_error(fmt::format("Could not resume state from file [{}]: {}", h5pp_file->getFilePath(), ex.what()));
        } catch(const except::load_error &ex) {
            throw except::load_error(fmt::format("Could not resume state from file [{}]: {}", h5pp_file->getFilePath(), ex.what()));
        } catch(const std::exception &ex) {
            throw std::runtime_error(fmt::format("Could not resume state from file [{}]: {}", h5pp_file->getFilePath(), ex.what()));
        }
    } else {
        run_default_task_list();
    }
}

void class_algorithm_finite::run_postprocessing() {
    tools::log->info("Running default postprocessing for {}", algo_name);
    auto t_pos = tools::common::profile::prof[algo_type]["t_pos"]->tic_token();
    write_to_file(StorageReason::CHECKPOINT, CopyPolicy::OFF);
    write_to_file(StorageReason::PROJ_STATE, CopyPolicy::OFF);
    write_to_file(StorageReason::FINISHED, CopyPolicy::FORCE);
    print_status_full();
    tools::log->info("Finished default postprocessing for {}", algo_name);
    tools::common::profile::print_profiling_all();
}

void class_algorithm_finite::move_center_point(std::optional<long> num_moves) {
    if(not num_moves.has_value()) {
        if(tensors.active_sites.empty())
            num_moves = 1;
        else {
            long num_sites = static_cast<long>(tensors.active_sites.size());
            if((tensors.state->get_direction() == 1 and tensors.active_sites.back() == tensors.get_length() - 1) or
               (tensors.state->get_direction() == -1 and tensors.active_sites.front() == 0)) {
                // In this case we have just updated from here to the edge. No point in updating
                // closer and closer to the edge. Just move until reaching the edge without flip
                num_moves = std::max<long>(1, num_sites - 1); // to the edge without flipping
            } else if(settings::strategy::multisite_mps_step == MultisiteMove::ONE)
                num_moves = 1ul;
            else if(settings::strategy::multisite_mps_step == MultisiteMove::MID) {
                num_moves = std::max<long>(1, num_sites / 2);
            } else if(settings::strategy::multisite_mps_step == MultisiteMove::MAX) {
                num_moves = std::max<long>(1, num_sites - 1); // Move so that the center point moves out of the active region
            } else
                throw std::logic_error("Could not determine how many sites to move");
        }
    }

    if(num_moves <= 0) throw std::runtime_error(fmt::format("Cannot move center point {} sites", num_moves.value()));
    tools::log->debug("Moving center point {} steps in direction {}", num_moves.value(), tensors.state->get_direction());
    tensors.clear_cache();
    tensors.clear_measurements();
    try {
        long moves = 0;
        while(num_moves > moves++) {
            if(chi_quench_steps > 0) chi_quench_steps--;
            if(tensors.position_is_outward_edge()) status.iter++;
            status.step += tensors.move_center_point(status.chi_lim);
            // Do not go past the edge if you aren't there already!
            // It's important to stay at the inward edge so we can do convergence checks and so on
            if(tensors.position_is_inward_edge()) break;
        }
        tensors.active_sites.clear();
        tensors.state->active_sites.clear();
        tensors.model->active_sites.clear();
        tensors.edges->active_sites.clear();
    } catch(std::exception &e) {
        tools::finite::print::dimensions(tensors);
        throw std::runtime_error("Failed to move center point: " + std::string(e.what()));
    }
    status.position  = tensors.state->get_position<long>();
    status.direction = tensors.state->get_direction();
    has_projected    = false;
}

void class_algorithm_finite::reduce_mpo_energy() {
    // Reduce mpo energy to avoid catastrophic cancellation
    // Note that this operation makes the Hamiltonian nearly singular,
    // which is tough for Lanczos/Arnoldi iterations to handle
    if(settings::precision::use_reduced_energy and tensors.position_is_inward_edge()) tensors.reduce_mpo_energy();
}

void class_algorithm_finite::update_bond_dimension_limit([[maybe_unused]] std::optional<long> tmp_bond_limit) {
    if(not tensors.position_is_inward_edge()) return;
    status.chi_lim_max                 = cfg_chi_lim_max();
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
    if(status.algorithm_has_stuck_for == 0 and status.chi_lim >= 16) {
        tools::log->info("Algorithm is not stuck yet. Kept current limit {}", status.chi_lim);
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
    bool state_is_bond_limited = tensors.state->is_bond_limited(status.chi_lim, 2 * settings::precision::svd_threshold);
    if(not state_is_bond_limited) {
        tools::log->info("State is not limited by its bond dimension. Kept current limit {}", status.chi_lim);
        return;
    }

    // Write current results before updating bond dimension
    write_to_file(StorageReason::CHI_UPDATE);
    if(settings::strategy::randomize_on_chi_update and status.chi_lim >= 32)
        randomize_state(ResetReason::CHI_UPDATE, StateInit::RANDOMIZE_PREVIOUS_STATE, std::nullopt, std::nullopt, status.chi_lim);

    tools::log->info("Updating bond dimension limit {} -> {}", status.chi_lim, status.chi_lim * 2);
    status.chi_lim *= 2;
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

void class_algorithm_finite::randomize_state(ResetReason reason, StateInit state_init, std::optional<StateInitType> state_type,
                                             std::optional<std::string> sector, std::optional<long> chi_lim, std::optional<bool> use_eigenspinors,
                                             std::optional<long> bitfield, std::optional<double> svd_threshold) {
    tools::log->info("Randomizing state [{}] to [{}] | Reason [{}] ...", tensors.state->get_name(), enum2str(state_init), enum2str(reason));
    if(reason == ResetReason::SATURATED) {
        if(status.num_resets >= settings::strategy::max_resets)
            return tools::log->warn("Skipped reset: num resets {} >= max resets {}", status.num_resets, settings::strategy::max_resets);
        else
            status.num_resets++; // Only increment if doing it for saturation reasons
    }
    auto t_rnd = tools::common::profile::prof[algo_type]["t_rnd"]->tic_token();
    if(not state_type) state_type = tensors.state->is_real() ? StateInitType::REAL : StateInitType::CPLX;
    if(not sector) sector = settings::strategy::target_sector;
    if(not use_eigenspinors) use_eigenspinors = settings::strategy::use_eigenspinors;
    if(not bitfield) bitfield = settings::input::bitfield;
    if(not svd_threshold and state_init == StateInit::RANDOMIZE_PREVIOUS_STATE) svd_threshold = 1e-2;
    if(not chi_lim) {
        chi_lim = cfg_chi_lim_init();
        if(not cfg_chi_lim_grow() and state_init == StateInit::RANDOMIZE_PREVIOUS_STATE)
            chi_lim = static_cast<long>(std::pow(2, std::floor(std::log2(tensors.state->find_largest_chi())))); // Nearest power of two from below
    }
    if(chi_lim.value() <= 0) throw std::runtime_error(fmt::format("Invalid chi_lim: {}", chi_lim.value()));
    svd::settings svd_settings;
    svd_settings.threshold = svd_threshold;
    tensors.randomize_state(state_init, sector.value(), chi_lim.value(), use_eigenspinors.value(), bitfield, std::nullopt, svd_settings);
    //    tensors.move_center_point_to_edge(chi_lim.value());
    clear_convergence_status();
    status.reset();
    status.iter       = 0;
    status.step       = 0;
    status.position   = tensors.state->get_position<long>();
    status.direction  = tensors.state->get_direction();
    num_perturbations = 0;
    num_chi_quenches  = 0;
    num_discards      = 0;
    num_dampings      = 0;
    stop_reason       = StopReason::NONE;
    if(cfg_chi_lim_grow()) status.chi_lim = chi_lim.value();
    if(reason == ResetReason::NEW_STATE) excited_state_number++;
    if(tensors.state->find_largest_chi() > chi_lim.value())
        //        tools::log->warn("Faulty truncation after randomize. Max found chi is {}, but chi limit is {}", tensors.state->find_largest_chi(),
        //        chi_lim.value());
        throw std::runtime_error(
            fmt::format("Faulty truncation after randomize. Max found chi is {}, but chi limit is {}", tensors.state->find_largest_chi(), chi_lim.value()));

    tensors.activate_sites(settings::precision::max_size_part_diag, 2); // Activate a pair of sites to make some measurements
    tools::log->info("Randomizing state [{}] to [{}] | Reason [{}] ... OK!", tensors.state->get_name(), enum2str(state_init), enum2str(reason));
    tools::log->info("-- Normalization            : {:.16f}", tools::finite::measure::norm(*tensors.state));
    tools::log->info("-- Spin components          : {:.6f}", fmt::join(tools::finite::measure::spin_components(*tensors.state), ", "));
    tools::log->info("-- Bond dimensions          : {}", tools::finite::measure::bond_dimensions(*tensors.state));
    tools::log->info("-- Energy per site          : {}", tools::finite::measure::energy_per_site(tensors));
    tools::log->info("-- Energy density           : {}",
                     tools::finite::measure::energy_normalized(tensors, status.energy_min_per_site, status.energy_max_per_site));
    tools::log->info("-- Energy variance          : {}", std::log10(tools::finite::measure::energy_variance(tensors)));
    tools::log->info("-- State labels             : {}", tensors.state->get_labels());
}

void class_algorithm_finite::try_projection() {
    if(not tensors.position_is_inward_edge()) return;
    if(has_projected) return;
    bool project_on_saturation = settings::strategy::project_on_saturation > 0 and
                                 status.algorithm_saturated_for > 0 and
                                 num::mod(status.algorithm_saturated_for-1, settings::strategy::project_on_saturation) == 0;

    if(settings::strategy::project_on_every_iter or project_on_saturation) {
        tools::log->info("Trying projection to {} | pos {}", settings::strategy::target_sector, tensors.get_position<long>());
        auto sector_sign  = tools::finite::mps::internal::get_sign(settings::strategy::target_sector);
        auto variance_old = tools::finite::measure::energy_variance(tensors);
        auto spincomp_old = tools::finite::measure::spin_components(*tensors.state);
        if(sector_sign != 0) {
            tensors.project_to_nearest_sector(settings::strategy::target_sector);
        } else {
            // We have a choice here.
            // If no sector sign has been given, and the spin component along the requested axis is near zero,
            // then we may inadvertently project to a sector opposite to the target state.
            // If that happened, we would get stuck in a local minima.
            // One reasonable thing to do here is to compare the variance of both projections,
            // and keep the one with lowest variance.
            // Of course, one problem is that if the spin component is already in one sector,
            // projecting to the other sector will zero the norm. So we can only make this
            // decision if the the |spin component| << 1. Maybe < 0.5 is enough?
            auto spin_component_along_requested_axis = tools::finite::measure::spin_component(*tensors.state, settings::strategy::target_sector);
            tools::log->info("Spin component along {} = {:.16f}", settings::strategy::target_sector, spin_component_along_requested_axis);
            if(std::abs(spin_component_along_requested_axis) < 0.5) {
                // Here we deem the spin component undecided enough to warrant a safe projection
                auto tensors_neg = tensors;
                auto tensors_pos = tensors;
                try {
                    tensors_neg.project_to_nearest_sector(fmt::format("-{}", settings::strategy::target_sector));
                } catch(const std::exception &ex) { tools::log->warn("Projection to -x failed: ", ex.what()); }

                try {
                    tensors_pos.project_to_nearest_sector(fmt::format("+{}", settings::strategy::target_sector));
                } catch(const std::exception &ex) { tools::log->warn("Projection to -x failed: ", ex.what()); }

                auto variance_neg = tools::finite::measure::energy_variance(tensors_neg);
                auto variance_pos = tools::finite::measure::energy_variance(tensors_pos);
                tools::log->info("Variance after projection to -{} = {:.6f}", settings::strategy::target_sector, std::log10(variance_neg));
                tools::log->info("Variance after projection to +{} = {:.6f}", settings::strategy::target_sector, std::log10(variance_pos));
                if(variance_neg < variance_pos)
                    tensors = tensors_neg;
                else
                    tensors = tensors_pos;
            } else {
                // Here the spin component is close to one sector. We just project to the nearest sector
                tensors.project_to_nearest_sector(settings::strategy::target_sector);
            }
            auto variance_new = tools::finite::measure::energy_variance(tensors);
            auto spincomp_new = tools::finite::measure::spin_components(*tensors.state);
            tools::log->info("Projection change: variance {:.6f} -> {:.6f}  | spin components {:.16f} -> {:.16f}", std::log10(variance_old),
                             std::log10(variance_new), fmt::join(spincomp_old, ", "), fmt::join(spincomp_new, ", "));
        }

        has_projected = true;
        write_to_file(StorageReason::PROJ_STATE, *tensors.state, CopyPolicy::OFF);
    }
}

void class_algorithm_finite::try_discard_small_schmidt() {
    if(settings::strategy::discard_schmidt_when_stuck == 0) return;
    if(not tensors.position_is_inward_edge()) return;
    if(num_discards >= max_discards) return;
    if(status.algorithm_has_stuck_for < 2) return;
    if(settings::strategy::discard_schmidt_when_stuck < 0)
        throw std::runtime_error(fmt::format("Expected positive discard threshold. Got: {:.16f}", settings::strategy::discard_schmidt_when_stuck));
    tools::log->info("Trying discard of smallest schmidt values: trials {}", num_discards);
    tensors.normalize_state(status.chi_lim, std::nullopt, NormPolicy::ALWAYS);
    clear_convergence_status();
    iter_discard = status.iter;
    num_discards++;
}

void class_algorithm_finite::try_bond_dimension_quench() {
    if(not settings::strategy::chi_quench_when_stuck) return;
    if(chi_quench_steps > 0) clear_convergence_status();
    if(not tensors.position_is_inward_edge()) return;
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
    if(status.algorithm_has_stuck_for == 0) {
        tools::log->info("Perturbation skipped: algorithm not stuck");
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
    // if(not state.position_is_inward_edge()) return;
    // If there are damping exponents to process, do so
    if(not damping_exponents.empty()) {
        tools::log->info("Setting damping exponent = {}", damping_exponents.back());
        tensors.damp_model_disorder(damping_exponents.back(), 0);
        damping_exponents.pop_back();
        if(damping_exponents.empty() and tensors.model->is_damped()) throw std::logic_error("Damping trial ended but the state is still damped");
    }
    if(not settings::strategy::damping_when_stuck) return;
    if(tensors.model->is_damped()) return;
    if(status.algorithm_has_stuck_for == 0) {
        tools::log->info("Damping skipped: algorithm not stuck");
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

void class_algorithm_finite::check_convergence_variance(std::optional<double> threshold, std::optional<double> saturation_sensitivity) {
    if(not tensors.position_is_inward_edge()) return;
    tools::log->trace("Checking convergence of variance mpo");
    if(not threshold) threshold = settings::precision::variance_convergence_threshold;
    if(not saturation_sensitivity) saturation_sensitivity = settings::precision::variance_saturation_sensitivity;
    var_mpo_iter.emplace_back(tools::finite::measure::energy_variance(tensors));
    auto report                       = check_saturation(var_mpo_iter, saturation_sensitivity.value());
    status.variance_mpo_converged_for = count_convergence(var_mpo_iter, threshold.value());
    if(report.has_computed) {
        status.variance_mpo_saturated_for = report.saturated_count;
        if(tools::log->level() >= spdlog::level::debug)
            tools::log->info("Energy variance convergence: last std {:7.4e} | saturated {} iters (since {})", report.Y_std.back(), report.saturated_count,
                             report.saturated_point);
        else if(tools::log->level() == spdlog::level::trace) {
            tools::log->info("Energy variance slope details:");
            tools::log->info(" -- sensitivity        = {:7.4e}", saturation_sensitivity.value());
            tools::log->info(" -- latest std         = {:7.4e}", report.Y_std.back());
            tools::log->info(" -- saturated point    = {} ", report.saturated_point);
            tools::log->info(" -- saturated count    = {} ", report.saturated_count);
            tools::log->info(" -- var saturated avg  = {:7.4e}", report.Y_avg);
            tools::log->info(" -- var history        = {:7.4e}", fmt::join(report.Y_vec, ", "));
        }
        tools::log->info("Energy variance slope details:");
        tools::log->info(" -- sensitivity        = {:7.4e}", saturation_sensitivity.value());
        tools::log->info(" -- latest std         = {:7.4e}", report.Y_std.back());
        tools::log->info(" -- saturated point    = {} ", report.saturated_point);
        tools::log->info(" -- saturated count    = {} ", report.saturated_count);
        tools::log->info(" -- var saturated avg  = {:7.4e}", report.Y_avg);
        tools::log->info(" -- var history        = {:7.4e}", fmt::join(report.Y_vec, ", "));
        tools::log->info(" -- log history        = {:7.4e}", fmt::join(report.Y_log, ", "));
        tools::log->info(" -- std history        = {:7.4e}", fmt::join(report.Y_std, ", "));
        tools::log->info(" -- ste history        = {:7.4e}", fmt::join(report.Y_ste, ", "));
    }
}

void class_algorithm_finite::check_convergence_entg_entropy(std::optional<double> saturation_sensitivity) {
    if(not tensors.position_is_inward_edge()) return;
    tools::log->trace("Checking convergence of entanglement");
    if(not saturation_sensitivity) saturation_sensitivity = settings::precision::entropy_saturation_sensitivity;
    auto                          entropies = tools::finite::measure::entanglement_entropies(*tensors.state);
    std::vector<SaturationReport> reports(entropies.size());

    for(size_t site = 0; site < entropies.size(); site++) {
        entropy_iter[site].emplace_back(entropies[site]);
        reports[site] = check_saturation(entropy_iter[site], saturation_sensitivity.value());
    }

    bool all_computed = std::all_of(reports.begin(), reports.end(), [](const SaturationReport &r) { return r.has_computed; });
    if(all_computed) {
        // Find the report which has the greatest change recently
        size_t idx_max_std = 0;
        double val_max_std = reports.front().Y_std.back();
        for(const auto &[i, r] : iter::enumerate(reports)) {
            auto &val = r.Y_std.back();
            if(std::isnan(val) or std::isinf(val)) continue;
            if(std::isnan(val_max_std) or std::isinf(val_max_std)) val_max_std = val;
            if(val >= val_max_std) {
                val_max_std = val;
                idx_max_std = i;
            }
        }
        auto &report                      = reports[idx_max_std];
        status.entanglement_saturated_for = report.saturated_count;
        std::vector<double> all_avergs;
        all_avergs.reserve(reports.size());
        for(auto &r : reports) all_avergs.push_back(r.Y_avg);
        if(tools::log->level() >= spdlog::level::debug)
            tools::log->info("Entanglement ent. convergence at site {}: last std {:7.4e} | saturated {} iters (since {})", idx_max_std, report.Y_std.back(),
                             report.saturated_count, report.saturated_point);
        else if(tools::log->level() == spdlog::level::trace) {
            tools::log->info("Entanglement slope details:");
            tools::log->info(" -- site               = {}", idx_max_std);
            tools::log->info(" -- sensitivity        = {:7.4e}", saturation_sensitivity.value());
            tools::log->info(" -- latest std         = {:7.4e}", reports[idx_max_std].Y_std.back());
            tools::log->info(" -- saturated point    = {} ", report.saturated_point);
            tools::log->info(" -- saturated count    = {} ", report.saturated_count);
            tools::log->info(" -- all averages       = {:.6f} ", fmt::join(all_avergs, ", "));
        }
        tools::log->info("Entanglement slope details:");
        tools::log->info(" -- site               = {}", idx_max_std);
        tools::log->info(" -- sensitivity        = {:7.4e}", saturation_sensitivity.value());
        tools::log->info(" -- latest std         = {:7.4e}", reports[idx_max_std].Y_std.back());
        tools::log->info(" -- saturated point    = {} ", report.saturated_point);
        tools::log->info(" -- saturated count    = {} ", report.saturated_count);
        tools::log->info(" -- all averages       = {:.6f} ", fmt::join(all_avergs, ", "));
        std::string sites;
        for(const auto &[site, s] : iter::enumerate(entropy_iter)) sites += fmt::format("{:^14}", site);
        tools::log->info(" -- sites              = {}", sites);
        for(size_t iter = 0; iter < entropy_iter.front().size(); iter++) {
            std::string vals;
            for(auto &&[site, s] : iter::enumerate(entropy_iter)) vals += fmt::format("{:12.10f}{}", s.at(iter), (site < entropy_iter.size() - 1 ? ", " : ""));
            tools::log->info(" -- ent[{:3}]           = {}", iter, vals);
        }
        //        for(auto &&[i, r] : iter::enumerate(reports))
        //            tools::log->info(" -- Ystd[{:2}]:{:3}   = {:8.2e} ",i,r.saturated_point, fmt::join(r.Y_std, ", "));
    }
    status.entanglement_converged_for = status.entanglement_saturated_for;
}

void class_algorithm_finite::check_convergence_spin_parity_sector(const std::string &target_sector, double threshold) {
    const std::array<std::string, 9> valid_sectors   = {"x", "+x", "-x", "y", "+y", "-y", "z", "+z", "-z"};
    bool                             sector_is_valid = std::find(valid_sectors.begin(), valid_sectors.end(), target_sector) != valid_sectors.end();
    if(sector_is_valid) {
        auto axis                         = tools::finite::mps::internal::get_axis(settings::strategy::target_sector);
        auto sign                         = tools::finite::mps::internal::get_sign(settings::strategy::target_sector);
        auto spin_component_along_axis    = tools::finite::measure::spin_component(*tensors.state, settings::strategy::target_sector);
        status.spin_parity_has_converged  = std::abs(std::abs(spin_component_along_axis) - 1) <= threshold ;
        if(status.spin_parity_has_converged and spin_component_along_axis * sign < 0)
            tools::log->warn("Spin parity has converged to {} = {:.16f} but requested sector was {}", axis, spin_component_along_axis, target_sector);
        tools::log->info("Spin component convergence: {} = {:.16f}", target_sector, spin_component_along_axis);
    } else
        status.spin_parity_has_converged = true; // Probably no sector was specified
}

void class_algorithm_finite::clear_convergence_status() {
    tools::log->trace("Clearing convergence status");
    for(auto &e : entropy_iter) { e.clear(); }
    var_mpo_iter.clear();
    var_mpo_step.clear();

    status.algorithm_has_finished      = false;
    status.algorithm_has_succeeded     = false;
    status.algorithm_has_to_stop       = false;
    status.algorithm_has_stuck_for     = 0;
    status.algorithm_saturated_for     = 0;
    status.algorithm_converged_for     = 0;
    status.entanglement_converged_for  = 0;
    status.entanglement_saturated_for  = 0;
    status.variance_mpo_converged_for  = 0;
    status.variance_mpo_saturated_for  = 0;
    status.chi_lim_has_reached_chi_max = false;
    status.spin_parity_has_converged   = false;
    has_projected                      = false;
    has_damped                         = false;
}

void class_algorithm_finite::setup_prefix(const StorageReason &storage_reason, StorageLevel &storage_level, const std::string &state_name,
                                          std::string &state_prefix, std::string &model_prefix, std::vector<std::string> &table_prefxs) {
    // Setup this save
    state_prefix  = fmt::format("{}/{}", algo_name, state_name); // May get modified
    model_prefix  = fmt::format("{}/{}", algo_name, "model");
    table_prefxs  = {fmt::format("{}/{}", state_prefix, "tables")}; // Common tables
    storage_level = StorageLevel::NONE;
    switch(storage_reason) {
        case StorageReason::FINISHED: {
            if(status.algorithm_has_succeeded)
                storage_level = settings::output::storage_level_good_state;
            else
                storage_level = settings::output::storage_level_fail_state;
            state_prefix += "/finished";
            table_prefxs.emplace_back(state_prefix); // Appends to its own table as well as the common ones
            break;
        }
        case StorageReason::SAVEPOINT: {
            storage_level = settings::output::storage_level_savepoint;
            if(stop_reason == StopReason::NONE) {
                if(num::mod(status.iter, settings::output::savepoint_frequency) != 0) storage_level = StorageLevel::NONE;
            }
            state_prefix += "/savepoint";
            if(settings::output::savepoint_keep_newest_only)
                state_prefix += "/iter_last";
            else
                state_prefix += fmt::format("/iter_{}", status.iter);
            table_prefxs = {state_prefix}; // Does not pollute common tables
            break;
        }
        case StorageReason::CHECKPOINT: {
            storage_level = settings::output::storage_level_checkpoint;
            if(stop_reason == StopReason::NONE) {
                if(num::mod(status.iter, settings::output::checkpoint_frequency) != 0) storage_level = StorageLevel::NONE;
            }
            state_prefix += "/checkpoint";
            if(settings::output::checkpoint_keep_newest_only)
                state_prefix += "/iter_last";
            else
                state_prefix += fmt::format("/iter_{}", status.iter);
            table_prefxs.emplace_back(state_prefix); // Appends to its own table as well as the common ones
            break;
        }
        case StorageReason::CHI_UPDATE: {
            storage_level = settings::output::storage_level_checkpoint;
            if(not cfg_chi_lim_grow()) storage_level = StorageLevel::NONE;
            // If we have updated chi we may want to write a projection too
            state_prefix += fmt::format("/checkpoint/chi_{}", status.chi_lim);
            table_prefxs = {state_prefix}; // Does not pollute common tables
            break;
        }
        case StorageReason::PROJ_STATE: {
            storage_level = settings::output::storage_level_proj_state;
            state_prefix += "/projection";
            table_prefxs = {state_prefix}; // Does not pollute common tables
            break;
        }
        case StorageReason::INIT_STATE: {
            storage_level = settings::output::storage_level_init_state;
            state_prefix += "/state_init";
            table_prefxs = {state_prefix}; // Does not pollute common tables
            break;
        }
        case StorageReason::EMIN_STATE: {
            storage_level = settings::output::storage_level_emin_state;
            break;
        }
        case StorageReason::EMAX_STATE: {
            storage_level = settings::output::storage_level_emax_state;
            break;
        }
        case StorageReason::MODEL: {
            storage_level = settings::output::storage_level_model;
            break;
        }
    }
}

void class_algorithm_finite::write_to_file(StorageReason storage_reason, std::optional<CopyPolicy> copy_policy) {
    write_to_file(storage_reason, *tensors.state, copy_policy);
}

void class_algorithm_finite::write_to_file(StorageReason storage_reason, const class_state_finite &state, std::optional<CopyPolicy> copy_policy) {
    // Setup this save
    StorageLevel             storage_level;
    std::string              state_prefix;
    std::string              model_prefix;
    std::vector<std::string> table_prefxs;
    setup_prefix(storage_reason, storage_level, state.get_name(), state_prefix, model_prefix, table_prefxs);

    switch(storage_reason) {
        case StorageReason::FINISHED: break;
        case StorageReason::SAVEPOINT:
        case StorageReason::CHECKPOINT: {
            if(stop_reason == StopReason::NONE)
                if(not state.position_is_inward_edge()) storage_level = StorageLevel::NONE;
            break;
        }
        case StorageReason::CHI_UPDATE: {
            if(not cfg_chi_lim_grow()) storage_level = StorageLevel::NONE;
            break;
        }
        case StorageReason::PROJ_STATE: {
            auto abs_spin_component = std::abs(tools::finite::measure::spin_component(state, settings::strategy::target_sector));
            if(std::abs(abs_spin_component - 1.0) > 1e-6) {
                auto state_projected = tensors.get_state_projected_to_nearest_sector(settings::strategy::target_sector, status.chi_lim);
                abs_spin_component   = std::abs(tools::finite::measure::spin_component(state_projected, settings::strategy::target_sector));
                if(std::abs(abs_spin_component - 1.0) > 1e-6)
                    throw std::runtime_error(fmt::format("Projection failed: spin {} = {:.16f}", settings::strategy::target_sector, abs_spin_component));
                return write_to_file(storage_reason, state_projected, copy_policy);
            }
            break;
        }
        case StorageReason::INIT_STATE:
        case StorageReason::EMIN_STATE:
        case StorageReason::EMAX_STATE:
        case StorageReason::MODEL: break;
    }
    if(storage_level == StorageLevel::NONE) return;
    if(state_prefix.empty()) throw std::runtime_error("State prefix is empty");
    tools::log->info("Writing to file: Reason [{}] | Level [{}] | state prefix [{}] | model prefix [{}]", enum2str(storage_reason), enum2str(storage_level),
                     state_prefix, model_prefix);

    // The file cam be kept open during writes
    h5pp_file->setKeepFileOpened();

    // Start saving tensors and metadata
    if(storage_reason == StorageReason::MODEL) {
        tools::finite::io::h5table::save_model(*h5pp_file, model_prefix + "/hamiltonian", storage_level, *tensors.model);
        tools::finite::io::h5dset::save_model(*h5pp_file, model_prefix + "/mpo", storage_level, *tensors.model);
    } else {
        tools::finite::io::h5dset::save_state(*h5pp_file, state_prefix, storage_level, state, status);
        tools::finite::io::h5dset::save_entgm(*h5pp_file, state_prefix, storage_level, state, status);
    }

    tools::common::io::h5attr::save_meta(*h5pp_file, storage_level, storage_reason, settings::model::model_type, settings::model::model_size, algo_type,
                                         state.get_name(), state_prefix, model_prefix, status);

    // The main results have now been written. Next we append data to tables
    for(const auto &table_prefix : table_prefxs) {
        if(storage_reason == StorageReason::MODEL) break;
        tools::finite::io::h5table::save_sim_status(*h5pp_file, table_prefix + "/status", storage_level, status);
        tools::finite::io::h5table::save_profiling(*h5pp_file, table_prefix + "/profiling", storage_level, status);
        tools::finite::io::h5table::save_mem_usage(*h5pp_file, table_prefix + "/mem_usage", storage_level, status);
        tools::finite::io::h5table::save_measurements(*h5pp_file, table_prefix + "/measurements", storage_level, tensors, status, algo_type);
    }
    h5pp_file->setKeepFileClosed();

    // Copy from temporary location to destination depending on given policy
    auto t_hdf = tools::common::profile::prof[algo_type]["t_hdf"]->tic_token();
    copy_from_tmp(storage_reason, copy_policy);
}

template<typename T>
void class_algorithm_finite::write_to_file(StorageReason storage_reason, const T &data, const std::string &name, std::optional<CopyPolicy> copy_policy) {
    // Setup this save
    StorageLevel             storage_level;
    std::string              state_prefix;
    std::string              model_prefix;
    std::vector<std::string> table_prefxs;
    setup_prefix(storage_reason, storage_level, tensors.state->get_name(), state_prefix, model_prefix, table_prefxs);

    std::string data_path = fmt::format("{}/{}", state_prefix, name);
    if(storage_reason == StorageReason::MODEL) data_path = fmt::format("{}/{}", model_prefix, name);

    tools::finite::io::h5dset::save_data(*h5pp_file, data, data_path, status);

    // Copy from temporary location to destination depending on given policy
    auto t_hdf = tools::common::profile::prof[algo_type]["t_hdf"]->tic_token();
    copy_from_tmp(storage_reason, copy_policy);
}

template void class_algorithm_finite::write_to_file(StorageReason storage_reason, const Eigen::Tensor<Scalar, 2> &data, const std::string &name,
                                                    std::optional<CopyPolicy> copy_policy);

void class_algorithm_finite::print_status_update() {
    if(num::mod(status.step, cfg_print_freq()) != 0) return;
    if(cfg_print_freq() == 0) return;
    if(tensors.state->position_is_outward_edge()) return;

    std::string report;
    //    report += fmt::format("{:<} ", algo_name);
    report += fmt::format("{:<} ", tensors.state->get_name());
    report += fmt::format("iter:{:<4} ", status.iter);
    report += fmt::format("step:{:<5} ", status.step);
    report += fmt::format("L:{} ", tensors.get_length());
    if(tensors.active_sites.empty())
        report += fmt::format("l:[{:>2} {:<2}] ", status.position, status.position);
    else if(tensors.state->get_direction() > 0)
        report += fmt::format("l:[{:>2} {:<2}] ", tensors.active_sites.front(), tensors.active_sites.back());
    else if(tensors.state->get_direction() < 0)
        report += fmt::format("l:[{:>2} {:<2}] ", tensors.active_sites.back(), tensors.active_sites.front());

    double energy = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : tools::finite::measure::energy_per_site(tensors);
    report += fmt::format("E/L:{:<20.16f} ", energy);

    if(algo_type == AlgorithmType::xDMRG) { report += fmt::format("ε:{:<6.4f} ", status.energy_dens); }
    report += fmt::format("Sₑ(l):{:<10.8f} ", tools::finite::measure::entanglement_entropy_current(*tensors.state));

    double variance = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : std::log10(tools::finite::measure::energy_variance(tensors));
    report += fmt::format("log₁₀σ²E:{:<10.6f} [{:<10.6f}] ", variance, std::log10(status.energy_variance_lowest));
    report += fmt::format("χ:{:<3}|{:<3}|", cfg_chi_lim_max(), status.chi_lim);
    size_t comma_width       = settings::strategy::multisite_mps_size_max <= 2 ? 0 : 2; // ", "
    size_t bond_single_width = static_cast<size_t>(std::log10(cfg_chi_lim_max())) + 1;
    size_t bond_string_width =
        2 + (bond_single_width + comma_width) * (settings::strategy::multisite_mps_size_max == 1 ? 1 : settings::strategy::multisite_mps_size_max - 1);
    std::string bond_string = fmt::format("{}", tools::finite::measure::bond_dimensions_merged(*tensors.state));
    report += fmt::format("{0:<{1}} ", bond_string, bond_string_width);

    if(last_optmode and last_optspace)
        report += fmt::format("opt:[{}|{}] ", enum2str(last_optmode.value()).substr(0, 3), enum2str(last_optspace.value()).substr(0, 3));
    report += fmt::format("log₁₀trnc:{:<8.4f} ", std::log10(tensors.state->get_truncation_error()));
    report += fmt::format("stk:{:<1} ", status.algorithm_has_stuck_for);
    report += fmt::format("sat:[σ² {:<1} Sₑ {:<1}] ", status.variance_mpo_saturated_for, status.entanglement_saturated_for);
    report += fmt::format("time:{:<} ", fmt::format("{:>6.2f}s", tools::common::profile::t_tot->get_measured_time()));
    report += fmt::format("mem[rss {:<.1f}|peak {:<.1f}|vm {:<.1f}]MB ", tools::common::profile::mem_rss_in_mb(), tools::common::profile::mem_hwm_in_mb(),
                          tools::common::profile::mem_vm_in_mb());
    tools::log->info(report);
}

void class_algorithm_finite::print_status_full() {
    tensors.rebuild_edges();
    tensors.do_all_measurements();
    tools::log->info("{:=^60}", "");
    tools::log->info("= {: ^56} =", " Full status: [" + algo_name + "][" + tensors.state->get_name() + "]");
    tools::log->info("{:=^60}", "");
    tools::log->info("Stop reason                        = {}", enum2str(stop_reason));
    tools::log->info("Sites                              = {}", tensors.get_length());
    tools::log->info("Position                           = {}", status.position);
    tools::log->info("Direction                          = {}", status.direction);
    tools::log->info("Iterations (full chain sweeps)     = {}", status.iter);
    tools::log->info("Steps (moves along the chain)      = {}", status.step);
    tools::log->info("Total time                         = {:<.1f} s = {:<.2f} min", tools::common::profile::t_tot->get_measured_time(),
                     tools::common::profile::t_tot->get_measured_time() / 60);
    double energy_per_site = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : tools::finite::measure::energy_per_site(tensors);
    tools::log->info("Energy per site E/L                = {:<.16f}", energy_per_site);
    if(algo_type == AlgorithmType::xDMRG)
        tools::log->info("Energy density (rescaled 0 to 1) ε = {:<6.4f}",
                         tools::finite::measure::energy_normalized(tensors, status.energy_min_per_site, status.energy_max_per_site));
    double variance = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : std::log10(tools::finite::measure::energy_variance(tensors));
    tools::log->info("Energy variance log₁₀ σ²(E)        = {:<.16f}", variance);
    tools::log->info("Bond dimension maximum χmax        = {}", cfg_chi_lim_max());
    tools::log->info("Bond dimensions χ                  = {}", tools::finite::measure::bond_dimensions(*tensors.state));
    tools::log->info("Bond dimension  χ (mid)            = {}", tools::finite::measure::bond_dimension_midchain(*tensors.state));
    tools::log->info("Entanglement entropies Sₑ          = {:.5f}", fmt::join(tools::finite::measure::entanglement_entropies(*tensors.state), ", "));
    tools::log->info("Entanglement entropy   Sₑ (mid)    = {:.5f}", tools::finite::measure::entanglement_entropy_midchain(*tensors.state), ", ");
    if(algo_type == AlgorithmType::fLBIT) {
        tools::log->info("Number entropies Sₙ                = {:.5f}", fmt::join(tools::finite::measure::number_entropies(*tensors.state), ", "));
        tools::log->info("Number entropy   Sₙ (mid)          = {:.5f}", tools::finite::measure::number_entropy_midchain(*tensors.state), ", ");
    }
    tools::log->info("Spin components                    = {:.5f}", fmt::join(tools::finite::measure::spin_components(*tensors.state), ", "));
    tools::log->info("Truncation Errors                  = {:.3e}", fmt::join(tensors.state->get_truncation_errors(), ", "));
    tools::log->info("Algorithm has succeeded            = {:<}", status.algorithm_has_succeeded);
    tools::log->info("Algorithm has saturated for        = {:<}", status.algorithm_saturated_for);
    tools::log->info("Algorithm has got stuck for        = {:<}", status.algorithm_has_stuck_for);
    tools::log->info("Algorithm has converged for        = {:<}", status.algorithm_converged_for);

    tools::log->info("σ²                                 = Converged : {:<4}  Saturated: {:<4}", status.variance_mpo_converged_for,
                     status.variance_mpo_saturated_for);
    tools::log->info("Sₑ                                 = Converged : {:<4}  Saturated: {:<4}", status.entanglement_converged_for,
                     status.entanglement_saturated_for);
    tools::log->info("{:=^60}", "");
}
