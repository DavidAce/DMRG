#include "AlgorithmFinite.h"
#include <config/settings.h>
#include <debug/exceptions.h>
#include <general/iter.h>
#include <h5pp/h5pp.h>
#include <math/num.h>
#include <tensors/edges/EdgesFinite.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/state/StateFinite.h>
#include <tid/tid.h>
#include <tools/common/h5.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/env.h>
#include <tools/finite/h5.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <tools/finite/print.h>

AlgorithmFinite::AlgorithmFinite(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type) : AlgorithmBase(std::move(h5ppFile_), algo_type) {
    tools::log->trace("Constructing class_algorithm_finite");
    tensors.initialize(settings::model::model_type, settings::model::model_size, 0);
    tensors.state->set_algorithm(algo_type);
    entropy_iter.resize(tensors.get_length() + 1);
}

// We need to make a destructor manually for the enclosing class "ModelFinite"
// that encloses "class_model_base". Otherwise unique_ptr will forcibly inline its
// own default deleter.
// This allows us to forward declare the abstract base class "class_model_base"
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
// class_algorithm_finite::~class_algorithm_finite() = default;

void AlgorithmFinite::run()
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
    tools::log->info("Starting {}", status.algo_type_sv());
    auto t_tot  = tid::get_unscoped("t_tot").tic_token();
    auto t_algo = tid::tic_scope(status.algo_type_sv());
    // We may want to resume this simulation.
    if(h5pp_file->linkExists("common/storage_level") and
       (settings::output::file_collision_policy == FileCollisionPolicy::RESUME or settings::output::file_collision_policy == FileCollisionPolicy::REVIVE)) {
        try {
            resume();
        } catch(const except::state_error &ex) {
            throw except::resume_error(fmt::format("Could not resume state from file [{}]: {}", h5pp_file->getFilePath(), ex.what()));
        } catch(const except::load_error &ex) {
            throw except::resume_error(fmt::format("Could not resume state from file [{}]: {}", h5pp_file->getFilePath(), ex.what()));
        } catch(const std::exception &ex) {
            throw std::runtime_error(fmt::format("Could not resume state from file [{}]: {}", h5pp_file->getFilePath(), ex.what()));
        }
    } else {
        run_default_task_list();
    }
}

void AlgorithmFinite::run_postprocessing() {
    tools::log->info("Running default postprocessing for {}", status.algo_type_sv());
    auto tic = tid::tic_scope("post");
    write_to_file(StorageReason::CHECKPOINT, CopyPolicy::OFF);
    write_to_file(StorageReason::PROJ_STATE, CopyPolicy::OFF);
    write_to_file(StorageReason::FINISHED, CopyPolicy::FORCE);
    print_status_full();
    tools::log->info("Finished default postprocessing for {}", status.algo_type_sv());
    tools::common::profile::print_profiling(status);
}

void AlgorithmFinite::move_center_point(std::optional<long> num_moves) {
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
            tools::log->trace("Moving center position step {}", status.step);
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

void AlgorithmFinite::reduce_mpo_energy() {
    if(not tensors.position_is_inward_edge()) return;
    if(not settings::precision::use_reduced_mpo_energy) return;
    // Reduce mpo energy to avoid catastrophic cancellation
    // Note that this operation makes the Hamiltonian nearly singular,
    // which is tough for Lanczos/Arnoldi iterations to handle
    tensors.reduce_mpo_energy(std::nullopt);
    // The reduction clears our squared mpo's. So we have to rebuild.
    rebuild_mpo_squared();
}

void AlgorithmFinite::rebuild_mpo_squared() {
    if(not tensors.position_is_inward_edge()) return;
    bool compress = settings::precision::use_compressed_mpo_squared_all;
    tensors.rebuild_mpo_squared(compress);
}

void AlgorithmFinite::update_variance_max_digits(std::optional<double> energy) {
    if(not tensors.position_is_inward_edge()) return;
    if(tensors.active_sites.empty()) return;
    if(not energy) energy = tools::finite::measure::energy(tensors);
    double energy_abs                 = std::abs(energy.value());
    double energy_pow                 = energy_abs * energy_abs;
    double digits10                   = std::numeric_limits<double>::digits10;
    double energy_top                 = settings::precision::use_reduced_mpo_energy ? energy_abs : energy_pow;
    double energy_exp                 = std::ceil(std::max(0.0, std::log10(energy_top))) + 1;
    double max_digits                 = std::floor(std::max(0.0, digits10 - energy_exp));
    status.energy_variance_max_digits = static_cast<size_t>(max_digits);
    status.energy_variance_prec_limit = std::pow(10.0, -max_digits);
    tools::log->debug("Estimated limit on energy variance precision: {:.3e}", status.energy_variance_prec_limit);
}

void AlgorithmFinite::update_bond_dimension_limit() {
    if(not tensors.position_is_inward_edge()) return;
    status.chi_lim_max                 = settings::chi_lim_max(status.algo_type);
    status.chi_lim_has_reached_chi_max = status.chi_lim >= status.chi_lim_max;
    if(settings::chi_lim_grow(status.algo_type) == ChiGrow::OFF) {
        status.chi_lim = status.chi_lim_max;
        return;
    }
    if(status.chi_lim_has_reached_chi_max) return;
    auto tic = tid::tic_scope("chi_grow");


    if constexpr (settings::debug){
        if(tools::log->level() == spdlog::level::trace) {
            double truncation_threshold = 2 * settings::precision::svd_threshold;
            size_t trunc_bond_count     = tensors.state->num_sites_truncated(truncation_threshold);
            size_t bond_at_lim_count    = tensors.state->num_bonds_reached_chi(status.chi_lim);
            tools::log->trace("Truncation threshold  : {:<.8e}", truncation_threshold);
            tools::log->trace("Truncation errors     : {}", tensors.state->get_truncation_errors());
            tools::log->trace("Bond dimensions       : {}", tools::finite::measure::bond_dimensions(*tensors.state));
            tools::log->trace("Truncated bond count  : {} ", trunc_bond_count);
            tools::log->trace("Bonds at limit  count : {} ", bond_at_lim_count);
            tools::log->trace("Entanglement entropies: {} ", tools::finite::measure::entanglement_entropies(*tensors.state));
        }
    }


    // If we got here we want to increase the bond dimension limit progressively during the simulation
    // Only increment the bond dimension if the following are all true
    //      * No experiments are on-going like perturbation
    //      * the state precision is limited by bond dimension
    // In addition, if chi_lim_grow == ChiGrow::ON_SATURATION we add the condition
    //      * the algorithm has got stuck

    // When schmidt values are highly truncated at every step the entanglement fluctuates a lot, so we should check both
    // variance and entanglement for saturation. Note that status.algorithm_saturaded_for uses an "and" condition.
    bool is_saturated       = status.entanglement_saturated_for > 0 or status.variance_mpo_saturated_for > 0;
    bool is_exp_ongoing     = tensors.model->is_perturbed();
    bool is_bond_limited    = tensors.state->is_bond_limited(status.chi_lim, 2 * settings::precision::svd_threshold);
    bool grow_on_saturation = settings::chi_lim_grow(status.algo_type) == ChiGrow::ON_SATURATION;
    if(is_exp_ongoing) {
        tools::log->info("State is undergoing perturbation -- cannot increase bond dimension yet");
        return;
    }
    if(grow_on_saturation and not is_saturated) {
        tools::log->info("Algorithm is not saturated yet. Kept current limit {}", status.chi_lim);
        return;
    }

    if(not is_bond_limited) {
        tools::log->info("State is not limited by its bond dimension. Kept current limit {}", status.chi_lim);
        return;
    }

    // If we got to this point we will update the bond dimension by a factor
    auto factor = settings::chi_lim_grow_factor(status.algo_type);
    if(factor <= 1.0) throw std::runtime_error(fmt::format("Error: chi_lim_grow_factor == {:.3f} | must be larger than one", factor));


    // Write current results before updating bond dimension
    write_to_file(StorageReason::CHI_UPDATE);
    if(settings::strategy::randomize_on_chi_update and status.chi_lim >= 32)
        randomize_state(ResetReason::CHI_UPDATE, StateInit::RANDOMIZE_PREVIOUS_STATE, std::nullopt, std::nullopt, status.chi_lim);

    double chi_prod = std::ceil(factor * static_cast<double>(status.chi_lim));
    long chi_new = std::min(static_cast<long>(chi_prod), status.chi_lim_max);
    tools::log->info("Updating bond dimension limit {} -> {}", status.chi_lim, chi_new);
    status.chi_lim                     = chi_new;
    status.chi_lim_has_reached_chi_max = status.chi_lim == status.chi_lim_max;

    // Last sanity check before leaving here
    if(status.chi_lim > status.chi_lim_max)
        throw std::runtime_error(fmt::format("chi_lim is larger than chi_lim_max! {} > {}", status.chi_lim, status.chi_lim_max));
}

void AlgorithmFinite::randomize_model() {
    tools::log->info("Randomizing model");
    auto tic = tid::tic_scope("init.rnd_model");
    tensors.randomize_model();
    clear_convergence_status();
}

void AlgorithmFinite::randomize_state(ResetReason reason, StateInit state_init, std::optional<StateInitType> state_type, std::optional<std::string> sector,
                                      std::optional<long> chi_lim, std::optional<bool> use_eigenspinors, std::optional<long> bitfield,
                                      std::optional<double> svd_threshold) {
    tools::log->info("Randomizing state [{}] to [{}] | Reason [{}] ...", tensors.state->get_name(), enum2sv(state_init), enum2sv(reason));
    auto t_rnd = tid::tic_scope("init.rnd_state");
    if(reason == ResetReason::SATURATED) {
        if(status.num_resets >= settings::strategy::max_resets)
            return tools::log->warn("Skipped reset: num resets {} >= max resets {}", status.num_resets, settings::strategy::max_resets);
        else
            status.num_resets++; // Only increment if doing it for saturation reasons
    }
    if(not state_type) state_type = tensors.state->is_real() ? StateInitType::REAL : StateInitType::CPLX;
    if(not sector) sector = settings::strategy::target_sector;
    if(not use_eigenspinors) use_eigenspinors = settings::strategy::use_eigenspinors;
    if(not bitfield) bitfield = settings::input::bitfield;
    if(not svd_threshold and state_init == StateInit::RANDOMIZE_PREVIOUS_STATE) svd_threshold = 1e-2;
    if(not chi_lim) {
        chi_lim = settings::chi_lim_init(status.algo_type);
        if(settings::chi_lim_grow(status.algo_type) == ChiGrow::OFF and state_init == StateInit::RANDOMIZE_PREVIOUS_STATE)
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
    status.algo_stop  = AlgorithmStop::NONE;
    if(settings::chi_lim_grow(status.algo_type) != ChiGrow::OFF) status.chi_lim = chi_lim.value();
    if(reason == ResetReason::NEW_STATE) excited_state_number++;
    if(tensors.state->find_largest_chi() > chi_lim.value())
        //        tools::log->warn("Faulty truncation after randomize. Max found chi is {}, but chi limit is {}", tensors.state->find_largest_chi(),
        //        chi_lim.value());
        throw std::runtime_error(
            fmt::format("Faulty truncation after randomize. Max found chi is {}, but chi limit is {}", tensors.state->find_largest_chi(), chi_lim.value()));

    tensors.activate_sites(settings::precision::max_size_part_diag, 2); // Activate a pair of sites to make some measurements
    tools::log->info("Randomizing state [{}] to [{}] | Reason [{}] ... OK!", tensors.state->get_name(), enum2sv(state_init), enum2sv(reason));
    tools::log->info("-- Normalization            : {:.16f}", tools::finite::measure::norm(*tensors.state));
    tools::log->info("-- Spin components          : {:.6f}", fmt::join(tools::finite::measure::spin_components(*tensors.state), ", "));
    tools::log->info("-- Bond dimensions          : {}", tools::finite::measure::bond_dimensions(*tensors.state));
    tools::log->info("-- Energy per site          : {}", tools::finite::measure::energy_per_site(tensors));
    tools::log->info("-- Energy density           : {}",
                     tools::finite::measure::energy_normalized(tensors, status.energy_min_per_site, status.energy_max_per_site));
    tools::log->info("-- Energy variance          : {}", std::log10(tools::finite::measure::energy_variance(tensors)));
    tools::log->info("-- State labels             : {}", tensors.state->get_labels());
}

void AlgorithmFinite::try_projection(std::optional<std::string> target_sector) {
    if(not tensors.position_is_inward_edge()) return;
    if(not target_sector and has_projected) return;
    if(status.variance_mpo_converged_for > 0 and status.spin_parity_has_converged) return; // No need

    bool project_on_saturation = settings::strategy::project_on_saturation > 0 and status.algorithm_saturated_for > 0 and
                                 num::mod(status.algorithm_saturated_for - 1, settings::strategy::project_on_saturation) == 0;
    bool project_on_every_iter   = settings::strategy::project_on_every_iter > 0 and num::mod(status.iter, settings::strategy::project_on_every_iter) == 0;
    bool project_to_given_sector = target_sector.has_value();

    if(project_on_every_iter or project_on_saturation or project_to_given_sector) {
        if(not target_sector) target_sector = settings::strategy::target_sector;
        tools::log->info("Trying projection to {}", target_sector.value());
        auto sector_sign  = tools::finite::mps::init::get_sign(target_sector.value());
        auto variance_old = tools::finite::measure::energy_variance(tensors);
        auto spincomp_old = tools::finite::measure::spin_components(*tensors.state);

        if(sector_sign != 0) {
            tensors.project_to_nearest_sector(target_sector.value());
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
            auto spin_component_along_requested_axis = tools::finite::measure::spin_component(*tensors.state, target_sector.value());
            tools::log->debug("Spin component along {} = {:.16f}", target_sector.value(), spin_component_along_requested_axis);
            if(std::abs(spin_component_along_requested_axis) < 0.5) {
                // Here we deem the spin component undecided enough to warrant a safe projection
                auto tensors_neg = tensors;
                auto tensors_pos = tensors;
                try {
                    tensors_neg.project_to_nearest_sector(fmt::format("-{}", target_sector.value()), std::nullopt);
                } catch(const std::exception &ex) { tools::log->warn("Projection to -x failed: ", ex.what()); }

                try {
                    tensors_pos.project_to_nearest_sector(fmt::format("+{}", target_sector.value()), std::nullopt);
                } catch(const std::exception &ex) { tools::log->warn("Projection to -x failed: ", ex.what()); }

                auto variance_neg = tools::finite::measure::energy_variance(tensors_neg);
                auto variance_pos = tools::finite::measure::energy_variance(tensors_pos);
                tools::log->debug("Variance after projection to -{} = {:.6f}", target_sector.value(), std::log10(variance_neg));
                tools::log->debug("Variance after projection to +{} = {:.6f}", target_sector.value(), std::log10(variance_pos));
                if(variance_neg < variance_pos)
                    tensors = tensors_neg;
                else
                    tensors = tensors_pos;
            } else {
                // Here the spin component is close to one sector. We just project to the nearest sector
                tensors.project_to_nearest_sector(target_sector.value(), std::nullopt);
            }
            auto variance_new = tools::finite::measure::energy_variance(tensors);
            auto spincomp_new = tools::finite::measure::spin_components(*tensors.state);
            tools::log->debug("Projection change: variance {:.6f} -> {:.6f}  | spin components {:.16f} -> {:.16f}", std::log10(variance_old),
                              std::log10(variance_new), fmt::join(spincomp_old, ", "), fmt::join(spincomp_new, ", "));
        }
        if(target_sector.value() == settings::strategy::target_sector) has_projected = true;
        write_to_file(StorageReason::PROJ_STATE, CopyPolicy::OFF);
    }
}

void AlgorithmFinite::try_full_expansion() {
    if(not tensors.position_is_inward_edge()) return;
    if(has_expanded) return;
    bool expand_on_saturation = settings::strategy::expand_on_saturation > 0 and status.algorithm_saturated_for > 0 and
                                num::mod(status.algorithm_saturated_for - 1, settings::strategy::expand_on_saturation) == 0;
    if(expand_on_saturation) {
        auto variance_old = tools::finite::measure::energy_variance(tensors);
        tools::log->info("Trying expansion to H|psi> | pos {}", tensors.get_position());
        tensors.apply_hamiltonian_on_state(status.chi_lim);
        auto variance_new = tools::finite::measure::energy_variance(tensors);
        has_expanded      = true;
        tools::log->info("Expansion change: variance {:.6f} -> {:.6f}", std::log10(variance_old), std::log10(variance_new));
    }
}

void AlgorithmFinite::try_bond_dimension_quench() {
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
    if(tools::finite::measure::energy_variance(tensors) < 10 * std::max(status.energy_variance_prec_limit, settings::precision::variance_convergence_threshold))
        return;
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

void AlgorithmFinite::try_hamiltonian_perturbation() {
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

void AlgorithmFinite::check_convergence_variance(std::optional<double> threshold, std::optional<double> saturation_sensitivity) {
    if(not tensors.position_is_inward_edge()) return;
    tools::log->trace("Checking convergence of variance mpo");
    if(not threshold) threshold = std::max(status.energy_variance_prec_limit, settings::precision::variance_convergence_threshold);
    if(not saturation_sensitivity) saturation_sensitivity = settings::precision::variance_saturation_sensitivity;
    var_mpo_iter.emplace_back(tools::finite::measure::energy_variance(tensors));
    auto report = check_saturation(var_mpo_iter, saturation_sensitivity.value());
    if(report.has_computed) {
        status.variance_mpo_converged_for = count_convergence(var_mpo_iter, threshold.value(), report.saturated_point);
        status.variance_mpo_saturated_for = report.saturated_count;
        if(tools::log->level() >= spdlog::level::debug)
            tools::log->debug("Energy variance convergence: saturated {} iters (since iter {})", report.saturated_count, report.saturated_point);
        if(tools::log->level() == spdlog::level::debug) {
            tools::log->debug("Energy variance slope details:");
            tools::log->debug(" -- sensitivity        = {:7.4e}", saturation_sensitivity.value());
            tools::log->debug(" -- threshold          = {:7.4e}", threshold.value());
            tools::log->debug(" -- saturated point    = {} ", report.saturated_point);
            tools::log->debug(" -- saturated count    = {} ", report.saturated_count);
            tools::log->debug(" -- converged count    = {} ", status.variance_mpo_converged_for);
            tools::log->debug(" -- var history        = {:7.4e}", fmt::join(report.Y_vec, ", "));
            tools::log->debug(" -- avg history        = {:7.4e}", fmt::join(report.Y_avg, ", "));
            tools::log->debug(" -- std history        = {:7.4e}", fmt::join(report.Y_std, ", "));
            tools::log->debug(" -- stn history        = {:7.4e}", fmt::join(report.Y_stn, ", "));
            tools::log->debug(" -- slp history        = {:7.4e}", fmt::join(report.Y_slp, ", "));
        }
    }
}

void AlgorithmFinite::check_convergence_entg_entropy(std::optional<double> saturation_sensitivity) {
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
        size_t last_saturated_site = 0;
        size_t last_saturated_idx  = 0;
        for(const auto &[site, r] : iter::enumerate(reports)) {
            if(r.saturated_point >= last_saturated_idx) {
                last_saturated_site = site;
                last_saturated_idx  = r.saturated_point;
            }
        }
        auto &report                      = reports[last_saturated_site];
        status.entanglement_saturated_for = report.saturated_count;
        if(tools::log->level() >= spdlog::level::debug)
            tools::log->debug("Entanglement ent. convergence at site {}: saturated {} iters (since {})", last_saturated_site, report.saturated_count,
                              report.saturated_point);
        if(tools::log->level() == spdlog::level::debug) {
            tools::log->debug("Entanglement slope details:");
            tools::log->debug(" -- site               = {}", last_saturated_site);
            tools::log->debug(" -- sensitivity        = {:7.4e}", saturation_sensitivity.value());
            tools::log->debug(" -- saturated point    = {} ", report.saturated_point);
            tools::log->debug(" -- saturated count    = {} ", report.saturated_count);
            tools::log->debug(" -- ent history        = {:7.4e}", fmt::join(report.Y_vec, ", "));
            tools::log->debug(" -- avg history        = {:7.4e}", fmt::join(report.Y_avg, ", "));
            tools::log->debug(" -- std history        = {:7.4e}", fmt::join(report.Y_std, ", "));
            tools::log->debug(" -- stn history        = {:7.4e}", fmt::join(report.Y_stn, ", "));
        }
    }
    status.entanglement_converged_for = status.entanglement_saturated_for;
}

void AlgorithmFinite::check_convergence_spin_parity_sector(std::string_view target_sector, double threshold) {
    constexpr std::array<std::string_view, 9> valid_sectors   = {"x", "+x", "-x", "y", "+y", "-y", "z", "+z", "-z"};
    bool                                      sector_is_valid = std::find(valid_sectors.begin(), valid_sectors.end(), target_sector) != valid_sectors.end();
    if(sector_is_valid) {
        auto axis                        = tools::finite::mps::init::get_axis(settings::strategy::target_sector);
        auto sign                        = tools::finite::mps::init::get_sign(settings::strategy::target_sector);
        auto spin_component_along_axis   = tools::finite::measure::spin_component(*tensors.state, settings::strategy::target_sector);
        status.spin_parity_has_converged = std::abs(std::abs(spin_component_along_axis) - 1) <= threshold;
        if(status.spin_parity_has_converged and spin_component_along_axis * sign < 0)
            tools::log->warn("Spin parity has converged to {} = {:.16f} but requested sector was {}", axis, spin_component_along_axis, target_sector);
        tools::log->debug("Spin component convergence: {} = {:.16f}", target_sector, spin_component_along_axis);
    } else
        status.spin_parity_has_converged = true; // Probably no sector was specified
}

void AlgorithmFinite::clear_convergence_status() {
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
}

void AlgorithmFinite::write_to_file(StorageReason storage_reason, std::optional<CopyPolicy> copy_policy) {
    tools::finite::h5::save::simulation(*h5pp_file, tensors, status, storage_reason, copy_policy);
}

void AlgorithmFinite::write_to_file(StorageReason storage_reason, const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges,
                                    std::optional<CopyPolicy> copy_policy) {
    tools::finite::h5::save::simulation(*h5pp_file, state, model, edges, status, storage_reason, copy_policy);
}

template<typename T>
void AlgorithmFinite::write_to_file(StorageReason storage_reason, const T &data, std::string_view name, std::optional<CopyPolicy> copy_policy) {
    tools::finite::h5::save::data(*h5pp_file, data, name, tensors.state->get_name(), status, storage_reason, copy_policy);
}

template void AlgorithmFinite::write_to_file(StorageReason storage_reason, const Eigen::Tensor<std::complex<double>, 2> &data, std::string_view name,
                                             std::optional<CopyPolicy> copy_policy);

void AlgorithmFinite::print_status_update() {
    if(num::mod(status.step, settings::print_freq(status.algo_type)) != 0) return;
    if(settings::print_freq(status.algo_type) == 0) return;
    if(tensors.state->position_is_outward_edge()) return;

    std::string report;
    //    report += fmt::format("{:<} ", status.algo_type_sv());
    report += fmt::format(FMT_STRING("{:<} "), tensors.state->get_name());
    report += fmt::format(FMT_STRING("iter:{:<4} "), status.iter);
    report += fmt::format(FMT_STRING("step:{:<5} "), status.step);
    report += fmt::format(FMT_STRING("L:{} "), tensors.get_length());
    std::string site_str;
    if(tensors.active_sites.empty()) site_str = fmt::format(FMT_STRING("{:^5}"), "?");
    if(tensors.active_sites.size() == 1) site_str = fmt::format(FMT_STRING("{:^5}"), tensors.active_sites.front());
    if(tensors.active_sites.size() >= 2) site_str = fmt::format(FMT_STRING("{:>2} {:>2}"), tensors.active_sites.front(), tensors.active_sites.back());
    if(tensors.state->get_direction() > 0) {
        report += fmt::format("l:|{}⟩ ", site_str);
    } else if(tensors.state->get_direction() < 0) {
        report += fmt::format("l:⟨{}| ", site_str);
    }

    double energy = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : tools::finite::measure::energy_per_site(tensors);
    report += fmt::format(FMT_STRING("E/L:{:<20.16f} "), energy);

    if(status.algo_type == AlgorithmType::xDMRG) { report += fmt::format(FMT_STRING("ε:{:<6.4f} "), status.energy_dens); }
    report += fmt::format(FMT_STRING("Sₑ({:>2}):{:<10.8f} "), tensors.state->get_position<long>(), tools::finite::measure::entanglement_entropy_current(*tensors.state));

    double variance = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : std::log10(tools::finite::measure::energy_variance(tensors));
    report += fmt::format(FMT_STRING("log₁₀σ²E:{:<10.6f} [{:<10.6f}] "), variance, std::log10(status.energy_variance_lowest));
    report += fmt::format(FMT_STRING("χ:{:<3}|{:<3}|"), settings::chi_lim_max(status.algo_type), status.chi_lim);
    size_t comma_width       = settings::strategy::multisite_mps_size_max <= 2 ? 0 : 2; // ", "
    size_t bracket_width     = 2;                                                       // The {} edges
    size_t bond_single_width = static_cast<size_t>(std::log10(settings::chi_lim_max(status.algo_type))) + 1;
    size_t bond_num_elements = settings::strategy::multisite_mps_size_max == 1 ? 1 : settings::strategy::multisite_mps_size_max - 1;
    size_t bond_string_width = bracket_width + (bond_single_width + comma_width) // The width of an element like " 54,"
                                                   * bond_num_elements;          // Number of bonds
    std::vector<long> bonds_merged = tools::finite::measure::bond_dimensions_merged(*tensors.state);
    if(bonds_merged.empty())
        report += fmt::format(FMT_STRING("{0:<{1}} "), " ", bond_string_width);
    else {
        std::string bonds_string = fmt::format("{}", bonds_merged);
        report += fmt::format(FMT_STRING("{0:<{1}} "), bonds_string, bond_string_width);
    }

    if(last_optmode and last_optspace)
        report += fmt::format(FMT_STRING("opt:[{}|{}] "), enum2sv(last_optmode.value()).substr(0, 3), enum2sv(last_optspace.value()).substr(0, 3));
    report += fmt::format(FMT_STRING("log₁₀trnc:{:<8.4f} "), std::log10(tensors.state->get_truncation_error()));
    report += fmt::format(FMT_STRING("stk:{:<1} "), status.algorithm_has_stuck_for);
    report += fmt::format(FMT_STRING("sat:[σ² {:<1} Sₑ {:<1}] "), status.variance_mpo_saturated_for, status.entanglement_saturated_for);
    report += fmt::format(FMT_STRING("time:{:<9} "), fmt::format("{:>7.1f}s", tid::get_unscoped("t_tot").get_time()));
    report += fmt::format(FMT_STRING("mem[rss {:<.1f}|peak {:<.1f}|vm {:<.1f}]MB "), tools::common::profile::mem_rss_in_mb(), tools::common::profile::mem_hwm_in_mb(),
                          tools::common::profile::mem_vm_in_mb());
    tools::log->info(report);
}

void AlgorithmFinite::print_status_full() {
    tensors.rebuild_edges();
    tensors.do_all_measurements();
    tools::log->info("{:=^60}", "");
    tools::log->info("= {: ^56} =", fmt::format("Full status [{}][{}]", status.algo_type_sv(), tensors.state->get_name()));
    tools::log->info("{:=^60}", "");
    tools::log->info("Stop reason                        = {}", status.algo_stop_sv());
    tools::log->info("Sites                              = {}", tensors.get_length());
    tools::log->info("Position                           = {}", status.position);
    tools::log->info("Direction                          = {}", status.direction);
    tools::log->info("Iterations (full chain sweeps)     = {}", status.iter);
    tools::log->info("Steps (moves along the chain)      = {}", status.step);
    tools::log->info("Total time                         = {:<.1f} s = {:<.2f} min", tid::get_unscoped("t_tot").get_time(),
                     tid::get_unscoped("t_tot").get_time() / 60);
    double energy_per_site = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : tools::finite::measure::energy_per_site(tensors);
    tools::log->info("Energy per site E/L                = {:<.16f}", energy_per_site);
    if(status.algo_type == AlgorithmType::xDMRG)
        tools::log->info("Energy density (rescaled 0 to 1) ε = {:<6.4f}",
                         tools::finite::measure::energy_normalized(tensors, status.energy_min_per_site, status.energy_max_per_site));
    double variance = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : std::log10(tools::finite::measure::energy_variance(tensors));
    tools::log->info("Energy variance log₁₀ σ²(E)        = {:<.16f}", variance);
    tools::log->info("Bond dimension maximum χmax        = {}", settings::chi_lim_max(status.algo_type));
    tools::log->info("Bond dimensions χ                  = {}", tools::finite::measure::bond_dimensions(*tensors.state));
    tools::log->info("Bond dimension  χ (mid)            = {}", tools::finite::measure::bond_dimension_midchain(*tensors.state));
    tools::log->info("Entanglement entropies Sₑ          = {:.5f}", fmt::join(tools::finite::measure::entanglement_entropies(*tensors.state), ", "));
    tools::log->info("Entanglement entropy   Sₑ (mid)    = {:.5f}", tools::finite::measure::entanglement_entropy_midchain(*tensors.state), ", ");
    if(status.algo_type == AlgorithmType::fLBIT) {
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