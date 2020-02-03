//
// Created by david on 2019-06-24.
//

#include "class_algorithm_finite.h"
#include <h5pp/h5pp.h>
#include <io/class_h5table_buffer.h>
#include <math/nmspc_math.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_finite.h>
#include <tools/nmspc_tools.h>

class_algorithm_finite::class_algorithm_finite(std::shared_ptr<h5pp::File> h5ppFile_, std::string sim_name, SimulationType sim_type, size_t num_sites)
    : class_algorithm_base(std::move(h5ppFile_), sim_name, sim_type)

{
    log->trace("Constructing class_algorithm_finite");
    state = std::make_unique<class_state_finite>();

    if(settings::output::storage_level >= StorageLevel::NORMAL) {
        log->trace("Constructing table buffers in finite base");
        h5tbuf_measurements = std::make_shared<class_h5table_buffer<class_h5table_measurements_finite>>(h5pp_file, sim_name + "/journal/measurements");
    }

    state->set_chi_lim(2); // Can't call chi_init() <-- it's a pure virtual function
    if(state->hasNaN()) throw std::runtime_error("State has NAN's before initializing it");

    tools::finite::mpo::initialize(*state, num_sites, settings::model::model_type);
    tools::finite::mps::initialize(*state, num_sites);
    tools::finite::mpo::randomize(*state);
    tools::finite::mps::randomize(*state, settings::model::initial_parity_sector, settings::model::state_number);
    tools::finite::mps::rebuild_environments(*state);

    tools::finite::debug::check_integrity(*state);

    S_mat.resize(state->get_length() + 1);
    X_mat.resize(state->get_length() + 1);

    tools::finite::print::print_hamiltonians(*state);
    tools::finite::print::print_state(*state);
    tools::finite::io::h5dset::write_model(*state, *h5pp_file, sim_name);
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
    if(not sim_on()) {
        return;
    }
    log->info("Starting {}", sim_name);
    tools::common::profile::t_tot.tic();
    if(h5pp_file and h5pp_file->getCreateMode() == h5pp::CreateMode::OPEN) {
        // This is case 1
        bool finOK_exists = h5pp_file->linkExists("common/finOK");
        bool mps_exists   = h5pp_file->linkExists(sim_name + "/state/mps");
        bool finOK        = false;
        if(finOK_exists) finOK = h5pp_file->readDataset<bool>("common/finOK");

        if(not finOK) {
            // Case 1 a -- run full simulation from scratch.
            log->trace("Case 1a");
            run_preprocessing();
            run_simulation();
        } else if(not mps_exists) {
            // Case 1 b
            log->trace("Case 1b");
            run_preprocessing();
            run_simulation();
        } else if(mps_exists) {
            // We can go ahead and load the state from output
            log->trace("Loading MPS from file");
            try {
                tools::finite::io::h5restore::load_from_hdf5(*h5pp_file, *state, sim_status, sim_name);
            } catch(std::exception &ex) {
                log->error("Failed to load from output: {}", ex.what());
                throw std::runtime_error("Failed to resume from file: " + std::string(ex.what()));
            } catch(...) {
                log->error("Unknown error when trying to resume from file.");
            }

            bool convergence_was_reached = h5pp_file->readDataset<bool>(sim_name + "/sim_status/simulation_has_converged");
            if(not convergence_was_reached) {
                // Case 1 c -- resume simulation, reset the number of sweeps first.
                log->trace("Case 1c");
                settings::xdmrg::max_sweeps += state->get_sweeps();
                run_simulation();

            } else {
                // Case 1 d -- not much else to do.. redo postprocessing for good measure.
                log->trace("Case 1d");
            }
        }
    } else {
        // This is case 2
        log->trace("Case 2");
        run_preprocessing();
        run_simulation();
    }
    tools::common::profile::t_tot.toc();
    run_postprocessing();
}

void class_algorithm_finite::run_preprocessing() {
    log->info("Running {} preprocessing (base)", sim_name);
    tools::common::profile::t_pre.tic();
    state->set_chi_max(chi_max());
    sim_status.chi_max = chi_max();
    update_bond_dimension_limit(chi_init());
    tools::common::profile::t_pre.toc();
    log->info("Finished {} preprocessing (base)", sim_name);
}

void class_algorithm_finite::single_DMRG_step(std::string ritz) {
    /*!
     * \fn void single_DMRG_step(std::string ritz)
     */
    log->trace("Starting single xDMRG step");
    tools::common::profile::t_sim.tic();
    Eigen::Tensor<Scalar, 4> theta = tools::finite::opt::find_ground_state(*state, ritz);
    tools::finite::opt::truncate_theta(theta, *state);
    state->clear_measurements();
    tools::common::profile::t_sim.toc();
    sim_status.wall_time = tools::common::profile::t_tot.get_age();
    sim_status.simu_time = tools::common::profile::t_sim.get_measured_time();
}

void class_algorithm_finite::run_postprocessing() {
    log->info("Running {} postprocessing", sim_name);
    tools::common::profile::t_pos.tic();
    tools::finite::debug::check_integrity(*state);
    state->clear_measurements();
    print_status_update();

    write_state(true);
    write_measurements(true);
    write_sim_status(true);
    write_profiling(true);

    auto state_projected = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::model::target_parity_sector);
    write_projection(state_projected, settings::model::target_parity_sector);

    print_status_full();
    tools::common::profile::t_pos.toc();
    tools::common::profile::print_profiling();
    log->info("Finished {} postprocessing", sim_name);
}

void class_algorithm_finite::move_center_point(std::optional<size_t> num_moves) {
    if(not num_moves.has_value()){
        if(state->active_sites.empty()) num_moves = 1ul;
        else if(settings::precision::move_sites_multidmrg == "one") num_moves = 1ul;
        else if(settings::precision::move_sites_multidmrg == "mid") num_moves = std::max(1ul, (state->active_sites.size() - 2ul)/2);
        else if(settings::precision::move_sites_multidmrg == "max") num_moves = std::max(1ul, state->active_sites.size() - 2ul);
        else{
            throw std::logic_error("Choose how many sites you should move! Expected one of {one,mid,max}, got [" + settings::precision::move_sites_multidmrg +"]");
        }
    }
    log->trace("Moving center point {} steps", num_moves.value());
    state->clear_cache();
    try {
        for(size_t i = 0; i < num_moves.value(); i++) {
            tools::finite::mps::move_center_point(*state);
            if(force_overlap_steps > 0) force_overlap_steps--;
        }
    } catch(std::exception &e) {
        tools::finite::print::print_state(*state);
        throw std::runtime_error("Failed to move center point: " + std::string(e.what()));
    }
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
        log->info("Setting initial bond dimension limit: {}", chi_init());
        state->set_chi_lim(chi_init());
        sim_status.chi_lim = chi_init();
        return;
    }

    sim_status.chi_lim_has_reached_chi_max = state->get_chi_lim() >= chi_max();
    if(not sim_status.chi_lim_has_reached_chi_max) {
        if(chi_grow()) {
            sim_status.chi_lim = state->get_chi_lim();

            // If we got here we want grow the bond dimension limit progressively during the simulation
            // Only increment the bond dimension if
            // * No experiments are on-going like perturbation or damping
            // * the simulation is stuck
            // * the state is limited by bond dimension
            if(not state->position_is_any_edge()) return;
            if(state->is_damped()) {
                log->info("State is undergoing disorder damping -- cannot increase bond dimension yet");
                return;
            }
            if(state->is_perturbed()) {
                log->info("State is undergoing perturbation -- cannot increase bond dimension yet");
                return;
            }
            if(not sim_status.simulation_has_got_stuck and state->get_chi_lim() >= 16){
                log->info("State is not stuck yet. Kept current limit {}", state->get_chi_lim());
                return;
            }

            if(log->level() <= spdlog::level::info){
                size_t trunc_bond_count  = state->num_sites_truncated(2 * settings::precision::svd_threshold);
                size_t bond_at_lim_count = state->num_bonds_at_limit();
                log->debug("Truncation errors     : {}", state->get_truncation_errors());
                log->debug("Bond dimensions       : {}", tools::finite::measure::bond_dimensions(*state));
                log->debug("Truncated bond count  : {} ", trunc_bond_count);
                log->debug("Bonds at limit  count : {} ", bond_at_lim_count);
                log->debug("Entanglement entropies: {} ", tools::finite::measure::entanglement_entropies(*state));
            }
            bool state_is_bond_limited = state->is_bond_limited(5 * settings::precision::svd_threshold);
            if(not state_is_bond_limited){
                log->info("State is not limited by its bond dimension. Kept current limit {}", state->get_chi_lim());
                return;
            }

            // Write current results before updating bond dimension
            write_state(true);
            write_measurements(true);
            write_sim_status(true);
            write_profiling(true);
            long chi_lim_new = std::min(state->get_chi_max(), state->get_chi_lim() * 2);
            log->info("Updating bond dimension limit {} -> {}", state->get_chi_lim(), chi_lim_new);
            state->set_chi_lim(chi_lim_new);
            clear_saturation_status();
            sim_status.chi_lim_has_reached_chi_max = state->get_chi_lim() == chi_max();
            if(sim_status.chi_lim_has_reached_chi_max and has_projected) has_projected = false;
            log->info("Projecting at site {} to direction {} after updating bond dimension to χ = {} ", state->get_position(),
                      settings::model::target_parity_sector, chi_lim_new);
            auto state_projected = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::model::target_parity_sector);
            write_projection(*state, settings::model::target_parity_sector);
            if(settings::model::projection_when_growing_chi) *state = state_projected;

            copy_from_tmp(true);
        } else {
            // Here the settings specify to just set the limit to maximum chi directly
            log->info("Setting bond dimension limit to maximum = {}", chi_max());
            state->set_chi_lim(chi_max());
        }
    } else {
        log->debug("Chi limit has reached max: {} -> Kept current bond dimension limit {}", chi_max(), state->get_chi_lim());
    }
    sim_status.chi_lim = state->get_chi_lim();
    if(state->get_chi_lim() > state->get_chi_max())
        throw std::runtime_error(fmt::format("chi_lim is larger than chi_max! {} > {}", state->get_chi_lim(), state->get_chi_max()));
}

void class_algorithm_finite::reset_to_initial_state() {
    log->trace("Resetting MPS to initial product state in parity sector: {}, state number {}", settings::model::initial_parity_sector,
               settings::model::state_number, settings::model::use_pauli_eigvecs);
    if(state->get_length() != (size_t) num_sites()) throw std::range_error("System size mismatch");
    // Initialize state
    tools::finite::mps::randomize(*state, settings::model::initial_parity_sector, settings::model::state_number, settings::model::use_pauli_eigvecs);
    clear_saturation_status();
    state->lowest_recorded_variance = 1;
    sim_status.iteration            = state->reset_sweeps();
}

void class_algorithm_finite::reset_to_random_state(const std::string &parity_sector) {
    log->trace("Resetting MPS to random product state in parity sector: {}", parity_sector);
    if(state->get_length() != (size_t) num_sites()) throw std::range_error("System size mismatch");
    // Randomize state
    tools::finite::mps::randomize(*state, parity_sector, -1, settings::model::use_pauli_eigvecs);
    clear_saturation_status();
    state->lowest_recorded_variance = 1;
    sim_status.iteration            = state->reset_sweeps();
}

void class_algorithm_finite::try_projection() {
    bool try_when_stuck =
        settings::model::projection_trial_when_stuck and sim_status.simulation_has_got_stuck and not has_projected and state->position_is_any_edge();

    bool try_every_sweep = settings::model::projection_on_every_sweep and sim_status.iteration >= 2 and state->position_is_any_edge();

    if(try_every_sweep or try_when_stuck) {
        log->debug("Trying projection to {}", settings::model::target_parity_sector);
        auto   state_projected    = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::model::target_parity_sector);
        double variance_original  = tools::finite::measure::energy_variance_per_site(*state);
        double variance_projected = tools::finite::measure::energy_variance_per_site(state_projected);

        has_projected = true;

        if(variance_projected < variance_original) {
            log->info("Projection: variance improved {:.8} -> {:.8}", std::log10(variance_original), std::log10(variance_projected));
            *state = state_projected;
        } else {
            log->info("Projection: variance would have worsened {:.8} -> {:.8}", std::log10(variance_original), std::log10(variance_projected));
        }
    }
}

void class_algorithm_finite::try_chi_quench() {
    if(not settings::model::quench_chi_when_stuck) return;
    if(not state->position_is_any_edge()) return;
    if(force_overlap_steps > 0) return;
    if(not sim_status.simulation_has_got_stuck) {
        tools::log->info("Chi quench skipped: simulation not stuck");
        return;
    }
    if(num_chi_quenches >= max_chi_quenches) {
        tools::log->info("Chi quench skipped: max number of chi quenches ({}) have been made already", num_chi_quenches);
        return;
    }
    size_t trunc_bond_count  = state->num_sites_truncated(5 * settings::precision::svd_threshold);
    size_t bond_at_lim_count = state->num_bonds_at_limit();
    log->debug("Truncation errors: {}", state->get_truncation_errors());
    log->debug("Bond dimensions  : {}", tools::finite::measure::bond_dimensions(*state));
    log->debug("Entanglement entr: {}", tools::finite::measure::entanglement_entropies(*state));
    log->debug("Truncated bond count: {} ", trunc_bond_count);
    log->debug("Bonds at limit  count: {} ", bond_at_lim_count);
    log->debug("threshold: {} ", std::pow(5 * settings::precision::svd_threshold,2));
    if(state->is_bond_limited(5 * settings::precision::svd_threshold)) {
        tools::log->info("Chi quench skipped: state is bond limited - prefer updating bond dimension");
        return;
    }

    tools::log->info("Chi quench started");
//    long smaller_chi_lim = std::min(16l, state->get_chi_lim() / 2);
//    tools::finite::mps::truncate_all_sites(*state, smaller_chi_lim, 1, 0);
    clear_saturation_status();
    force_overlap_steps = 1 * (state->get_length() - 2);
    num_chi_quenches++;
    if(force_overlap_steps > 2 * state->get_length() or force_overlap_steps < 0)
        throw std::runtime_error("Force overlap out of range: " + std::to_string(force_overlap_steps));
}

void class_algorithm_finite::check_convergence_variance(double threshold, double slope_threshold) {
    // Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    if(not state->position_is_any_edge()) {
        return;
    }
    log->debug("Checking convergence of variance mpo");
    threshold       = std::isnan(threshold) ? settings::precision::variance_convergence_threshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::variance_slope_threshold : slope_threshold;
    auto report = check_saturation_using_slope(V_mpo_vec, X_mpo_vec, tools::finite::measure::energy_variance(*state), sim_status.iteration, 1, slope_threshold);
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
        log->debug("Variance slope details:");
        log->debug(" -- relative slope    = {} %", report.slope);
        log->debug(" -- tolerance         = {} %", slope_threshold);
        log->debug(" -- last var average  = {} ", report.avgY);
        log->debug(" -- check from        = {} ", report.check_from);
        log->debug(" -- var history       = {} ", V_mpo_vec);
        log->debug(" -- slope history     = {} ", V_mpo_slopes);
        log->debug(" -- has saturated     = {} ", sim_status.variance_mpo_has_saturated);
        log->debug(" -- has saturated for = {} ", sim_status.variance_mpo_saturated_for);
        log->debug(" -- has converged     = {} ", sim_status.variance_mpo_has_converged);
        log->debug(" -- has converged for = {} ", converged_count);
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
    log->debug("Checking convergence of entanglement");

    slope_threshold                         = std::isnan(slope_threshold) ? settings::precision::entropy_slope_threshold : slope_threshold;
    auto                          entropies = tools::finite::measure::entanglement_entropies(*state);
    std::vector<SaturationReport> reports(entropies.size());

    for(size_t site = 0; site < entropies.size(); site++) {
        reports[site] = check_saturation_using_slope(S_mat[site], X_mat[site], entropies[site], sim_status.iteration, 1, slope_threshold);
    }
    bool all_computed                     = std::all_of(reports.begin(), reports.end(), [](const SaturationReport r) { return r.has_computed; });
    sim_status.entanglement_has_saturated = false;
    if(all_computed) {
        // idx_max_slope is the index to the site with maximum slope
        size_t idx_max_slope =
            std::distance(reports.begin(), std::max_element(reports.begin(), reports.end(),
                                                            [](const SaturationReport &r1, const SaturationReport &r2) { return r1.slope < r2.slope; }));
        // idx_max_slope is the index to the site with maximum slope
        //        size_t idx_min_satur = std::distance(reports.begin(),
        //                                             std::min_element(reports.begin(),reports.end(),
        //                                   [](const SaturationReport &r1, const SaturationReport &r2)
        //                                   {return r1.saturated_for < r2.saturated_for;}));

        S_slopes.push_back(reports[idx_max_slope].slope);
        auto   last_nonsaturated_ptr = std::find_if(S_slopes.rbegin(), S_slopes.rend(), [slope_threshold](auto const &val) { return val > slope_threshold; });
        size_t saturated_count       = (size_t) std::distance(S_slopes.rbegin(), last_nonsaturated_ptr);

        sim_status.entanglement_has_saturated = S_slopes.back() < slope_threshold;
        sim_status.entanglement_saturated_for = saturated_count;
        std::vector<double> all_avergs;
        std::vector<double> all_slopes;
        for(auto &r : reports) all_avergs.push_back(r.avgY);
        for(auto &r : reports) all_slopes.push_back(r.slope);
        log->debug("Max slope of entanglement entropy at site {}: {:.8f} %", idx_max_slope, S_slopes.back());
        log->debug("Entanglement slope details of worst slope:");
        log->debug(" -- site              = {}", idx_max_slope);
        log->debug(" -- relative slope    = {} %", reports[idx_max_slope].slope);
        log->debug(" -- tolerance         = {} %", slope_threshold);
        log->debug(" -- check from        = {} ", reports[idx_max_slope].check_from);
        log->debug(" -- ent history       = {} ", S_mat[idx_max_slope]);
        log->debug(" -- slope history     = {} ", S_slopes);
        log->debug(" -- has saturated     = {} ", sim_status.entanglement_has_saturated);
        log->debug(" -- has saturated for = {} (site {} )", sim_status.entanglement_saturated_for, saturated_count);
        log->debug(" -- all averages      = {} ", all_avergs);
        log->debug(" -- all slopes        = {} ", all_slopes);
        if(reports[idx_max_slope].slope > slope_threshold and sim_status.entanglement_has_saturated) throw std::logic_error("Not supposed to be saturated!!");
    }
    sim_status.entanglement_has_converged = sim_status.entanglement_has_saturated;
}

void class_algorithm_finite::clear_saturation_status() {
    log->trace("Clearing saturation status");
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

void class_algorithm_finite::write_state(bool result) {
    if(settings::output::storage_level == StorageLevel::NONE) {
        return;
    }
    // There are two places to write states, in "/results" and "/journal".
    // "/journal" is written once per sweep, and "/results" only when there has been
    // convergence, or updating bond dimension, or similar.
    // The idea is to keep a more detailed journal if storage level is >= NORMAL,
    // but the main results should be written always
    if(result) {
        // This means that we are writing an important result:
        // Either the simulation has converged successfully or
        // it has finalized some stage, like saturated at the
        // current bond dimension.

        tools::finite::io::h5dset::write_all_state(*state, *h5pp_file, sim_name);
        if(store_wave_function()) {
            //  Write the wavefunction (this is only defined for short enough state ( L < 14 say)
            h5pp_file->writeDataset(tools::finite::measure::mps_wavefn(*state), sim_name + "/state/psi");
        }
        if(settings::output::storage_level >= StorageLevel::FULL and chi_grow()) {
            // In full mode we keep a copy every time we update
            std::string prefix = sim_name + "/results/iter_" + std::to_string(sim_status.iteration) + "_step_" + std::to_string(sim_status.step) + "_chi_" +
                                 std::to_string(sim_status.chi_lim);
            tools::finite::io::h5dset::write_all_state(*state, *h5pp_file, prefix);
            if(store_wave_function()) {
                //  Write the wavefunction (this is only defined for short enough state ( L < 14 say)
                h5pp_file->writeDataset(tools::finite::measure::mps_wavefn(*state), prefix + "/state/psi");
            }
        }
    }

    if(not state->position_is_any_edge()) {
        return;
    }
    if(math::mod(sim_status.iteration, write_freq()) != 0) {
        return;
    } // Check that we write according to the frequency given
    tools::finite::io::h5dset::write_all_state(*state, *h5pp_file, sim_name);
    if(settings::output::storage_level >= StorageLevel::FULL) {
        std::string prefix = sim_name + "/journal/iter_" + std::to_string(sim_status.iteration);
        tools::finite::io::h5dset::write_all_state(*state, *h5pp_file, prefix);
    }
}

void class_algorithm_finite::write_measurements(bool result) {
    if(not h5pp_file) return;
    if(settings::output::storage_level == StorageLevel::NONE) return;

    // There are two places to write measurements, in "/results" and "/journal".
    // "/journal" is written once per sweep, and "/results" only when there has been
    // convergence, or updating bond dimension, or similar.
    // The idea is to keep a more detailed journal if storage level is >= NORMAL,
    // but the main results should be written always, if we're writing anything,
    // and to write it immediately, without buffer.

    if(result) {
        // This means that we are writing an important result:
        // Either the simulation has converged successfully or
        // it has finalized some stage, like saturated at the
        // current bond dimension.
        class_h5table_buffer<class_h5table_measurements_finite> h5tbuf_measurements_results(h5pp_file, sim_name + "/results/measurements");
        tools::finite::io::h5table::write_measurements(*state, sim_status, h5tbuf_measurements_results);
        tools::finite::io::h5dset::write_array_measurements(*state, *h5pp_file, sim_name + "/results");
        if(settings::output::storage_level >= StorageLevel::NORMAL and chi_grow()) {
            // Write even more results
            std::string prefix = sim_name + "/results/chi_" + std::to_string(sim_status.chi_lim);
            tools::finite::io::h5dset::write_array_measurements(*state, *h5pp_file, prefix);
        }
    }

    // For the journalning we do the usual checks
    // Except we can forcefully insert an entry to the journal if it's a result
    // and we're actually keeping a journal. We do no journaling on LIGHT mode
    if(h5tbuf_measurements == nullptr) {
        return;
    }
    if(settings::output::storage_level <= StorageLevel::LIGHT) {
        return;
    }
    if(result and not state->position_is_any_edge()) {
        // To avoid duplicate lines, we write results to journal if we are not on an edge.
        // At the edge we call this function with result = false anyway
        tools::finite::io::h5table::write_measurements(*state, sim_status, *h5tbuf_measurements);
        tools::finite::io::h5dset::write_array_measurements(*state, *h5pp_file, sim_name + "/journal");

    } else {
        // This is the most common journaling. This is not a result, just an entry in the journal
        if(not state->position_is_any_edge()) {
            return;
        }
        if(math::mod(sim_status.iteration, write_freq()) != 0) {
            return;
        } // Check that we write according to the frequency given
        tools::finite::io::h5table::write_measurements(*state, sim_status, *h5tbuf_measurements);
        tools::finite::io::h5dset::write_array_measurements(*state, *h5pp_file, sim_name + "/journal");
        if(settings::output::storage_level >= StorageLevel::FULL) {
            // In full mode we keep everything at each sweep
            std::string prefix = sim_name + "/journal/iter_" + std::to_string(sim_status.iteration);
            tools::finite::io::h5dset::write_array_measurements(*state, *h5pp_file, prefix);
        }
    }
}

void class_algorithm_finite::write_sim_status(bool result) {
    if(settings::output::storage_level == StorageLevel::NONE) {
        return;
    }
    // There are two places to write simulation status, in "/results" and "/journal".
    // "/journal" is written once per sweep, and "/results" only when there has been
    // convergence, or updating bond dimension, or similar.
    // The idea is to keep a more detailed journal if storage level is >= NORMAL,
    // but the main results should be written always, if we're writing anything,
    // and to write it immediately, without buffer.
    if(result) {
        // This means that we are writing an important result:
        // Either the simulation has converged successfully or
        // it has finalized some stage, like saturated at the
        // current bond dimension.
        class_h5table_buffer<class_h5table_simulation_status> h5tbuf_sim_status_results(h5pp_file, sim_name + "/results/sim_status");
        tools::finite::io::h5table::write_sim_status(sim_status, h5tbuf_sim_status_results);
    }

    // For the journalning we do the usual checks
    // Except we can forcefully insert an entry to the journal if it's a result
    // and we're actually keeping a journal. We do no journaling on LIGHT mode
    if(h5tbuf_sim_status == nullptr) {
        return;
    }
    if(settings::output::storage_level <= StorageLevel::LIGHT) {
        return;
    }
    if(result and not state->position_is_any_edge()) {
        // To avoid duplicate lines, we write results to journal if we are not on an edge.
        // At the edge we call this function with result = false anyway
        tools::finite::io::h5table::write_sim_status(sim_status, *h5tbuf_sim_status);
    } else {
        // This is the most common journaling. This is not a result, just an entry in the journal
        if(not state->position_is_any_edge()) {
            return;
        }
        if(math::mod(sim_status.iteration, write_freq()) != 0) {
            return;
        } // Check that we write according to the frequency given
        tools::finite::io::h5table::write_sim_status(sim_status, *h5tbuf_sim_status);
    }
}

void class_algorithm_finite::write_profiling(bool result) {
    if(not settings::profiling::on) return;
    if(settings::output::storage_level == StorageLevel::NONE) return;
    if(h5tbuf_profiling == nullptr) return;
    if(not result) {
        if(not state->position_is_any_edge()) return;
        if(settings::output::storage_level <= StorageLevel::LIGHT) return;
        if(math::mod(sim_status.iteration, write_freq()) != 0) return; // Check that we write according to the frequency given
    }
    tools::finite::io::h5table::write_profiling(sim_status, *h5tbuf_profiling);
}

void class_algorithm_finite::write_projection(const class_state_finite &state_projected, std::string parity_sector) {
    if(settings::output::storage_level == StorageLevel::NONE) {
        return;
    }
    if(parity_sector == "none") return;
    std::string prefix = sim_name + "/projections/" + parity_sector;
    tools::finite::io::h5dset::write_all_state(state_projected, *h5pp_file, prefix);
    class_h5table_buffer<class_h5table_measurements_finite> h5tbuf_measurements_projection(h5pp_file, prefix + "/results/measurements");
    class_h5table_buffer<class_h5table_simulation_status>   h5tbuf_sim_status_projection(h5pp_file, prefix + "/results/sim_status");
    tools::finite::io::h5table::write_sim_status(sim_status, h5tbuf_sim_status_projection);
    tools::finite::io::h5table::write_measurements(state_projected, sim_status, h5tbuf_measurements_projection);
    tools::finite::io::h5dset::write_array_measurements(state_projected, *h5pp_file, prefix + "/results");
}

void class_algorithm_finite::copy_from_tmp(bool result) {
    if(settings::output::storage_level == StorageLevel::NONE) return;
    if(result) tools::common::io::h5tmp::copy_from_tmp(h5pp_file->getFilePath());
    if(not state->position_is_any_edge()) return;
    if(math::mod(sim_status.iteration, settings::output::copy_from_temp_freq) != 0) return; // Check that we write according to the frequency given
    tools::common::io::h5tmp::copy_from_tmp(h5pp_file->getFilePath());
}

void class_algorithm_finite::print_status_update() {
    if(math::mod(sim_status.step, print_freq()) != 0) return;
    if(print_freq() == 0) return;

    using namespace tools::finite::measure;
    std::string report;
    report += fmt::format("{:<} ", sim_name);
    report += fmt::format("iter: {:<4} ", sim_status.iteration);
    report += fmt::format("step: {:<5} ", sim_status.step);
    report += fmt::format("L: {} l: {:<2} ", state->get_length(), state->get_position());
    report += fmt::format("E/L: {:<20.16f} ", energy_per_site(*state));
    if(sim_type == SimulationType::xDMRG) {
        report + fmt::format("ε: {:<6.4f} ", sim_status.energy_dens);
    }
    report += fmt::format("Sₑ(l): {:<10.8f} ", entanglement_entropy_current(*state));
    report += fmt::format("log₁₀ σ²(E)/L: {:<10.6f} [{:<10.6f}] ", std::log10(energy_variance_per_site(*state)),
                         std::log10(state->lowest_recorded_variance / state->get_length()));
    report += fmt::format("χmax: {:<3} χlim: {:<3} χ: {:<3} ", chi_max(), state->get_chi_lim(), bond_dimension_current(*state));
    report += fmt::format("log₁₀ trunc: {:<10.4f} ", std::log10(state->get_truncation_error(state->get_position())));
    report += fmt::format("stk: {:<1} ", sim_status.simulation_has_stuck_for);
    report += fmt::format("sat: [σ² {:<1} Sₑ {:<1}] ", sim_status.variance_mpo_saturated_for, sim_status.entanglement_saturated_for);
    report += fmt::format("con: {:<5} ", sim_status.simulation_has_converged);
    report += fmt::format("time: {:<8.2f}s ", tools::common::profile::t_tot.get_age());
    report += fmt::format("mem MB: [Rss {:<.1f} Peak {:<.1f} Vm {:<.1f}] ", process_memory_in_mb("VmRSS"), process_memory_in_mb("VmHWM"),
                         process_memory_in_mb("VmPeak"));
    log->info(report);
}

void class_algorithm_finite::print_status_full() {
    using namespace tools::finite::measure;
    log->info("{:=^60}", "");
    log->info("= {: ^56} =", "Final results [" + sim_name + "]");
    log->info("{:=^60}", "");
    log->info("--- Final results  --- {} ---", sim_name);
    log->info("Sites                              = {}", state->get_length());
    log->info("Iterations (sweeps)                = {}", sim_status.iteration);
    log->info("Steps                              = {}", sim_status.step);
    log->info("Simulation time                    = {:<.1f} s = {:<.2f} min", tools::common::profile::t_tot.get_age(),
              tools::common::profile::t_tot.get_age() / 60);
    log->info("Energy per site E/L                = {:<.16f}", energy_per_site(*state));
    if(sim_type == SimulationType::xDMRG) {
        log->info("Energy density (rescaled 0 to 1) ε = {:<6.4f}", energy_normalized(*state, sim_status));
    }
    log->info("Variance per site log₁₀ σ²(E)/L    = {:<.16f}", std::log10(energy_variance_per_site(*state)));
    log->info("Bond dimension maximum χmax        = {}", chi_max());
    log->info("Bond dimensions χ                  = {}", bond_dimensions(*state));
    log->info("Entanglement entropies Sₑ          = {}", entanglement_entropies(*state));
    log->info("Truncation Errors                  = {}", state->get_truncation_errors());
    log->info("Simulation converged               = {:<}", sim_status.simulation_has_converged);
    log->info("Simulation saturated               = {:<}", sim_status.simulation_has_saturated);
    log->info("Simulation succeeded               = {:<}", sim_status.simulation_has_succeeded);
    log->info("Simulation got stuck               = {:<}", sim_status.simulation_has_got_stuck);
    log->info("σ² slope                           = {:<8.4f} %   Converged : {:<8}  Saturated: {:<8}", V_mpo_slopes.back(),
              sim_status.variance_mpo_has_converged, sim_status.variance_mpo_has_saturated);
    log->info("Sₑ slope                           = {:<8.4f} %   Converged : {:<8}  Saturated: {:<8}", S_slopes.back(), sim_status.entanglement_has_converged,
              sim_status.entanglement_has_saturated);
    log->info("Memory RSS                         = {:<.1f} MB", process_memory_in_mb("VmRSS"));
    log->info("Memory Peak                        = {:<.1f} MB", process_memory_in_mb("VmHWM"));
    log->info("Memory Vm                          = {:<.1f} MB", process_memory_in_mb("VmPeak"));
}
