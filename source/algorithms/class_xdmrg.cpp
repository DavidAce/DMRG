//
// Created by david on 2018-02-09.
//

#include "class_xdmrg.h"
#include "class_fdmrg.h"
#include <config/nmspc_settings.h>
#include <math/num.h>
#include <math/rnd.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/debug.h>
#include <tools/finite/io.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt.h>
#include <tools/finite/opt_tensor.h>

class_xdmrg::class_xdmrg(std::shared_ptr<h5pp::File> h5ppFile_) : class_algorithm_finite(std::move(h5ppFile_), AlgorithmType::xDMRG) {
    tools::log->trace("Constructing class_xdmrg");
}

void class_xdmrg::resume() {
    // Resume can imply many things
    // 1) Resume a simulation which terminated prematurely
    // 2) Resume a previously successful simulation. This may be desireable if the config
    //    wants something that is not present in the file.
    //      a) A certain number of states
    //      b) A state inside of a particular energy window
    //      c) The ground or "roof" states
    // To guide the behavior, we check the setting ResumePolicy.

    auto state_prefix = tools::common::io::h5resume::find_resumable_state(*h5pp_file, algo_type);
    if(state_prefix.empty()) throw std::runtime_error("Could not resume: no valid state candidates found for resume");
    tools::log->info("Resuming state [{}]", state_prefix);
    tools::finite::io::h5resume::load_tensors(*h5pp_file, state_prefix, tensors, status);

    // Our first task is to decide on a state name for the newly loaded state
    // The simplest is to inferr it from the state prefix itself
    auto name   = tools::common::io::h5resume::extract_state_name(state_prefix);
    auto number = tools::common::io::h5resume::extract_state_number(state_prefix);
    if(number) excited_state_number = number.value();

    // Initialize a custom task list
    std::list<xdmrg_task> task_list;

    if(not status.algorithm_has_finished) {
        // This could be a checkpoint state
        // Simply "continue" the algorithm until convergence
        task_list = {xdmrg_task::FIND_EXCITED_STATE, xdmrg_task::POST_DEFAULT};
        run_task_list(task_list);
    }

    // If we reached this point the current state has finished for one reason or another.
    // We may still have some more things to do, e.g. the config may be asking for more states
    // Note that if max_states = 4 we are asking for 4 unbiased states, the 0th does not count (it's a seed!)
    // Example:
    //      max_states = 4
    //      excited_state_number = 1
    //      missing_state_number = 4 - 1 = 3
    auto missing_state_number = settings::xdmrg::max_states >= excited_state_number ? settings::xdmrg::max_states - excited_state_number : 0;
    for(size_t new_state_num = 0; new_state_num < missing_state_number; new_state_num++) {
        task_list.emplace_back(xdmrg_task::NEXT_RANDOMIZE_INTO_ENTANGLED_STATE);
        task_list.emplace_back(xdmrg_task::FIND_EXCITED_STATE);
        task_list.emplace_back(xdmrg_task::POST_DEFAULT);
    }
}

void class_xdmrg::run_default_task_list() {
    std::list<xdmrg_task> default_task_list = {
        xdmrg_task::INIT_DEFAULT,
        xdmrg_task::FIND_EXCITED_STATE,
        xdmrg_task::POST_DEFAULT,
    };

    // Insert requested number of excited states
    for(size_t num = 1; num < settings::xdmrg::max_states; num++) {
        default_task_list.emplace_back(xdmrg_task::NEXT_RANDOMIZE_FROM_CURRENT_STATE);
        default_task_list.emplace_back(xdmrg_task::FIND_EXCITED_STATE);
        default_task_list.emplace_back(xdmrg_task::POST_DEFAULT);
    }

    run_task_list(default_task_list);
    if(not default_task_list.empty()) {
        for(auto &task : default_task_list) tools::log->critical("Unfinished task: {}", enum2str(task));
        throw std::runtime_error("Simulation ended with unfinished tasks");
    }
}

void class_xdmrg::run_task_list(std::list<xdmrg_task> &task_list) {
    while(not task_list.empty()) {
        auto task = task_list.front();
        switch(task) {
            case xdmrg_task::INIT_RANDOMIZE_MODEL: randomize_model(); break;
            case xdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE: randomize_state(ResetReason::INIT, StateType::RANDOM_PRODUCT_STATE); break;
            case xdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE: randomize_state(ResetReason::INIT, StateType::RANDOM_ENTANGLED_STATE); break;
            case xdmrg_task::INIT_RANDOMIZE_FROM_CURRENT_STATE: randomize_state(ResetReason::INIT, StateType::RANDOMIZE_GIVEN_STATE); break;
            case xdmrg_task::INIT_RANDOMIZE_INTO_STATE_IN_WIN:
                randomize_into_state_in_energy_window(ResetReason::INIT, settings::strategy::initial_state);
                break;
            case xdmrg_task::INIT_BOND_DIM_LIMITS: init_bond_dimension_limits(); break;
            case xdmrg_task::INIT_ENERGY_LIMITS: init_energy_limits(); break;
            case xdmrg_task::INIT_WRITE_MODEL: write_to_file(StorageReason::MODEL); break;
            case xdmrg_task::INIT_CLEAR_STATUS: status.clear(); break;
            case xdmrg_task::INIT_DEFAULT:
                run_preprocessing();
                break;
                //            case xdmrg_task::NEXT_TRUNCATE_ALL_SITES: truncate_all_sites(); break;
            case xdmrg_task::NEXT_RANDOMIZE_INTO_PRODUCT_STATE: randomize_state(ResetReason::NEW_STATE, StateType::RANDOM_PRODUCT_STATE); break;
            case xdmrg_task::NEXT_RANDOMIZE_INTO_ENTANGLED_STATE: randomize_state(ResetReason::NEW_STATE, StateType::RANDOM_ENTANGLED_STATE); break;
            case xdmrg_task::NEXT_RANDOMIZE_FROM_CURRENT_STATE: randomize_state(ResetReason::NEW_STATE, StateType::RANDOMIZE_GIVEN_STATE); break;
            case xdmrg_task::NEXT_RANDOMIZE_INTO_STATE_IN_WIN: randomize_state(ResetReason::NEW_STATE, settings::strategy::initial_state); break;
            case xdmrg_task::FIND_ENERGY_RANGE: find_energy_range(); break;
            case xdmrg_task::FIND_EXCITED_STATE:
                state_name = fmt::format("state_{}", excited_state_number);
                run_algorithm();
                break;
            case xdmrg_task::POST_WRITE_RESULT: write_to_file(StorageReason::FINISHED); break;
            case xdmrg_task::POST_PRINT_RESULT: print_status_full(); break;
            case xdmrg_task::POST_DEFAULT: run_postprocessing(); break;
        }
        task_list.pop_front();
    }
}

void class_xdmrg::init_energy_limits(std::optional<double> energy_density_target, std::optional<double> energy_density_window) {
    if(not energy_density_target) energy_density_target = settings::xdmrg::energy_density_target;
    if(not energy_density_window) energy_density_window = settings::xdmrg::energy_density_window;
    if(energy_density_target.value() < 0.0 or energy_density_target.value() > 1.0)
        throw std::runtime_error(
            fmt::format("Error setting energy density target: Expected value in range [0 - 1.0], got: [{:.8f}]", energy_density_target.value()));
    if(energy_density_window.value() < 0.0 or energy_density_window.value() > 0.5)
        throw std::runtime_error(
            fmt::format("Error setting energy density window: Expected value in range [0 - 0.5], got: [{:.8f}]", energy_density_window.value()));
    status.energy_dens_target = energy_density_target.value();
    status.energy_dens_window = energy_density_window.value();

    // Set energy boundaries. This function is supposed to run after find_energy_range!
    if(status.energy_max_per_site == status.energy_min_per_site)
        throw std::runtime_error(fmt::format("Could not set energy limits because energy_max_per_site == {} and energy_min_per_site == {}\n"
                                             "Try running find_energy_range() first",
                                             status.energy_max_per_site, status.energy_min_per_site));
    status.energy_tgt_per_site  = status.energy_min_per_site + status.energy_dens_target * (status.energy_max_per_site - status.energy_min_per_site);
    status.energy_ulim_per_site = status.energy_tgt_per_site + status.energy_dens_window * (status.energy_max_per_site - status.energy_min_per_site);
    status.energy_llim_per_site = status.energy_tgt_per_site - status.energy_dens_window * (status.energy_max_per_site - status.energy_min_per_site);
    tools::log->info("Energy minimum     (per site) = {:.8f}", status.energy_min_per_site);
    tools::log->info("Energy maximum     (per site) = {:.8f}", status.energy_max_per_site);
    tools::log->info("Energy target      (per site) = {:.8f}", status.energy_tgt_per_site);
    tools::log->info("Energy lower limit (per site) = {:.8f}", status.energy_llim_per_site);
    tools::log->info("Energy upper limit (per site) = {:.8f}", status.energy_ulim_per_site);
}

void class_xdmrg::run_preprocessing() {
    tools::log->info("Running {} preprocessing", algo_name);
    tools::common::profile::t_pre->tic();
    status.clear();
    randomize_model(); // First use of random!
    init_bond_dimension_limits();
    find_energy_range();
    init_energy_limits();
    if(settings::xdmrg::energy_density_window != 0.5) randomize_into_state_in_energy_window(ResetReason::INIT, settings::strategy::initial_state);
    else
        randomize_state(ResetReason::INIT, settings::strategy::initial_state);
    write_to_file(StorageReason::MODEL);
    tools::common::profile::t_pre->toc();
    tools::log->info("Finished {} preprocessing", algo_name);
}

void class_xdmrg::run_algorithm() {
    if(state_name.empty()) fmt::format("state_{}", excited_state_number);
    tools::common::profile::reset_for_run_algorithm();
    tools::log->info("Starting {} simulation of model [{}] for state [{}]", algo_name, enum2str(settings::model::model_type), state_name);
    tools::common::profile::t_sim->tic();
    while(true) {
        single_xDMRG_step();
        print_status_update();
        write_to_file();
        check_convergence();
        update_bond_dimension_limit(); // Will update bond dimension if the state precision is being limited by bond dimension

        try_projection();
        try_bond_dimension_quench();
        try_disorder_damping();
        try_hamiltonian_perturbation();
        reduce_mpo_energy();

        // It's important not to perform the last move.
        // That last state would not get optimized
        if(tensors.state->position_is_any_edge() and not tensors.model->is_perturbed() and not tensors.model->is_damped()) {
            if(status.iter >= settings::xdmrg::max_iters) {
                stop_reason = StopReason::MAX_ITERS;
                break;
            }
            if(status.algorithm_has_succeeded) {
                stop_reason = StopReason::SUCCEEDED;
                break;
            }
            if(status.algorithm_has_to_stop) {
                stop_reason = StopReason::SATURATED;
                break;
            }
            if(status.num_resets > settings::strategy::max_resets) {
                stop_reason = StopReason::MAX_RESET;
                break;
            }
            if(settings::strategy::randomize_early and excited_state_number == 0 and tensors.state->find_largest_chi() >= 32 and
               tools::finite::measure::energy_variance(tensors) < 1e-4) {
                stop_reason = StopReason::RANDOMIZE;
                break;
            }
        }
        // Update record holder
        if(tensors.position_is_any_edge() or tensors.measurements.energy_variance_per_site) {
            tools::log->trace("Updating variance record holder");
            auto var = tools::finite::measure::energy_variance_per_site(tensors);
            if(var < status.lowest_recorded_variance_per_site) status.lowest_recorded_variance_per_site = var;
        }
        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);
        move_center_point();
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", algo_name, state_name, enum2str(stop_reason));
    status.algorithm_has_finished = true;
    tools::common::profile::t_sim->toc();
}

void class_xdmrg::single_xDMRG_step() {
    tools::log->debug("Starting xDMRG step {} | iter {} | pos {} | dir {}", status.step, status.iter, status.position, status.direction);

    using namespace tools::finite;
    using namespace tools::finite::opt;

    //    IDEA:
    //        Try converging to some state from a product state.
    //        We know this is a biased state because it was generated  from a product state of low entanglement.
    //        After converging, apply randomly either identity or pauli_x MPO operators on some sites, to generate
    //        a new initial guess for some other state. Make sure the parity sector is honored. Converging again now
    //        should erase the any biasing that may have occurred due to the selection of initial product state.

    // Setup optpmization mode and space.
    //      - Start from the most general condition
    //      - Override later with narrower conditions

    // Set the fastest mode by default
    opt::OptMode  optMode  = opt::OptMode::VARIANCE;
    opt::OptSpace optSpace = opt::OptSpace::DIRECT;
    opt::OptType  optType  = opt::OptType::CPLX;

    // The first decision is easy. Real or complex optimization
    if(tensors.is_real()) optType = OptType::REAL;

    if(status.chi_lim <= 12) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

    if(status.chi_lim <= 8 or status.iter < 2 or status.step < 2*tensors.get_length()) {
        // TODO: We may not want to do this on states post randomization
        optMode  = OptMode::OVERLAP;
        optSpace = OptSpace::SUBSPACE_ONLY;
    }

    // Setup strong overrides to normal conditions, e.g.,
    //      - for experiments like perturbation or chi quench
    //      - when the algorithm has already converged

    if(status.variance_mpo_has_converged) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::DIRECT;
    }

    // Setup the maximum problem size
    long max_problem_size = 0;
    switch(optSpace) {
        case OptSpace::DIRECT: max_problem_size = settings::precision::max_size_direct; break;
        case OptSpace::SUBSPACE_ONLY:
        case OptSpace::SUBSPACE_AND_DIRECT: max_problem_size = settings::precision::max_size_part_diag; break;
    }

    // Make sure not to use SUBSPACE if the 2-site problem size is huge
    // When this happens, we can't use SUBSPACE even with two sites,
    // so we might as well give it to DIRECT, which handles larger problems.
    if(tensors.state->size_2site() > max_problem_size) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::DIRECT;
    }

    // We can optionally make many optimization trials with different
    // number of sites. I.e., if a 2 site optimization did not give
    // any significant improvement we may try with 4 or 8 sites.
    // Here we define a list of trials, where each entry determines
    // how many sites to try.
    std::vector<size_t> max_num_sites_list;
    if(status.algorithm_has_stuck_for > 0) max_num_sites_list = {4};
    else if(status.algorithm_has_stuck_for > 1 or chi_quench_steps > 0)
        max_num_sites_list = {4, settings::strategy::multisite_max_sites};
    else
        max_num_sites_list = {2};
    // Make sure the site list is sane by sorting and filtering out repeated/invalid entries.
    std::sort(max_num_sites_list.begin(), max_num_sites_list.end());
    std::unique(max_num_sites_list.begin(), max_num_sites_list.end());
    std::remove_if(max_num_sites_list.begin(), max_num_sites_list.end(), [](auto &elem) { return elem > settings::strategy::multisite_max_sites; });

    if(max_num_sites_list.empty()) max_num_sites_list = {2};
    if(max_num_sites_list.empty()) throw std::runtime_error("No sites selected for multisite xDMRG");

    tools::log->debug("Possible multisite step sizes: {}", max_num_sites_list);
    double                  variance_old_per_site = 1;
    std::vector<opt_tensor> results;
    for(auto &max_num_sites : max_num_sites_list) {
        if(optMode == opt::OptMode::OVERLAP and optSpace == opt::OptSpace::DIRECT) throw std::logic_error("OVERLAP mode and DIRECT space are incompatible");

        auto old_num_sites = tensors.active_sites.size();
        auto old_prob_size = tensors.active_problem_size();
        tensors.activate_sites(max_problem_size, max_num_sites);

        // Check that we are not about to solve the same problem again
        if(not results.empty() and tensors.active_sites.size() == old_num_sites and tensors.active_problem_size() == old_prob_size) {
            // If we reached this point we have exhausted the number of sites available
            tools::log->debug("Can't activate more sites");
            break;
        }
        variance_old_per_site = measure::energy_variance_per_site(tensors); // Should just take value from cache
        results.emplace_back(opt::find_excited_state(tensors, status, optMode, optSpace, optType));

        // We can now decide if we are happy with the result or not.
        double percent = 100 * results.back().get_variance_per_site() / variance_old_per_site;
        if(percent < 99) {
            tools::log->debug("State improved during {} optimization. Variance: new/old: {:.4f} %", optSpace, percent);
            break;
        } else {
            tools::log->debug("State did not improve during {} optimization. Variance new/old: {:.4f} %", optSpace, percent);
            continue;
        }
    }
    // Sort the results in order of increasing variance
    std::sort(results.begin(), results.end(), [](const opt_tensor &lhs, const opt_tensor &rhs) { return lhs.get_variance() < rhs.get_variance(); });

    for(auto &candidate : results)
        tools::log->debug("Candidate: {} | sites [{:>2}-{:<2}] | variance {:.16f} | energy {:.16f} | overlap {:.16f} | norm {:.16f} | time {:.4f} ms",
                          candidate.get_name(), candidate.get_sites().front(), candidate.get_sites().back(), std::log10(candidate.get_variance_per_site()),
                          candidate.get_energy_per_site(), candidate.get_overlap(), candidate.get_norm(), 1000 * candidate.get_time());

    // Take the best result
    const auto &winner   = results.front();
    tensors.active_sites = results.front().get_sites();
    tensors.sync_active_sites();
    if(std::log10(variance_old_per_site) / std::log10(winner.get_variance_per_site()) < 0.999) tensors.state->tag_active_sites_have_been_updated(true);

    // Truncate even more if doing chi quench
    //   if(chi_quench_steps > 0) chi_lim = chi_lim_quench_trail;

    // Do the truncation with SVD
    tensors.merge_multisite_tensor(winner.get_tensor(), status.chi_lim);
    if constexpr(settings::debug) {
        auto variance_per_site_before_svd = winner.get_variance_per_site();
        auto variance_per_site_after_svd  = tools::finite::measure::energy_variance_per_site(tensors);
        tools::log->trace("Variance check before SVD: {:.16f}", std::log10(variance_per_site_before_svd));
        tools::log->trace("Variance check after  SVD: {:.16f}", std::log10(variance_per_site_after_svd));
        tools::log->debug("Variance loss due to  SVD: {:.16f}", (variance_per_site_after_svd - variance_per_site_before_svd) / variance_per_site_after_svd);
    }

    debug::check_integrity(*tensors.state);
    status.energy_dens =
        (tools::finite::measure::energy_per_site(tensors) - status.energy_min_per_site) / (status.energy_max_per_site - status.energy_min_per_site);

    status.wall_time = tools::common::profile::t_tot->get_measured_time();
    status.algo_time = tools::common::profile::t_sim->get_measured_time();
}

void class_xdmrg::check_convergence() {
    tools::common::profile::t_con->tic();
    if(tensors.state->position_is_any_edge()) {
        check_convergence_variance();
        check_convergence_entg_entropy();
    }

    // TODO: Move this reset block away from here

    //    status.energy_dens = (tools::finite::measure::energy_per_site(*tensors.state) - status.energy_min_per_site ) / (status.energy_max_per_site -
    //    status.energy_min_per_site);
    bool outside_of_window = std::abs(status.energy_dens - status.energy_dens_target) > status.energy_dens_window;
    if(status.iter > 2 and tensors.state->position_is_any_edge()) {
        if(outside_of_window and
           (status.variance_mpo_has_saturated or status.variance_mpo_has_converged or tools::finite::measure::energy_variance_per_site(tensors) < 1e-4)) {
            double      old_energy_dens_window = status.energy_dens_window;
            double      new_energy_dens_window = std::min(energy_window_growth_factor * status.energy_dens_window, 0.5);
            std::string reason = fmt::format("energy {:.16f} saturated outside of energy window {} ± {}", tools::finite::measure::energy_per_site(tensors),
                                             status.energy_dens_target, status.energy_dens_window);
            tools::log->info("Increasing energy window: {} --> {}", old_energy_dens_window, new_energy_dens_window);
            status.energy_dens_window = new_energy_dens_window;
            randomize_into_state_in_energy_window(ResetReason::SATURATED, settings::strategy::initial_state, settings::strategy::target_sector);
        }
        //        else
        //        if( not     tensors.state->all_sites_updated()
        //            and     status.algorithm_has_got_stuck
        //            and     tools::finite::measure::energy_variance_per_site(*tensors.state) > 1e-4)
        //        {
        //            status.energy_dens_window = std::min(energy_window_growth_factor*status.energy_dens_window, 0.5);
        //            std::string reason = fmt::format("could not update all sites. Energy density: {}, Energy window: {} --> {}",
        //                     status.energy_dens, status.energy_dens_window, std::min(energy_window_growth_factor*status.energy_dens_window, 0.5)
        //                     );
        //                     randomize_into_state_in_energy_window(ResetReason::SATURATED,settings::strategy::target_sector);
        //        }
    }

    status.algorithm_has_converged = status.variance_mpo_has_converged and status.entanglement_has_converged;
    status.algorithm_has_saturated = status.variance_mpo_saturated_for >= min_saturation_iters and status.entanglement_saturated_for >= min_saturation_iters;
    status.algorithm_has_succeeded = status.algorithm_has_converged and status.algorithm_has_saturated;
    status.algorithm_has_got_stuck = status.algorithm_has_saturated and not status.algorithm_has_succeeded;
    if(tensors.state->position_is_any_edge()) status.algorithm_has_stuck_for = status.algorithm_has_got_stuck ? status.algorithm_has_stuck_for + 1 : 0;
    status.algorithm_has_to_stop = status.algorithm_has_stuck_for >= max_stuck_iters;

    if(tensors.state->position_is_any_edge()) {
        tools::log->debug("Simulation has converged: {}", status.algorithm_has_converged);
        tools::log->debug("Simulation has saturated: {}", status.algorithm_has_saturated);
        tools::log->debug("Simulation has succeeded: {}", status.algorithm_has_succeeded);
        tools::log->debug("Simulation has got stuck: {}", status.algorithm_has_got_stuck);
        tools::log->debug("Simulation has stuck for: {}", status.algorithm_has_stuck_for);
        tools::log->debug("Simulation has to stop  : {}", status.algorithm_has_to_stop);
    }

    tools::common::profile::t_con->toc();
}

void class_xdmrg::randomize_into_state_in_energy_window(ResetReason reason, StateType state_type, std::optional<std::string> sector) {
    tools::log->info("Resetting to state in energy window -- reason: {}", enum2str(reason));
    tools::log->info("Searching for state in normalized energy range: {} +- {}", status.energy_dens_target, status.energy_dens_window);

    status.num_resets++;
    if(reason == ResetReason::SATURATED and status.num_resets > settings::strategy::max_resets) {
        tools::log->info("Not allowed more resets due to saturation: num resets {} > max resets {}", status.num_resets, settings::strategy::max_resets);
        return;
    }

    int  counter           = 0;
    bool outside_of_window = true;
    tensors.activate_sites(settings::precision::max_size_full_diag, 2);
    while(true) {
        randomize_state(ResetReason::FIND_WINDOW, state_type, sector, -1); // Do not use the bitfield: set to -1
        status.energy_dens = tools::finite::measure::energy_normalized(tensors, status.energy_min_per_site, status.energy_max_per_site);
        outside_of_window  = std::abs(status.energy_dens - status.energy_dens_target) >= status.energy_dens_window;
        tools::log->info("New energy density: {:.16f} | window {} | outside of window: {}", status.energy_dens, status.energy_dens_window, outside_of_window);
        if(not outside_of_window) break;
        counter++;
        if(counter >= 200) throw std::runtime_error(fmt::format("Failed to find initial state in energy window after {}. retries: ", counter));
        if(counter % 10 == 0 and energy_window_growth_factor != 1.0) {
            double old_energy_dens_window = status.energy_dens_window;
            double new_energy_dens_window = std::min(energy_window_growth_factor * status.energy_dens_window, 0.5);

            tools::log->info("Can't find state in energy window.  Increasing energy window: {} --> {}", old_energy_dens_window, new_energy_dens_window);
            status.energy_dens_window = new_energy_dens_window;
        }
    }
    tools::log->info("Energy initial (per site) = {:.16f} | density = {:.8f} | retries = {}", tools::finite::measure::energy_per_site(tensors),
                     status.energy_dens, counter);
    clear_convergence_status();
    init_energy_limits(std::nullopt, status.energy_dens_window);
    tools::log->info("Number of product state resets: {}", status.num_resets);
}

void class_xdmrg::find_energy_range() {
    tools::log->trace("Finding energy range");
    // Here we define a set of tasks for fdmrg in order to produce the lowest and highest energy eigenstates,
    // We don't want it to randomize its own model, so we implant our current model before running the tasks.

    std::list<fdmrg_task> gs_tasks = {fdmrg_task::INIT_CLEAR_STATUS, fdmrg_task::INIT_BOND_DIM_LIMITS,
                                      fdmrg_task::INIT_WRITE_MODEL,  fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE,
                                      fdmrg_task::FIND_GROUND_STATE, fdmrg_task::POST_WRITE_RESULT};
    std::list<fdmrg_task> hs_tasks = {fdmrg_task::INIT_CLEAR_STATUS, fdmrg_task::INIT_BOND_DIM_LIMITS, fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE,
                                      fdmrg_task::FIND_HIGHEST_STATE, fdmrg_task::POST_WRITE_RESULT};
    class_fdmrg           fdmrg(h5pp_file);
    *fdmrg.tensors.model = *tensors.model; // Copy the model
    // Find loewst energy state
    tools::log = Logger::setLogger(std::string(enum2str(algo_type)) + "-gs", settings::console::verbosity, settings::console::timestamp);
    fdmrg.run_task_list(gs_tasks);
    status.energy_min_per_site = tools::finite::measure::energy_per_site(fdmrg.tensors);

    // Find highest energy state
    tools::log = Logger::setLogger(std::string(enum2str(algo_type)) + "-hs", settings::console::verbosity, settings::console::timestamp);
    fdmrg.run_task_list(hs_tasks);
    status.energy_max_per_site = tools::finite::measure::energy_per_site(fdmrg.tensors);

    tools::log = Logger::getLogger(std::string(enum2str(algo_type)));
}

bool   class_xdmrg::cfg_algorithm_is_on() { return settings::xdmrg::on; }
long   class_xdmrg::cfg_chi_lim_max() { return settings::xdmrg::chi_lim_max; }
size_t class_xdmrg::cfg_print_freq() { return settings::xdmrg::print_freq; }
bool   class_xdmrg::cfg_chi_lim_grow() { return settings::xdmrg::chi_lim_grow; }
long   class_xdmrg::cfg_chi_lim_init() { return settings::xdmrg::chi_lim_init; }
bool   class_xdmrg::cfg_store_wave_function() { return settings::xdmrg::store_wavefn; }
