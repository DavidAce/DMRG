//
// Created by david on 2018-02-09.
//

#include "class_xdmrg.h"
#include "class_fdmrg.h"
#include <config/nmspc_settings.h>
#include <general/nmspc_exceptions.h>
#include <general/nmspc_iter.h>
#include <math/linalg/tensor.h>
#include <math/num.h>
#include <math/rnd.h>
#include <physics/nmspc_quantum_mechanics.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/io.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/multisite.h>
#include <tools/finite/opt.h>
#include <tools/finite/opt_mps.h>

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
    if(state_prefix.empty()) throw except::state_error("no valid state candidates found for resume");
    tools::log->info("Resuming state [{}]", state_prefix);
    tools::finite::io::h5resume::load_simulation(*h5pp_file, state_prefix, tensors, status, algo_type);

    // Our first task is to decide on a state name for the newly loaded state
    // The simplest is to infer it from the state prefix itself
    auto name   = tools::common::io::h5resume::extract_state_name(state_prefix);
    auto number = tools::common::io::h5resume::extract_state_number(state_prefix);
    if(number) {
        excited_state_number = number.value();
        tensors.state->set_name(fmt::format("state_{}", excited_state_number));
    } else if(not name.empty())
        tensors.state->set_name(name);

    // Initialize a custom task list
    std::deque<xdmrg_task> task_list;

    if(status.algorithm_has_succeeded)
        task_list = {xdmrg_task::POST_PRINT_RESULT};
    else
        task_list = {xdmrg_task::INIT_CLEAR_CONVERGENCE, xdmrg_task::FIND_EXCITED_STATE,
                     xdmrg_task::POST_DEFAULT}; // Probably a savepoint. Simply "continue" the algorithm until convergence

    // If we reached this point the current state has finished for one reason or another.
    // We may still have some more things to do, e.g. the config may be asking for more states
    // Example 1:
    //      max_states = 4
    //      excited_state_number = 1
    //      excited_states = excited_state_number + 1 = 2  <-- Due to counting from 0
    //      missing_states = max_states - excited_states  = 2
    // Example 2:
    //      max_states = 2
    //      excited_state_number = 1
    //      excited_states = excited_state_number + 1 = 2  <-- Due to counting from 0
    //      missing_states = max_states - excited_states  = 0
    auto excited_states = excited_state_number + 1;
    auto missing_states = std::max(0ul, settings::xdmrg::max_states - excited_states);
    for(size_t new_state_num = 0; new_state_num < missing_states; new_state_num++) {
        task_list.emplace_back(xdmrg_task::PROF_RESET);
        switch(settings::strategy::secondary_states) {
            case StateInit::RANDOM_PRODUCT_STATE: task_list.emplace_back(xdmrg_task::NEXT_RANDOMIZE_INTO_STATE_IN_WIN); break;
            case StateInit::RANDOM_ENTANGLED_STATE: task_list.emplace_back(xdmrg_task::NEXT_RANDOMIZE_INTO_ENTANGLED_STATE); break;
            case StateInit::PRODUCT_STATE_ALIGNED: throw std::runtime_error("TODO! Product state aligned initialization not implemented yet"); break;
            case StateInit::PRODUCT_STATE_NEEL: throw std::runtime_error("TODO! Product state neel initialization not implemented yet"); break;
            case StateInit::RANDOMIZE_PREVIOUS_STATE: task_list.emplace_back(xdmrg_task::NEXT_RANDOMIZE_PREVIOUS_STATE); break;
        }
        task_list.emplace_back(xdmrg_task::FIND_EXCITED_STATE);
        task_list.emplace_back(xdmrg_task::POST_DEFAULT);
    }

    tools::log->info("Resuming task list:");
    for(const auto &task : task_list) tools::log->info(" -- {}", enum2str(task));

    run_task_list(task_list);
}

void class_xdmrg::run_default_task_list() {
    std::deque<xdmrg_task> default_task_list = {
        xdmrg_task::INIT_DEFAULT,
        xdmrg_task::FIND_EXCITED_STATE,
        xdmrg_task::POST_DEFAULT,
    };

    // Insert requested number of excited states
    for(size_t num = 1; num < settings::xdmrg::max_states; num++) {
        default_task_list.emplace_back(xdmrg_task::PROF_RESET);
        switch(settings::strategy::secondary_states) {
            case StateInit::RANDOM_PRODUCT_STATE: default_task_list.emplace_back(xdmrg_task::NEXT_RANDOMIZE_INTO_STATE_IN_WIN); break;
            case StateInit::RANDOM_ENTANGLED_STATE: default_task_list.emplace_back(xdmrg_task::NEXT_RANDOMIZE_INTO_ENTANGLED_STATE); break;
            case StateInit::RANDOMIZE_PREVIOUS_STATE: default_task_list.emplace_back(xdmrg_task::NEXT_RANDOMIZE_PREVIOUS_STATE); break;
            case StateInit::PRODUCT_STATE_ALIGNED: throw std::runtime_error("TODO! Product state aligned initialization not implemented yet"); break;
            case StateInit::PRODUCT_STATE_NEEL: throw std::runtime_error("TODO! Product state neel initialization not implemented yet"); break;
        }
        default_task_list.emplace_back(xdmrg_task::FIND_EXCITED_STATE);
        default_task_list.emplace_back(xdmrg_task::POST_DEFAULT);
    }
    run_task_list(default_task_list);
}

void class_xdmrg::run_task_list(std::deque<xdmrg_task> &task_list) {
    while(not task_list.empty()) {
        auto task = task_list.front();
        switch(task) {
            case xdmrg_task::INIT_RANDOMIZE_MODEL: randomize_model(); break;
            case xdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE: randomize_state(ResetReason::INIT, StateInit::RANDOM_PRODUCT_STATE); break;
            case xdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE: randomize_state(ResetReason::INIT, StateInit::RANDOM_ENTANGLED_STATE); break;
            case xdmrg_task::INIT_RANDOMIZE_FROM_CURRENT_STATE: randomize_state(ResetReason::INIT, StateInit::RANDOMIZE_PREVIOUS_STATE); break;
            case xdmrg_task::INIT_RANDOMIZE_INTO_STATE_IN_WIN:
                randomize_into_state_in_energy_window(ResetReason::INIT, settings::strategy::initial_state);
                break;
            case xdmrg_task::INIT_BOND_DIM_LIMITS: init_bond_dimension_limits(); break;
            case xdmrg_task::INIT_ENERGY_LIMITS: init_energy_limits(); break;
            case xdmrg_task::INIT_WRITE_MODEL: write_to_file(StorageReason::MODEL); break;
            case xdmrg_task::INIT_CLEAR_STATUS: status.clear(); break;
            case xdmrg_task::INIT_CLEAR_CONVERGENCE: clear_convergence_status(); break;
            case xdmrg_task::INIT_DEFAULT:
                run_preprocessing();
                break;
                //            case xdmrg_task::NEXT_TRUNCATE_ALL_SITES: truncate_all_sites(); break;
            case xdmrg_task::NEXT_RANDOMIZE_INTO_PRODUCT_STATE: randomize_state(ResetReason::NEW_STATE, StateInit::RANDOM_PRODUCT_STATE); break;
            case xdmrg_task::NEXT_RANDOMIZE_INTO_ENTANGLED_STATE: randomize_state(ResetReason::NEW_STATE, StateInit::RANDOM_ENTANGLED_STATE); break;
            case xdmrg_task::NEXT_RANDOMIZE_PREVIOUS_STATE: randomize_state(ResetReason::NEW_STATE, StateInit::RANDOMIZE_PREVIOUS_STATE); break;
            case xdmrg_task::NEXT_RANDOMIZE_INTO_STATE_IN_WIN: randomize_state(ResetReason::NEW_STATE, settings::strategy::initial_state); break;
            case xdmrg_task::FIND_ENERGY_RANGE: find_energy_range(); break;
            case xdmrg_task::FIND_EXCITED_STATE:
                tensors.state->set_name(fmt::format("state_{}", excited_state_number));
                run_algorithm();
                break;
            case xdmrg_task::POST_WRITE_RESULT: write_to_file(StorageReason::FINISHED); break;
            case xdmrg_task::POST_PRINT_RESULT: print_status_full(); break;
            case xdmrg_task::POST_PRINT_PROFILING: tools::common::profile::print_profiling(algo_type); break;
            case xdmrg_task::POST_DEFAULT: run_postprocessing(); break;
            case xdmrg_task::PROF_RESET: tools::common::profile::reset_profiling(algo_type); break;
        }
        task_list.pop_front();
    }
    if(not task_list.empty()) {
        for(auto &task : task_list) tools::log->critical("Unfinished task: {}", enum2str(task));
        throw std::runtime_error("Simulation ended with unfinished tasks");
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
    auto t_pre = tools::common::profile::prof[algo_type]["t_pre"]->tic_token();
    status.clear();
    randomize_model(); // First use of random!
    init_bond_dimension_limits();
    find_energy_range();
    init_energy_limits();
    if(settings::xdmrg::energy_density_window != 0.5)
        randomize_into_state_in_energy_window(ResetReason::INIT, settings::strategy::initial_state);
    else
        randomize_state(ResetReason::INIT, settings::strategy::initial_state);
    write_to_file(StorageReason::MODEL);
    tools::log->info("Finished {} preprocessing", algo_name);
}

void class_xdmrg::run_algorithm() {
    if(tensors.state->get_name().empty()) tensors.state->set_name(fmt::format("state_{}", excited_state_number));
    tools::log->info("Starting {} simulation of model [{}] for state [{}]", algo_name, enum2str(settings::model::model_type), tensors.state->get_name());
    auto t_sim  = tools::common::profile::prof[algo_type]["t_sim"]->tic_token();
    stop_reason = StopReason::NONE;
    while(true) {
        tools::log->trace("Starting step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);
        single_xDMRG_step();
        check_convergence();
        write_to_file();
        print_status_update();
        print_profiling_lap();

        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);

        // It's important not to perform the last move, so we break now: that last state would not get optimized
        if(stop_reason != StopReason::NONE) break;

        // Prepare for next step

        // Updating bond dimension must go first since it decides based on truncation error, but a projection+normalize resets truncation.
        update_bond_dimension_limit(); // Will update bond dimension if the state precision is being limited by bond dimension
        reduce_mpo_energy();
        try_discard_small_schmidt();
        try_bond_dimension_quench();
        try_disorder_damping();
        try_hamiltonian_perturbation();
        try_projection();
        try_full_expansion();
        move_center_point();
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", algo_name, tensors.state->get_name(), enum2str(stop_reason));
    status.algorithm_has_finished = true;
}

std::vector<class_xdmrg::OptConf> class_xdmrg::get_opt_conf_list() {
    tools::log->trace("Configuring xDMRG optimization trial");
    std::vector<OptConf> configs;

    /*
     *
     *  First trial
     *
     */
    OptConf c1;
    c1.label = "c1";

    // The first decision is easy. Real or complex optimization
    if(tensors.is_real()) c1.optType = OptType::REAL;

    // Normally we do 2-site dmrg, unless settings specifically ask for 1-site
    c1.max_sites = std::min(2ul, settings::strategy::multisite_mps_size_def);

    // Next we setup the mode at the early stages of the simulation
    // Note that we make stricter requirements as we go down the if-list

    if(status.iter <= 1) {
        // On the zeroth iteration, the rest of the chain is probably just randomized, so don't focus on growing the bond dimension yet.
        status.chi_lim = std::max<long>(tensors.state->find_largest_chi(), cfg_chi_lim_init());
    }

    if(status.iter < settings::xdmrg::olap_iters + settings::xdmrg::vsub_iters) {
        // If early in the simulation, and the bond dimension is small enough we use subspace optimization
        c1.optMode       = OptMode::VARIANCE;
        c1.optSpace      = OptSpace::SUBSPACE_AND_DIRECT;
        c1.max_sites     = settings::strategy::multisite_mps_size_init;
        c1.second_chance = false;
        if(settings::xdmrg::chi_lim_vsub > 0) status.chi_lim = settings::xdmrg::chi_lim_vsub;
    }

    if(status.iter < settings::xdmrg::olap_iters) {
        // Very early in the simulation it is worth just following the overlap to get the overall structure of the final state
        c1.optMode       = OptMode::OVERLAP;
        c1.optSpace      = OptSpace::SUBSPACE_ONLY;
        c1.max_sites     = settings::strategy::multisite_mps_size_init;
        c1.second_chance = false;
        if(settings::xdmrg::chi_lim_olap > 0) status.chi_lim = settings::xdmrg::chi_lim_olap;
    }
    if(num_discards > 0 and status.iter < iter_discard + 2) {
        // When discarded the bond dimensions are highly truncated. We can then afford subspace optimization
        c1.optMode       = OptMode::VARIANCE;
        c1.optSpace      = OptSpace::SUBSPACE_AND_DIRECT;
        c1.max_sites     = settings::strategy::multisite_mps_size_max;
        c1.second_chance = false;
        //        if(settings::xdmrg::chi_lim_vsub > 0) status.chi_lim = settings::xdmrg::chi_lim_vsub;
    }

    // Setup strong overrides to normal conditions, e.g.,
    //      - for experiments like perturbation or chi quench
    //      - when the algorithm has already converged

    if(status.variance_mpo_converged_for > 0) {
        // No need to do expensive operations -- just finish
        c1.optMode       = OptMode::VARIANCE;
        c1.optSpace      = OptSpace::DIRECT;
        c1.second_chance = false;
    }
    if(tensors.measurements.energy_variance.has_value() and
       tensors.measurements.energy_variance.value()< settings::precision::variance_convergence_threshold){
        // No need to do expensive operations -- just finish
        c1.optMode       = OptMode::VARIANCE;
        c1.optSpace      = OptSpace::DIRECT;
        c1.second_chance = false;
    }


    if(tensors.state->size_1site() > settings::precision::max_size_part_diag) {
        // Make sure not to use SUBSPACE if the 1-site problem size is huge
        // When this happens, we can't use SUBSPACE even with two sites,
        // so we might as well give it to DIRECT, which handles larger problems.
        c1.optMode = OptMode::VARIANCE;
        if(c1.optSpace == OptSpace::SUBSPACE_AND_DIRECT) c1.optSpace = OptSpace::DIRECT;
        if(c1.optSpace == OptSpace::SUBSPACE_ONLY) c1.optSpace = OptSpace::DIRECT;
    }

    // If we are doing 1-site dmrg, then we better use subspace expansion
    if(c1.max_sites == 1 or settings::strategy::multisite_mps_size_def == 1)
        c1.alpha_expansion = std::min(0.1, status.energy_variance_lowest); // Usually a good value to start with

    if(settings::strategy::expand_subspace_when_stuck and (status.algorithm_has_stuck_for > 0 or sub_expansion_alpha.has_value()) and
       num_expansion_iters < max_expansion_iters) {
        if(sub_expansion_alpha) c1.alpha_expansion = sub_expansion_alpha;
        if(not c1.alpha_expansion) c1.alpha_expansion = std::min(0.1, status.energy_variance_lowest);
        // Update alpha
        auto report = check_saturation(var_mpo_step, settings::precision::variance_saturation_sensitivity / 10);
        tools::log->info("Determining alpha: report computed {} | saturated {} | expansion iters {} ({})", report.has_computed, report.has_saturated,
                         num_expansion_iters, max_expansion_iters);
        if(report.has_computed) {
            double factor_up = std::pow(5.0, 1.0 / static_cast<double>(tensors.get_length()));
            double factor_dn = std::pow(0.1, 1.0 / static_cast<double>(tensors.get_length()));
            if(report.has_saturated) {
                c1.alpha_expansion = c1.alpha_expansion.value() * factor_up;
                c1.alpha_expansion = std::min(max_expansion_alpha, c1.alpha_expansion.value());
            } else {
                auto alpha_dn = c1.alpha_expansion.value() * factor_dn;
                if(alpha_dn > status.energy_variance_lowest)
                    c1.alpha_expansion = alpha_dn;
                else
                    c1.alpha_expansion = std::nullopt; // Back to normal
            }
        } else {
            tools::log->info("Could not enable alpha: Report wasn't computed");
        }
        sub_expansion_alpha = c1.alpha_expansion;
        if(tensors.position_is_inward_edge()) num_expansion_iters++;
    } else if(num_expansion_iters >= max_expansion_iters) {
        tools::log->info("Could not enable alpha: More than num_expansion iters ({}) >=  max_expansion_iters ({})", num_expansion_iters, max_expansion_iters);
        c1.alpha_expansion = std::nullopt; // Back to normal
    }
    sub_expansion_alpha = c1.alpha_expansion;

    // Setup the maximum problem size here
    switch(c1.optSpace) {
        case OptSpace::KRYLOV_ENERGY:
        case OptSpace::KRYLOV_VARIANCE:
        case OptSpace::DIRECT: c1.max_problem_size = settings::precision::max_size_direct; break;
        case OptSpace::SUBSPACE_ONLY:
        case OptSpace::SUBSPACE_AND_DIRECT: c1.max_problem_size = settings::precision::max_size_part_diag; break;
    }

    // We can make trials with different number of sites.
    // Eg if the simulation is stuck we may try with more sites.
    if(status.algorithm_has_stuck_for > max_stuck_iters / 2) c1.max_sites = settings::strategy::multisite_mps_size_max;
    if(status.algorithm_has_succeeded) c1.max_sites = c1.min_sites; // No need to do expensive operations -- just finish

    c1.chosen_sites = tools::finite::multisite::generate_site_list(*tensors.state, c1.max_problem_size, c1.max_sites, c1.min_sites);
    c1.problem_dims = tools::finite::multisite::get_dimensions(*tensors.state, c1.chosen_sites);
    c1.problem_size = tools::finite::multisite::get_problem_size(*tensors.state, c1.chosen_sites);

    configs.emplace_back(c1);
    if(not c1.second_chance) return configs;
    /*
     *
     *  Second trial
     *
     */
    OptConf c2 = c1; // Start with c1 as a baseline
    c2.label   = "c2";
    // NOTES
    // 1) OVERLAP does not get a second chance: It is supposed to pick best overlap, not improve variance
    // 2)
    //
    if(c2.optSpace == OptSpace::SUBSPACE_ONLY and c2.optMode == OptMode::VARIANCE and c1.max_sites < settings::strategy::multisite_mps_size_max) {
        // I.e. if we did a SUBSPACE run that did not result in better variance, try a few more sites
        c2.max_sites    = settings::strategy::multisite_mps_size_max;
        c2.chosen_sites = tools::finite::multisite::generate_site_list(*tensors.state, c2.max_problem_size, c2.max_sites, c2.min_sites);
        c2.problem_dims = tools::finite::multisite::get_dimensions(*tensors.state, c2.chosen_sites);
        c2.problem_size = tools::finite::multisite::get_problem_size(*tensors.state, c2.chosen_sites);
        if(c2.chosen_sites != c1.chosen_sites) configs.emplace_back(c2);
    } else if(c2.optSpace != OptSpace::DIRECT and c2.optMode == OptMode::VARIANCE) {
        // I.e. if we did a SUBSPACE VARIANCE run that did not improve, and we are already trying max sites, switch to DIRECT
        c2.optInit  = c1.alpha_expansion ? OptInit::CURRENT_STATE : OptInit::LAST_RESULT; // Use the result from c1
        c2.optSpace = OptSpace::DIRECT;
        c2.optMode  = OptMode::VARIANCE;
        configs.emplace_back(c2);
    }
    else if(c2.optSpace == OptSpace::DIRECT and status.algorithm_has_stuck_for > 0 and settings::strategy::krylov_opt_when_stuck) {
        c2.optSpace      = OptSpace::KRYLOV_VARIANCE;
        c2.optMode         = OptMode::VARIANCE;
        c2.max_sites       = std::min(2ul, settings::strategy::multisite_mps_size_def);
        c2.second_chance   = false;
        c2.alpha_expansion = std::nullopt;
        // Use the same number of sites if stuck and not more than two, so that we can reuse the lbfgs progress as a good starting point for arpack.
        if(status.algorithm_has_stuck_for >= 2) c2.max_sites = std::min(2ul, c1.max_sites);
//
        if (c1.max_sites == c2.max_sites and not c1.alpha_expansion) c2.optInit = OptInit::LAST_RESULT;
        else  c2.optInit = OptInit::CURRENT_STATE;

        c2.chosen_sites = tools::finite::multisite::generate_site_list(*tensors.state, c2.max_problem_size, c2.max_sites, c2.min_sites);
        c2.problem_dims = tools::finite::multisite::get_dimensions(*tensors.state, c2.chosen_sites);
        c2.problem_size = tools::finite::multisite::get_problem_size(*tensors.state, c2.chosen_sites);
        configs.emplace_back(c2);
    }


    if(not c2.second_chance) return configs;


    return configs;
}

void class_xdmrg::single_xDMRG_step(std::vector<class_xdmrg::OptConf> optConf) {
    using namespace tools::finite;
    using namespace tools::finite::opt;

    if(optConf.empty()) optConf = get_opt_conf_list();
    std::vector<opt_mps> results;
    std::optional<std::vector<class_mps_site>> mps_original = std::nullopt;
    variance_before_step = std::nullopt;
    tools::log->debug("Starting xDMRG iter {} | step {} | pos {} | dir {} | confs {}", status.iter, status.step, status.position, status.direction, optConf.size());
    for(const auto &[idx_conf, conf] : iter::enumerate(optConf)) {
        if(conf.optMode == OptMode::OVERLAP and conf.optSpace == OptSpace::DIRECT) throw std::logic_error("[OVERLAP] mode and [DIRECT] space are incompatible");
        tensors.activate_sites(conf.chosen_sites);
        if(tensors.active_sites.empty()) continue;
        if(not variance_before_step) variance_before_step = measure::energy_variance(tensors); // Should just take value from cache

        // Use subspace expansion if alpha_expansion is set
        // Note that this changes the mps and edges adjacent to "tensors.active_sites"
        if(conf.alpha_expansion) {
            auto pos_expanded = tensors.expand_subspace(std::nullopt, status.chi_lim); // nullopt implies a pos query
            if(not mps_original) mps_original = tensors.state->get_mps_sites(pos_expanded);
            tensors.expand_subspace(conf.alpha_expansion, status.chi_lim);
        }

        tools::log->debug("Running conf {}/{}: {} | init {} | mode {} | space {} | type {} | sites {} | dims {} = {} | alpha {:.3e}", idx_conf + 1,
                          optConf.size(), conf.label, enum2str(conf.optInit), enum2str(conf.optMode), enum2str(conf.optSpace), enum2str(conf.optType),
                          conf.chosen_sites, conf.problem_dims, conf.problem_size,
                          (conf.alpha_expansion ? conf.alpha_expansion.value() : std::numeric_limits<double>::quiet_NaN()));
        switch(conf.optInit) {
            case OptInit::CURRENT_STATE: {
                results.emplace_back(opt::find_excited_state(tensors, status, conf.optMode, conf.optSpace, conf.optType));
                results.back().validate_result();
                break;
            }
            case OptInit::LAST_RESULT: {
                if(results.empty()) throw std::logic_error("There are no previous results to select an initial state");
                results.emplace_back(opt::find_excited_state(tensors, results.back(), status, conf.optMode, conf.optSpace, conf.optType));
                results.back().validate_result();
                break;
            }
        }

        if(conf.alpha_expansion){
            results.back().set_alpha(conf.alpha_expansion);
            auto pos_expanded   = tensors.expand_subspace(std::nullopt, status.chi_lim); // nullopt implies a pos query
            results.back().mps_backup = tensors.state->get_mps_sites(pos_expanded); // Backup the mps sites that this run was compatible with
        }
        // Reset the mps to the original if they were backed up earlier
        // This is so that the next conf starts with unchanged mps as neighbors
        if (mps_original) {
            tensors.state->set_mps_sites(mps_original.value());
            tensors.rebuild_edges(); // Make sure the edges are up to date
        }



        // We can now decide if we are happy with the result or not.
        results.back().set_relchange(results.back().get_variance() / variance_before_step.value());

        /* clang-format off */
        if(results.back().get_relchange() > 1 )     conf.optExit |= OptExit::FAIL_WORSENED;
        if(results.back().get_relchange() > 0.999 ) conf.optExit |= OptExit::FAIL_NOCHANGE;
        if(results.back().get_grad_norm() > 0.100 ) conf.optExit |= OptExit::FAIL_GRADIENT;
        results.back().set_optexit(conf.optExit);

        tools::log->debug("Optimization [{}|{}]: {}. Variance change {:.4f} --> {:.4f} ({:.3f} %)",
                          enum2str(conf.optMode), enum2str(conf.optSpace), flag2str(conf.optExit),
                          std::log10(variance_before_step.value()),std::log10(results.back().get_variance()), results.back().get_relchange()*100);

        if(conf.optExit == OptExit::SUCCESS) break;
        if(not has_flag(conf.optExit, OptExit::FAIL_NOCHANGE)) break; // Good enough
        /* clang-format on */
    }

    if(not results.empty()) {
        if(tools::log->level() <= spdlog::level::debug)
            for(auto &candidate : results)
                tools::log->debug(
                    "Candidate: {} | sites [{:>2}-{:<2}] | alpha {:8.2e} ({} iters) | variance {:<12.6f} | energy {:<12.6f} | overlap {:.16f} | norm {:.16f} "
                    "| [{}][{}] | time {:.4f} ms",
                    candidate.get_name(), candidate.get_sites().front(), candidate.get_sites().back(), candidate.get_alpha(), num_expansion_iters,
                    std::log10(candidate.get_variance()), candidate.get_energy_per_site(), candidate.get_overlap(), candidate.get_norm(),
                    enum2str(candidate.get_optmode()), enum2str(candidate.get_optspace()), 1000 * candidate.get_time());

        // Sort the results in order of increasing variance
        auto comp_variance = [](const opt_mps &lhs, const opt_mps &rhs) {
            return lhs.get_variance() < rhs.get_variance() or lhs.get_optspace() == OptSpace::KRYLOV_ENERGY or lhs.get_optspace() == OptSpace::KRYLOV_VARIANCE;
        };
        std::sort(results.begin(), results.end(), comp_variance);

        // Take the best result
        const auto &winner = results.front();
        last_optspace      = winner.get_optspace();
        last_optmode       = winner.get_optmode();
        tensors.activate_sites(winner.get_sites());
        tensors.state->tag_active_sites_normalized(false);
        if(not winner.mps_backup.empty()) {
            // Restore the neighboring mps that are compatible with this winner
            tensors.state->set_mps_sites(winner.mps_backup);
            tensors.rebuild_edges();
        }

    // Do the truncation with SVD
        tensors.merge_multisite_tensor(winner.get_tensor(), status.chi_lim);
        if(tools::log->level() <= spdlog::level::debug) {
            auto truncation_errors = tensors.state->get_truncation_errors_active();
            for(auto &t : truncation_errors) t = std::log10(t);
            tools::log->debug("Truncation errors: {:.4f}", fmt::join(truncation_errors, ", "));
        }

        if constexpr(settings::debug) {
            auto variance_before_svd = winner.get_variance();
            auto variance_after_svd  = tools::finite::measure::energy_variance(tensors);
            tools::log->info("Variance check before SVD: {:.16f}", std::log10(variance_before_svd));
            tools::log->info("Variance check after  SVD: {:.16f}", std::log10(variance_after_svd));
            tools::log->info("Variance delta due to SVD: {:.16f}", 100 * variance_after_svd / variance_before_svd);
        }

        // Update current energy density ε
        status.energy_dens =
            (tools::finite::measure::energy_per_site(tensors) - status.energy_min_per_site) / (status.energy_max_per_site - status.energy_min_per_site);

        if(not tensors.active_sites.empty()) {
            tools::log->trace("Updating variance record holder");
            auto var = tools::finite::measure::energy_variance(tensors);
            if(var < status.energy_variance_lowest) status.energy_variance_lowest = var;
            var_mpo_step.emplace_back(var);
        }
    }

    status.wall_time = tools::common::profile::t_tot->get_measured_time();
    status.algo_time = tools::common::profile::prof[algo_type]["t_sim"]->get_measured_time();
}

void class_xdmrg::check_convergence() {
    if(not tensors.position_is_inward_edge()) return;
    auto t_con = tools::common::profile::prof[algo_type]["t_con"]->tic_token();

    // TODO: Move this reset block away from here
    //    bool outside_of_window = std::abs(status.energy_dens - status.energy_dens_target) > status.energy_dens_window;
    //    if(status.iter > 2 and tensors.position_is_inward_edge()) {
    //        if(outside_of_window and
    //           (status.variance_mpo_has_saturated or status.variance_mpo_has_converged or tools::finite::measure::energy_variance_per_site(tensors) < 1e-4)) {
    //            double      old_energy_dens_window = status.energy_dens_window;
    //            double      new_energy_dens_window = std::min(energy_window_growth_factor * status.energy_dens_window, 0.5);
    //            std::string reason = fmt::format("energy {:.16f} saturated outside of energy window {} ± {}",
    //            tools::finite::measure::energy_per_site(tensors),
    //                                             status.energy_dens_target, status.energy_dens_window);
    //            tools::log->info("Increasing energy window: {} --> {}", old_energy_dens_window, new_energy_dens_window);
    //            status.energy_dens_window = new_energy_dens_window;
    //            randomize_into_state_in_energy_window(ResetReason::SATURATED, settings::strategy::initial_state, settings::strategy::target_sector);
    //        }
    //    }
    update_variance_max_digits();
    check_convergence_variance();
    check_convergence_entg_entropy();
    check_convergence_spin_parity_sector(settings::strategy::target_sector);
//    std::max(status.variance_mpo_saturated_for, status.entanglement_saturated_for) > max_saturation_iters or
    if(status.variance_mpo_saturated_for > 0 and status.entanglement_saturated_for > 0)
        status.algorithm_saturated_for++;
    else
        status.algorithm_saturated_for = 0;

    if(status.variance_mpo_converged_for > 0 and status.entanglement_converged_for > 0 and status.spin_parity_has_converged)
        status.algorithm_converged_for++;
    else
        status.algorithm_converged_for = 0;

    if(status.algorithm_saturated_for > 0 and status.algorithm_converged_for == 0)
        status.algorithm_has_stuck_for++;
    else
        status.algorithm_has_stuck_for = 0;

    if(status.algorithm_converged_for == 0 and status.variance_mpo_saturated_for * status.entanglement_saturated_for == 0 and status.algorithm_has_stuck_for != 0)
        throw std::logic_error("Should have zeroed");
    if(status.algorithm_converged_for == 0 and status.variance_mpo_saturated_for * status.entanglement_saturated_for > 0 and status.algorithm_has_stuck_for == 0)
        throw std::logic_error("Should not have zeroed");

    status.algorithm_has_succeeded = status.algorithm_converged_for >= min_converged_iters and status.algorithm_saturated_for >= min_saturation_iters;
    status.algorithm_has_to_stop   = status.algorithm_has_stuck_for >= max_stuck_iters;

    tools::log->info("Algorithm report: converged {} | saturated {} | stuck {} | succeeded {} | has to stop {} | var prec limit {:.6f}",
                     status.algorithm_converged_for, status.algorithm_saturated_for, status.algorithm_has_stuck_for, status.algorithm_has_succeeded,
                     status.algorithm_has_to_stop, std::log10(status.energy_variance_prec_limit));

    stop_reason = StopReason::NONE;
    if(status.iter >= settings::xdmrg::min_iters and not tensors.model->is_perturbed() and not tensors.model->is_damped()) {
        if(status.iter >= settings::xdmrg::max_iters) stop_reason = StopReason::MAX_ITERS;
        if(status.algorithm_has_succeeded) stop_reason = StopReason::SUCCEEDED;
        if(status.algorithm_has_to_stop) stop_reason = StopReason::SATURATED;
        if(status.num_resets > settings::strategy::max_resets) stop_reason = StopReason::MAX_RESET;
        if(status.entanglement_saturated_for > 0 and settings::xdmrg::finish_if_entanglm_saturated) stop_reason = StopReason::SATURATED;
        if(status.variance_mpo_saturated_for > 0 and settings::xdmrg::finish_if_variance_saturated) stop_reason = StopReason::SATURATED;
        if(settings::strategy::randomize_early and excited_state_number == 0 and tensors.state->find_largest_chi() >= 32 and
           tools::finite::measure::energy_variance(tensors) < 1e-4)
            stop_reason = StopReason::RANDOMIZE;
    }
}

void class_xdmrg::randomize_into_state_in_energy_window(ResetReason reason, StateInit state_type, std::optional<std::string> sector) {
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
        randomize_state(ResetReason::FIND_WINDOW, state_type, std::nullopt, sector, -1); // Do not use the bitfield: set to -1
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

    std::deque<fdmrg_task> gs_tasks = {
        fdmrg_task::INIT_CLEAR_STATUS, fdmrg_task::INIT_BOND_DIM_LIMITS, fdmrg_task::INIT_WRITE_MODEL,    fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE,
        fdmrg_task::FIND_GROUND_STATE, fdmrg_task::POST_WRITE_RESULT,    fdmrg_task::POST_PRINT_PROFILING};
    std::deque<fdmrg_task> hs_tasks = {fdmrg_task::INIT_CLEAR_STATUS,  fdmrg_task::INIT_BOND_DIM_LIMITS, fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE,
                                       fdmrg_task::FIND_HIGHEST_STATE, fdmrg_task::POST_WRITE_RESULT,    fdmrg_task::POST_PRINT_PROFILING};
    // Find lowest energy state
    {
        class_fdmrg fdmrg_gs(h5pp_file);
        *fdmrg_gs.tensors.model = *tensors.model; // Copy the model
        tools::log = tools::Logger::setLogger(std::string(enum2str(algo_type)) + "-gs", settings::console::verbosity, settings::console::timestamp);
        fdmrg_gs.run_task_list(gs_tasks);
        status.energy_min_per_site = tools::finite::measure::energy_per_site(fdmrg_gs.tensors);
    }

    // Find highest energy state
    {
        class_fdmrg fdmrg_hs(h5pp_file);
        *fdmrg_hs.tensors.model = *tensors.model; // Copy the model
        tools::log = tools::Logger::setLogger(std::string(enum2str(algo_type)) + "-hs", settings::console::verbosity, settings::console::timestamp);
        fdmrg_hs.run_task_list(hs_tasks);
        status.energy_max_per_site = tools::finite::measure::energy_per_site(fdmrg_hs.tensors);
    }

    // Reset our logger
    tools::log = tools::Logger::getLogger(std::string(enum2str(algo_type)));

    // Reset the profiling data
    tools::common::profile::reset_profiling(AlgorithmType::fDMRG);

    // Set the default profile back to xDMRG because the constructor of class_fdmrg changed it
    tools::common::profile::set_default_prof(AlgorithmType::xDMRG);
}

void class_xdmrg::update_time_step() {
    tools::log->trace("Updating time step");
    status.delta_t = std::complex<double>(1e-6, 0);
}

void class_xdmrg::create_hamiltonian_gates() {
    tools::log->info("Creating Hamiltonian gates");
    ham_gates_1body.clear();
    ham_gates_2body.clear();
    ham_gates_3body.clear();
    // Create the hamiltonian gates with n-site terms
    auto list_1site = num::range<size_t>(0, settings::model::model_size - 0, 1);
    auto list_2site = num::range<size_t>(0, settings::model::model_size - 1, 1);
    auto list_3site = num::range<size_t>(0, settings::model::model_size - 2, 1);
    //    for(auto pos : list_1site) ham_gates_1body.emplace_back(qm::Gate(tensors.model->get_multisite_ham({pos}, {1}), {pos},
    //    tensors.state->get_spin_dims({pos}))); for(auto pos : list_2site) ham_gates_2body.emplace_back(qm::Gate(tensors.model->get_multisite_ham({pos, pos +
    //    1}, {2}), {pos, pos + 1},tensors.state->get_spin_dims({pos, pos + 1}))); for(auto pos : list_3site)
    //    ham_gates_3body.emplace_back(qm::Gate(tensors.model->get_multisite_ham({pos, pos + 1, pos + 2}, {3}), {pos, pos + 1, pos + 2},
    //    tensors.state->get_spin_dims({pos, pos + 1, pos + 2})));
    //
    for(auto pos : list_1site) ham_gates_1body.emplace_back(qm::Gate(tensors.model->get_multisite_ham({pos}, {1}), {pos}, tensors.state->get_spin_dims({pos})));
    for(auto pos : list_2site)
        ham_gates_2body.emplace_back(
            qm::Gate(tensors.model->get_multisite_ham({pos, pos + 1}, {2}), {pos, pos + 1}, tensors.state->get_spin_dims({pos, pos + 1})));
    for(auto pos : list_3site)
        ham_gates_3body.emplace_back(qm::Gate(tensors.model->get_multisite_ham({pos, pos + 1, pos + 2}, {3}), {pos, pos + 1, pos + 2},
                                              tensors.state->get_spin_dims({pos, pos + 1, pos + 2})));
}

void class_xdmrg::create_time_evolution_gates() {
    // Create the time evolution operators
    if(std::abs(status.delta_t) == 0) update_time_step();
    tools::log->info("Creating time evolution gates");
    time_gates_1site = qm::timeEvolution::get_time_evolution_gates(status.delta_t, ham_gates_1body);
    time_gates_2site = qm::timeEvolution::get_time_evolution_gates(status.delta_t, ham_gates_2body);
    time_gates_3site = qm::timeEvolution::get_time_evolution_gates(status.delta_t, ham_gates_3body);
}

bool   class_xdmrg::cfg_algorithm_is_on() { return settings::xdmrg::on; }
long   class_xdmrg::cfg_chi_lim_max() { return settings::xdmrg::chi_lim_max; }
size_t class_xdmrg::cfg_print_freq() { return settings::xdmrg::print_freq; }
bool   class_xdmrg::cfg_chi_lim_grow() { return settings::xdmrg::chi_lim_grow; }
long   class_xdmrg::cfg_chi_lim_init() { return settings::xdmrg::chi_lim_init; }
bool   class_xdmrg::cfg_store_wave_function() { return settings::xdmrg::store_wavefn; }
