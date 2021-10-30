#include "xdmrg.h"
#include "fdmrg.h"
#include <config/settings.h>
#include <debug/exceptions.h>
#include <general/iter.h>
#include <io/fmt.h>
#include <math/linalg/tensor.h>
#include <math/num.h>
#include <math/rnd.h>
#include <qm/time.h>
#include <tensors/edges/EdgesFinite.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/state/StateFinite.h>
#include <tid/tid.h>
#include <tools/common/h5.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/h5.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/multisite.h>
#include <tools/finite/opt.h>
#include <tools/finite/opt_meta.h>
#include <tools/finite/opt_mps.h>
#include <tools/finite/print.h>

xdmrg::xdmrg(std::shared_ptr<h5pp::File> h5ppFile_) : AlgorithmFinite(std::move(h5ppFile_), AlgorithmType::xDMRG) {
    tools::log->trace("Constructing class_xdmrg");
}

void xdmrg::resume() {
    // Resume can imply many things
    // 1) Resume a simulation which terminated prematurely
    // 2) Resume a previously successful simulation. This may be desireable if the config
    //    wants something that is not present in the file.
    //      a) A certain number of states
    //      b) A state inside of a particular energy window
    //      c) The ground or "roof" states
    // To guide the behavior, we check the setting ResumePolicy.

    auto resumable_states = tools::common::h5::resume::find_resumable_states(*h5file, status.algo_type);
    if(resumable_states.empty()) throw except::state_error("no resumable states were found");

    for(const auto &state_prefix : resumable_states) {
        tools::log->info("Resuming state [{}]", state_prefix);
        tools::finite::h5::load::simulation(*h5file, state_prefix, tensors, status, status.algo_type);

        // Our first task is to decide on a state name for the newly loaded state
        // The simplest is to infer it from the state prefix itself
        auto name   = tools::common::h5::resume::extract_state_name(state_prefix);
        auto number = tools::common::h5::resume::extract_state_number(state_prefix);
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
        run_task_list(task_list);
    }

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

    // Initialize another custom task list
    std::deque<xdmrg_task> task_list;
    auto                   excited_states = excited_state_number + 1;
    auto                   missing_states = std::max(0ul, settings::xdmrg::max_states - excited_states);
    for(size_t new_state_num = 0; new_state_num < missing_states; new_state_num++) {
        task_list.emplace_back(xdmrg_task::TIMER_RESET);
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
    for(const auto &task : task_list) tools::log->info(" -- {}", enum2sv(task));

    run_task_list(task_list);
}

void xdmrg::run_default_task_list() {
    std::deque<xdmrg_task> default_task_list = {
        xdmrg_task::INIT_DEFAULT,
        xdmrg_task::FIND_EXCITED_STATE,
        xdmrg_task::POST_DEFAULT,
    };

    // Insert requested number of excited states
    for(size_t num = 1; num < settings::xdmrg::max_states; num++) {
        default_task_list.emplace_back(xdmrg_task::TIMER_RESET);
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

void xdmrg::run_task_list(std::deque<xdmrg_task> &task_list) {
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
            case xdmrg_task::POST_PRINT_TIMERS: tools::common::timer::print_timers(); break;
            case xdmrg_task::POST_DEFAULT: run_postprocessing(); break;
            case xdmrg_task::TIMER_RESET: tid::reset("xDMRG"); break;
        }
        task_list.pop_front();
    }
    if(not task_list.empty()) {
        for(auto &task : task_list) tools::log->critical("Unfinished task: {}", enum2sv(task));
        throw std::runtime_error("Simulation ended with unfinished tasks");
    }
}

void xdmrg::init_energy_limits(std::optional<double> energy_density_target, std::optional<double> energy_density_window) {
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

void xdmrg::run_preprocessing() {
    tools::log->info("Running {} preprocessing", status.algo_type_sv());
    auto t_pre = tid::tic_scope("pre");
    status.clear();
    randomize_model(); // First use of random!
    tools::finite::print::model(*tensors.model);
    init_bond_dimension_limits();
    find_energy_range();
    init_energy_limits();
    if(settings::xdmrg::energy_density_window != 0.5)
        randomize_into_state_in_energy_window(ResetReason::INIT, settings::strategy::initial_state);
    else
        randomize_state(ResetReason::INIT, settings::strategy::initial_state);
    write_to_file(StorageReason::MODEL);
    tools::log->info("Finished {} preprocessing", status.algo_type_sv());
}

void xdmrg::run_algorithm() {
    if(tensors.state->get_name().empty()) tensors.state->set_name(fmt::format("state_{}", excited_state_number));
    tools::log->info("Starting {} simulation of model [{}] for state [{}]", status.algo_type_sv(), enum2sv(settings::model::model_type),
                     tensors.state->get_name());
    auto t_run       = tid::tic_scope("run");
    status.algo_stop = AlgorithmStop::NONE;

    while(true) {
        tools::log->trace("Starting step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);
        single_xDMRG_step();
        print_status_update();
        check_convergence();
        write_to_file();

        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);

        // It's important not to perform the last move, so we break now: that last state would not get optimized
        if(status.algo_stop != AlgorithmStop::NONE) break;

        // Prepare for next step

        // Updating bond dimension must go first since it decides based on truncation error, but a projection+normalize resets truncation.
        update_bond_dimension_limit();   // Will update bond dimension if the state precision is being limited by bond dimension
        update_expansion_factor_alpha(); // Will update the subspace expansion factor
        reduce_mpo_energy();
        try_bond_dimension_quench();
        try_hamiltonian_perturbation();
        try_projection();
        try_full_expansion();
        move_center_point();
        status.wall_time = tid::get_unscoped("t_tot").get_time();
        status.algo_time = t_run->get_time();
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", status.algo_type_sv(), tensors.state->get_name(), status.algo_stop_sv());
    status.algorithm_has_finished = true;
}

std::vector<xdmrg::OptConf> xdmrg::get_opt_conf_list() {
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
    if(c1.optType == OptType::CPLX) throw std::runtime_error("");
    // Normally we do 2-site dmrg, unless settings specifically ask for 1-site
    c1.max_sites = std::min(2ul, settings::strategy::multisite_mps_size_def);

    // Set the default gradient precision
    c1.max_grad_tolerance = settings::precision::max_grad_tolerance;

    // When stuck we can require better accuracy on the gradient
    if(status.algorithm_has_stuck_for > 0) c1.max_grad_tolerance = 1e-2 * settings::precision::max_grad_tolerance;

    // Next we setup the mode at the early stages of the simulation
    // Note that we make stricter requirements as we go down the if-list

    if(status.iter <= 1) {
        // On the zeroth iteration, the rest of the chain is probably just randomized, so don't focus on growing the bond dimension yet.
        status.chi_lim = std::max<long>(tensors.state->find_largest_chi(), settings::chi_lim_init(status.algo_type));
    }

    if(status.iter < settings::xdmrg::olap_iters + settings::xdmrg::vsub_iters) {
        // If early in the simulation, and the bond dimension is small enough we use subspace optimization
        c1.optMode   = OptMode::VARIANCE;
        c1.optSpace  = OptSpace::SUBSPACE;
        c1.max_sites = settings::strategy::multisite_mps_size_max;
        c1.retry     = true;
        if(settings::xdmrg::chi_lim_vsub > 0) status.chi_lim = settings::xdmrg::chi_lim_vsub;
    }

    if(status.iter < settings::xdmrg::olap_iters) {
        // Very early in the simulation it is worth just following the overlap to get the overall structure of the final state
        c1.optMode   = OptMode::OVERLAP;
        c1.optSpace  = OptSpace::SUBSPACE;
        c1.max_sites = settings::strategy::multisite_mps_size_max;
        c1.retry     = false;
        if(settings::xdmrg::chi_lim_olap > 0) status.chi_lim = settings::xdmrg::chi_lim_olap;
    }

    if(settings::strategy::krylov_opt_when_stuck and status.algorithm_has_stuck_for > 0) {
        c1.optMode   = OptMode::VARIANCE;
        c1.optSpace  = OptSpace::KRYLOV;
        c1.max_sites = settings::strategy::multisite_mps_size_def;
        c1.retry     = true;
    }

    // Setup strong overrides to normal conditions, e.g.,
    //      - for experiments like perturbation or chi quench
    //      - when the algorithm has already converged
    //    auto &evar = tensors.measurements.energy_variance;
    //    if(status.variance_mpo_converged_for > 0 or (evar.has_value() and evar.value() < settings::precision::variance_convergence_threshold)) {
    //        // No need to do expensive operations -- just finish
    //        c1.optMode  = OptMode::VARIANCE;
    //        c1.optSpace = OptSpace::DIRECT;
    //        c1.retry    = false;
    //    }

    if(tensors.state->size_1site() > settings::precision::max_size_part_diag) {
        // Make sure to avoid SUBSPACE if the 1-site problem size is huge
        // When this happens, we can't use SUBSPACE even with two sites,
        // so we might as well give it to DIRECT, which handles larger problems.
        c1.optMode = OptMode::VARIANCE;
        if(c1.optSpace == OptSpace::SUBSPACE) c1.optSpace = OptSpace::DIRECT;
    }

    // Setup the maximum problem size here
    switch(c1.optSpace) {
        case OptSpace::SUBSPACE: c1.max_problem_size = settings::precision::max_size_part_diag; break;
        case OptSpace::DIRECT:
        case OptSpace::KRYLOV: c1.max_problem_size = settings::precision::max_size_direct; break;
    }

    // We can make trials with different number of sites.
    // Eg if the simulation is stuck we may try with more sites.
    if(status.algorithm_has_stuck_for > settings::strategy::max_stuck_iters / 2) c1.max_sites = settings::strategy::multisite_mps_size_max;
    if(status.algorithm_has_succeeded) c1.max_sites = c1.min_sites; // No need to do expensive operations -- just finish

    c1.chosen_sites = tools::finite::multisite::generate_site_list(*tensors.state, c1.max_problem_size, c1.max_sites, c1.min_sites);
    c1.problem_dims = tools::finite::multisite::get_dimensions(*tensors.state, c1.chosen_sites);
    c1.problem_size = tools::finite::multisite::get_problem_size(*tensors.state, c1.chosen_sites);

    if(status.sub_expansion_alpha > 0) {
        // If we are doing 1-site dmrg, then we better use subspace expansion
        if(c1.chosen_sites.size() == 1 and status.sub_expansion_alpha > 0) c1.alpha_expansion = status.sub_expansion_alpha;
        // If we are stuck and enabled subspace expansion when stuck
        if(settings::strategy::expand_subspace_when_stuck and status.algorithm_has_stuck_for > 0 and status.sub_expansion_alpha > 0)
            c1.alpha_expansion = status.sub_expansion_alpha;
    }

    configs.emplace_back(c1);
    if(not c1.retry) return configs;
    /*
     *
     *  Second trial
     *
     *  NOTES
     *      1) OVERLAP does not get a second chance: It is supposed to pick best overlap, not improve variance
     *      2) By default, the result from c1 is used as a starting point for c2.
     *         However, if you change the number of sites, i.e. c2.max_sites,  then you need to start from scratch with the "current state"
     */

    OptConf c2 = c1; // Start with c1 as a baseline
    c2.label   = "c2";
    c2.optInit = OptInit::LAST_RESULT;
    if(c1.optSpace == OptSpace::SUBSPACE and c1.optMode == OptMode::VARIANCE) {
        // I.e. if we did a SUBSPACE run that did not result in better variance, try direct optimization
        // This usually helps to fine-tune the result from subspace, which may have had low subspace quality.
        c2.optWhen  = OptWhen::PREV_FAIL_WORSENED; // Don't worry about the gradient
        c2.optMode  = OptMode::VARIANCE;
        c2.optSpace = OptSpace::DIRECT;
        c2.retry    = false;
        configs.emplace_back(c2);
    } else if(c1.optSpace == OptSpace::DIRECT) {
        // If we did a DIRECT optimization that terminated with gradient too high, try KRYLOV
        c2.optWhen  = OptWhen::PREV_FAIL_GRADIENT | OptWhen::PREV_FAIL_WORSENED;
        c2.optSpace = OptSpace::KRYLOV;
        c2.retry    = false;
        configs.emplace_back(c2);
    } else if(c1.optSpace == OptSpace::KRYLOV) {
        c2.optWhen  = OptWhen::PREV_FAIL_GRADIENT | OptWhen::PREV_FAIL_WORSENED;
        c2.optMode  = OptMode::VARIANCE;
        c2.optSpace = OptSpace::DIRECT;
        c2.retry    = false;
        configs.emplace_back(c2);
    }

    return configs;
}

bool xdmrg::try_again(const std::vector<tools::finite::opt::opt_mps> &results, const xdmrg::OptConf &conf) {
    if(results.empty()) return true;
    bool should = conf.should_proceed(results.back().get_optexit());
    tools::log->trace("Should try again: {} | when {} | previous exit {}", should, flag2str(conf.optWhen), flag2str(results.back().get_optexit()));
    return should;
}

void xdmrg::single_xDMRG_step() {
    using namespace tools::finite;
    using namespace tools::finite::opt;
    auto                                t_step   = tid::tic_scope("step");
    auto                                confList = get_opt_conf_list();
    std::vector<opt_mps>                results;
    std::optional<std::vector<MpsSite>> mps_original = std::nullopt;
    variance_before_step                             = std::nullopt;
    tools::log->debug("Starting xDMRG iter {} | step {} | pos {} | dir {} | confs {}", status.iter, status.step, status.position, status.direction,
                      confList.size());
    for(const auto &[idx_conf, meta] : iter::enumerate(confList)) {
        if(meta.optMode == OptMode::OVERLAP and meta.optSpace != OptSpace::SUBSPACE)
            throw std::logic_error(fmt::format("[OVERLAP] mode can only run in [SUBSPACE]. Got: [{}|{}]", enum2sv(meta.optMode), enum2sv(meta.optSpace)));

        if(not try_again(results, meta)) break;

        // Try activating the sites asked for;
        tensors.activate_sites(meta.chosen_sites);
        if(tensors.active_sites.empty()) continue;

        // Hold the variance before the optimization step for comparison
        if(not variance_before_step) variance_before_step = measure::energy_variance(tensors); // Should just take value from cache

        // Use subspace expansion if alpha_expansion is set
        // Note that this changes the mps and edges adjacent to "tensors.active_sites"
        if(meta.alpha_expansion) {
            auto pos_expanded = tensors.expand_subspace(std::nullopt, status.chi_lim); // nullopt implies a pos query
            if(not mps_original) mps_original = tensors.state->get_mps_sites(pos_expanded);
            tensors.expand_subspace(meta.alpha_expansion, status.chi_lim);
        }

        // Announce the current configuration for optimization
        tools::log->debug("Running meta {}/{}: {} | init {} | mode {} | space {} | type {} | sites {} | dims {} = {} | alpha {:.3e}", idx_conf + 1,
                          confList.size(), meta.label, enum2sv(meta.optInit), enum2sv(meta.optMode), enum2sv(meta.optSpace), enum2sv(meta.optType),
                          meta.chosen_sites, meta.problem_dims, meta.problem_size,
                          (meta.alpha_expansion ? meta.alpha_expansion.value() : std::numeric_limits<double>::quiet_NaN()));
        // Run the optimization
        switch(meta.optInit) {
            case OptInit::CURRENT_STATE: {
                results.emplace_back(opt::find_excited_state(tensors, status, meta));
                break;
            }
            case OptInit::LAST_RESULT: {
                if(results.empty()) throw std::logic_error("There are no previous results to select an initial state");
                results.emplace_back(opt::find_excited_state(tensors, results.back(), status, meta));
                break;
            }
        }

        // Save the neighboring mps that this result is compatible with,
        // so that we may use them if this result turns out to be the winner
        if(meta.alpha_expansion) {
            results.back().set_alpha(meta.alpha_expansion);
            auto pos_expanded         = tensors.expand_subspace(std::nullopt, status.chi_lim); // nullopt implies a pos query
            results.back().mps_backup = tensors.state->get_mps_sites(pos_expanded);            // Backup the mps sites that this run was compatible with
        }

        // Reset the mps to the original if they were backed up earlier
        // This is so that the next meta starts with unchanged mps as neighbors
        if(mps_original) {
            tensors.state->set_mps_sites(mps_original.value());
            tensors.rebuild_edges(); // Make sure the edges are up to date
        }

        // We can now decide if we are happy with the result or not.
        results.back().set_relchange(results.back().get_variance() / variance_before_step.value());

        /* clang-format off */
        meta.optExit = OptExit::SUCCESS;
        double grad_tol = meta.max_grad_tolerance ? meta.max_grad_tolerance.value() : 1e0;
        if(results.back().get_max_grad()       > grad_tol               ) meta.optExit |= OptExit::FAIL_GRADIENT;
        if(results.back().get_relchange()      > 1.001                  ) meta.optExit |= OptExit::FAIL_WORSENED;
        else if(results.back().get_relchange() > 0.999                  ) meta.optExit |= OptExit::FAIL_NOCHANGE;
        results.back().set_optexit(meta.optExit);
        /* clang-format on */

        tools::log->debug(FMT_STRING("Optimization [{}|{}]: {}. Variance change {:8.2e} --> {:8.2e} ({:.3f} %)"), enum2sv(meta.optMode), enum2sv(meta.optSpace),
                          flag2str(meta.optExit), variance_before_step.value(), results.back().get_variance(), results.back().get_relchange() * 100);
        if(results.back().get_relchange() > 1000) tools::log->error("Variance increase by over 1000x: Something is very wrong");
    }

    if(not results.empty()) {
        if(tools::log->level() <= spdlog::level::debug)
            for(const auto &r : results)
                tools::log->debug(FMT_STRING("Candidate: {:<24} | E/L {:<10.6f}| σ²H {:<8.2e}| sites [{:>2}-{:<2}] | "
                                             "alpha {:8.2e} | "
                                             "overlap {:.16f} | "
                                             "{:20} | {} | time {:.4f} ms"),
                                  r.get_name(), r.get_energy_per_site(), r.get_variance(), r.get_sites().front(), r.get_sites().back(), r.get_alpha(),
                                  r.get_overlap(), fmt::format("[{}][{}]", enum2sv(r.get_optmode()), enum2sv(r.get_optspace())), flag2str(r.get_optexit()),
                                  1000 * r.get_time());

        // Sort the results in order of increasing variance
        auto comp_variance = [](const opt_mps &lhs, const opt_mps &rhs) {
            return lhs.get_variance() < rhs.get_variance();
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
        tensors.merge_multisite_mps(winner.get_tensor(), status.chi_lim);
        tensors.rebuild_edges(); // This will only do work if edges were modified, which is the case in 1-site dmrg.
        if(tools::log->level() <= spdlog::level::debug) {
            auto truncation_errors = tensors.state->get_truncation_errors_active();
            tools::log->debug("Truncation errors: {:8.2e}", fmt::join(truncation_errors, ", "));
        }

        if constexpr(settings::debug) {
            auto variance_before_svd = winner.get_variance();
            auto variance_after_svd  = tools::finite::measure::energy_variance(tensors);
            tools::log->debug("Variance check before SVD: {:8.2e}", variance_before_svd);
            tools::log->debug("Variance check after  SVD: {:8.2e}", variance_after_svd);
            tools::log->debug("Variance change from  SVD: {:.16f}%", 100 * variance_after_svd / variance_before_svd);
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
        if constexpr(settings::debug) tensors.assert_validity();
    }
}

void xdmrg::check_convergence() {
    if(not tensors.position_is_inward_edge()) return;
    auto t_con = tid::tic_scope("conv");

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

    if(status.algorithm_converged_for == 0 and status.variance_mpo_saturated_for * status.entanglement_saturated_for == 0 and
       status.algorithm_has_stuck_for != 0)
        throw std::logic_error("Should have zeroed");
    if(status.algorithm_converged_for == 0 and status.variance_mpo_saturated_for * status.entanglement_saturated_for > 0 and
       status.algorithm_has_stuck_for == 0)
        throw std::logic_error("Should not have zeroed");

    status.algorithm_has_succeeded = status.algorithm_converged_for >= settings::strategy::min_converged_iters and
                                     status.algorithm_saturated_for >= settings::strategy::min_saturation_iters;
    status.algorithm_has_to_stop = status.algorithm_has_stuck_for >= settings::strategy::max_stuck_iters;

    tools::log->info("Algorithm report: converged {} | saturated {} | stuck {} | succeeded {} | has to stop {} | var prec limit {:8.2e}",
                     status.algorithm_converged_for, status.algorithm_saturated_for, status.algorithm_has_stuck_for, status.algorithm_has_succeeded,
                     status.algorithm_has_to_stop, status.energy_variance_prec_limit);

    status.algo_stop = AlgorithmStop::NONE;
    if(status.iter >= settings::xdmrg::min_iters and not tensors.model->is_perturbed()) {
        if(status.iter >= settings::xdmrg::max_iters) status.algo_stop = AlgorithmStop::MAX_ITERS;
        if(status.algorithm_has_succeeded) status.algo_stop = AlgorithmStop::SUCCEEDED;
        if(status.algorithm_has_to_stop) status.algo_stop = AlgorithmStop::SATURATED;
        if(status.num_resets > settings::strategy::max_resets) status.algo_stop = AlgorithmStop::MAX_RESET;
        if(status.entanglement_saturated_for > 0 and settings::xdmrg::finish_if_entanglm_saturated) status.algo_stop = AlgorithmStop::SATURATED;
        if(status.variance_mpo_saturated_for > 0 and settings::xdmrg::finish_if_variance_saturated) status.algo_stop = AlgorithmStop::SATURATED;
        if(settings::strategy::randomize_early and excited_state_number == 0 and tensors.state->find_largest_chi() >= 32 and
           tools::finite::measure::energy_variance(tensors) < 1e-4)
            status.algo_stop = AlgorithmStop::RANDOMIZE;
    }
}

void xdmrg::randomize_into_state_in_energy_window(ResetReason reason, StateInit state_type, std::optional<std::string> sector) {
    tools::log->info("Resetting to state in energy window -- reason: {}", enum2sv(reason));
    tools::log->info("Searching for state in normalized energy range: {} +- {}", status.energy_dens_target, status.energy_dens_window);

    status.num_resets++;
    if(reason == ResetReason::SATURATED and status.num_resets > settings::strategy::max_resets) {
        tools::log->info("Not allowed more resets due to saturation: num resets {} > max resets {}", status.num_resets, settings::strategy::max_resets);
        return;
    }
    auto t_rnd             = tid::tic_scope("rnd_state");
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

void xdmrg::find_energy_range() {
    tools::log->trace("Finding energy range");
    auto t_init = tid::tic_scope("init");
    // Here we define a set of tasks for fdmrg in order to produce the lowest and highest energy eigenstates,
    // We don't want it to randomize its own model, so we implant our current model before running the tasks.

    std::deque<fdmrg_task> gs_tasks = {fdmrg_task::INIT_CLEAR_STATUS, fdmrg_task::INIT_BOND_DIM_LIMITS,
                                       fdmrg_task::INIT_WRITE_MODEL,  fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE,
                                       fdmrg_task::FIND_GROUND_STATE, fdmrg_task::POST_WRITE_RESULT};

    std::deque<fdmrg_task> hs_tasks = {fdmrg_task::INIT_CLEAR_STATUS, fdmrg_task::INIT_BOND_DIM_LIMITS, fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE,
                                       fdmrg_task::FIND_HIGHEST_STATE, fdmrg_task::POST_WRITE_RESULT};
    // Find lowest energy state
    {
        auto  t_gs = tid::tic_scope("fDMRG");
        fdmrg fdmrg_gs(h5file);
        *fdmrg_gs.tensors.model = *tensors.model; // Copy the model
        tools::log              = tools::Logger::setLogger(status.algo_type_str() + "-gs", settings::console::loglevel, settings::console::timestamp);
        fdmrg_gs.run_task_list(gs_tasks);
        status.energy_min_per_site = tools::finite::measure::energy_per_site(fdmrg_gs.tensors);
    }

    // Find highest energy state
    {
        auto  t_hs = tid::tic_scope("fDMRG");
        fdmrg fdmrg_hs(h5file);
        *fdmrg_hs.tensors.model = *tensors.model; // Copy the model
        tools::log              = tools::Logger::setLogger(status.algo_type_str() + "-hs", settings::console::loglevel, settings::console::timestamp);
        fdmrg_hs.run_task_list(hs_tasks);
        status.energy_max_per_site = tools::finite::measure::energy_per_site(fdmrg_hs.tensors);
    }

    // Reset our logger
    tools::log = tools::Logger::getLogger(status.algo_type_str());
}

void xdmrg::update_time_step() {
    tools::log->trace("Updating time step");
    status.delta_t = std::complex<double>(1e-6, 0);
}

void xdmrg::create_hamiltonian_gates() {
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

    for(auto pos : list_1site) {
        auto sites = num::range<size_t>(pos, pos + 1);
        auto nbody = {1ul};
        auto spins = tensors.state->get_spin_dims(sites);
        ham_gates_1body.emplace_back(qm::Gate(tensors.model->get_multisite_ham(sites, nbody), sites, spins));
    }
    for(auto pos : list_2site) {
        auto sites = num::range<size_t>(pos, pos + 2);
        auto nbody = {2ul};
        auto spins = tensors.state->get_spin_dims(sites);
        ham_gates_2body.emplace_back(qm::Gate(tensors.model->get_multisite_ham(sites, nbody), sites, spins));
    }
    for(auto pos : list_3site) {
        auto sites = num::range<size_t>(pos, pos + 3);
        auto nbody = {3ul};
        auto spins = tensors.state->get_spin_dims(sites);
        ham_gates_3body.emplace_back(qm::Gate(tensors.model->get_multisite_ham(sites, nbody), sites, spins));
    }
}

void xdmrg::create_time_evolution_gates() {
    // Create the time evolution operators
    if(std::abs(status.delta_t) == 0) update_time_step();
    tools::log->info("Creating time evolution gates");
    time_gates_1site = qm::time::get_time_evolution_gates(status.delta_t, ham_gates_1body);
    time_gates_2site = qm::time::get_time_evolution_gates(status.delta_t, ham_gates_2body);
    time_gates_3site = qm::time::get_time_evolution_gates(status.delta_t, ham_gates_3body);
}
