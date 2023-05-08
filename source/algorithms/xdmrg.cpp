#include "xdmrg.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "fdmrg.h"
#include "general/iter.h"
#include "io/fmt.h"
#include "math/num.h"
#include "math/rnd.h"
#include "qm/time.h"
#include "tensors/edges/EdgesFinite.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include "tools/common/prof.h"
#include "tools/finite/h5.h"
#include "tools/finite/measure.h"
#include "tools/finite/multisite.h"
#include "tools/finite/opt.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"
#include "tools/finite/print.h"

xdmrg::xdmrg(std::shared_ptr<h5pp::File> h5ppFile_) : AlgorithmFinite(std::move(h5ppFile_), AlgorithmType::xDMRG) {
    tools::log->trace("Constructing class_xdmrg");
    tensors.state->set_name(fmt::format("state_{}", excited_state_number));
}

void xdmrg::resume() {
    // Resume can imply many things
    // 1) Resume a simulation which terminated prematurely
    // 2) Resume a previously successful simulation. This may be desireable if the config
    //    wants something that is not present in the file.
    //      a) A certain number of states
    //      b) A state inside a particular energy window
    //      c) The ground or "roof" states
    // To guide the behavior, we check the setting ResumePolicy.

    auto state_prefixes = tools::common::h5::resume::find_state_prefixes(*h5file, status.algo_type, "state_");
    if(state_prefixes.empty()) throw except::state_error("no resumable states were found");
    for(const auto &state_prefix : state_prefixes) {
        tools::log->info("Resuming state [{}]", state_prefix);
        try {
            tools::finite::h5::load::simulation(*h5file, state_prefix, tensors, status, status.algo_type);
        } catch(const except::load_error &le) { continue; }

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
}

void xdmrg::run_default_task_list() {
    std::deque<xdmrg_task> default_task_list = {
        xdmrg_task::INIT_DEFAULT,
        xdmrg_task::FIND_EXCITED_STATE,
        xdmrg_task::POST_DEFAULT,
    };

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
            case xdmrg_task::INIT_BOND_LIMITS: init_bond_dimension_limits(); break;
            case xdmrg_task::INIT_TRNC_LIMITS: init_truncation_error_limits(); break;
            case xdmrg_task::INIT_ENERGY_LIMITS: init_energy_limits(); break;
            case xdmrg_task::INIT_WRITE_MODEL: write_to_file(StorageEvent::MODEL); break;
            case xdmrg_task::INIT_CLEAR_STATUS: status.clear(); break;
            case xdmrg_task::INIT_CLEAR_CONVERGENCE: clear_convergence_status(); break;
            case xdmrg_task::INIT_DEFAULT: run_preprocessing(); break;

            case xdmrg_task::FIND_ENERGY_RANGE: find_energy_range(); break;
            case xdmrg_task::FIND_EXCITED_STATE:
                tensors.state->set_name(fmt::format("state_{}", excited_state_number));
                run_algorithm();
                break;
            case xdmrg_task::POST_WRITE_RESULT: write_to_file(StorageEvent::LAST_STATE); break;
            case xdmrg_task::POST_PRINT_RESULT: print_status_full(); break;
            case xdmrg_task::POST_PRINT_TIMERS: tools::common::timer::print_timers(); break;
            case xdmrg_task::POST_FES_ANALYSIS: run_fes_analysis(); break;
            case xdmrg_task::POST_DEFAULT: run_postprocessing(); break;
            case xdmrg_task::TIMER_RESET: tid::reset("xDMRG"); break;
        }
        task_list.pop_front();
    }
    if(not task_list.empty()) {
        for(auto &task : task_list) tools::log->critical("Unfinished task: {}", enum2sv(task));
        throw except::runtime_error("Simulation ended with unfinished tasks");
    }
}

void xdmrg::init_energy_limits(std::optional<double> energy_density_target, std::optional<double> energy_density_window) {
    if(not energy_density_target) energy_density_target = settings::xdmrg::energy_density_target;
    if(not energy_density_window) energy_density_window = settings::xdmrg::energy_density_window;
    if(energy_density_target.value() < 0.0 or energy_density_target.value() > 1.0)
        throw except::runtime_error(
            fmt::format("Error setting energy density target: Expected value in range [0 - 1.0], got: [{:.8f}]", energy_density_target.value()));
    if(energy_density_window.value() < 0.0 or energy_density_window.value() > 0.5)
        throw except::runtime_error(
            fmt::format("Error setting energy density window: Expected value in range [0 - 0.5], got: [{:.8f}]", energy_density_window.value()));
    status.energy_dens_target = energy_density_target.value();
    status.energy_dens_window = energy_density_window.value();

    // Set energy boundaries. This function is supposed to run after find_energy_range!
    if(status.energy_max == status.energy_min)
        throw except::runtime_error("Could not set energy limits because energy_max == {} and energy_min == {}\n"
                                    "Try running find_energy_range() first",
                                    status.energy_max, status.energy_min);
    status.energy_tgt  = status.energy_min + status.energy_dens_target * (status.energy_max - status.energy_min);
    status.energy_ulim = status.energy_tgt + status.energy_dens_window * (status.energy_max - status.energy_min);
    status.energy_llim = status.energy_tgt - status.energy_dens_window * (status.energy_max - status.energy_min);
    tools::log->info("Energy minimum     = {:.8f}", status.energy_min);
    tools::log->info("Energy maximum     = {:.8f}", status.energy_max);
    tools::log->info("Energy target      = {:.8f}", status.energy_tgt);
    tools::log->info("Energy lower limit = {:.8f}", status.energy_llim);
    tools::log->info("Energy upper limit = {:.8f}", status.energy_ulim);
}

void xdmrg::run_preprocessing() {
    tools::log->info("Running {} preprocessing", status.algo_type_sv());
    auto t_pre = tid::tic_scope("pre");
    status.clear();
    randomize_model(); // First use of random!
    tools::finite::print::model(*tensors.model);
    init_bond_dimension_limits();
    init_truncation_error_limits();
    find_energy_range();
    init_energy_limits();
    if(settings::xdmrg::energy_density_window != 0.5)
        randomize_into_state_in_energy_window(ResetReason::INIT, settings::strategy::initial_state);
    else
        randomize_state(ResetReason::INIT, settings::strategy::initial_state);
    write_to_file(StorageEvent::MODEL);
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
        update_state();
        print_status();
        check_convergence();
        write_to_file();

        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);

        // It's important not to perform the last move, so we break now: that last state would not get optimized
        if(status.algo_stop != AlgorithmStop::NONE) break;

        // Prepare for next step

        // Updating bond dimension must go first since it decides based on truncation error, but a projection+normalize resets truncation.
        update_bond_dimension_limit();   // Will update bond dimension if the state precision is being limited by bond dimension
        update_truncation_error_limit(); // Will update truncation error limit if the state is being truncated
        update_expansion_factor_alpha(); // Will update the subspace expansion factor
        try_projection();                // Tries to project the state to the nearest global spin parity sector along settings::strategy::target_axis
        try_parity_shift();              // This shifts the variance of the opposite spin parity sector, to resolve degeneracy/spectral pairing
        shift_mpo_energy();              // Subtracts the current energy per site E/L from each MPO.
        try_moving_sites();              // Tries to overcome an entanglement barrier by moving sites around the lattice, to optimize non-nearest neighbors
        try_residual_optimization();
        move_center_point();             // Moves the center point AC to the next site and increments status.iter and status.step
        status.wall_time = tid::get_unscoped("t_tot").get_time();
        status.algo_time = t_run->get_time();
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", status.algo_type_sv(), tensors.state->get_name(), status.algo_stop_sv());
    status.algorithm_has_finished = true;
    //    tools::finite::measure::parity_components(*tensors.state, qm::spin::half::sz);
}

void xdmrg::try_residual_optimization() {
    if(not tensors.position_is_inward_edge_left()) return;
    if(status.energy_variance_lowest > 1e-2) return;
    auto t_resopt = tid::tic_scope("residual_optimization");
    auto t_var    = tid::tic_scope("variance");
    auto variance = tools::finite::measure::energy_variance(tensors);
    t_var.toc();
    auto t_res    = tid::tic_scope("residual");
    auto residual = tools::finite::measure::residual_norm_full(*tensors.state, *tensors.model);
    t_res.toc();
    tools::log->info("Variance          : {:.3e} | t = {:.3e}", variance, t_var->get_time());
    tools::log->info("Standard deviation: {:.3e} | t = {:.3e}", std::sqrt(variance), t_var->get_time());
    tools::log->info("Residual norm full: {:.3e} | t = {:.3e}", residual, t_res->get_time());
}

void xdmrg::run_fes_analysis() {
    if(settings::strategy::fes_rate == 0) return;
    tools::log = tools::Logger::setLogger(status.algo_type_str() + "-fes", settings::console::loglevel, settings::console::timestamp);
    tools::log->info("Starting {} finite entanglement scaling analysis with bond size step {} of model [{}] for state [{}]", status.algo_type_sv(),
                     settings::strategy::fes_rate, enum2sv(settings::model::model_type), tensors.state->get_name());
    auto t_fes = tid::tic_scope("fes");
    if(not status.algorithm_has_finished) throw except::logic_error("Finite entanglement scaling analysis can only be done after a finished simulation");
    clear_convergence_status();
    status.fes_is_running = true;
    status.bond_max       = settings::xdmrg::bond_max; // Set to highest
    status.bond_lim       = settings::xdmrg::bond_max; // Set to highest
    tensors.move_center_point_to_inward_edge();
    tensors.activate_sites({tensors.get_position()});
    tensors.rebuild_edges();
    auto bond_max_dim = static_cast<long>(std::pow(2.0, tensors.get_length<double>() / 2));
    status.bond_lim   = std::min(bond_max_dim, status.bond_max);
    while(true) {
        tools::log->trace("Starting xDMRG FES step {}, iter {}, pos {}, dir {}, bond_lim {}, trnc_lim {:.2e}", status.step, status.iter, status.position,
                          status.direction, status.bond_lim, status.trnc_lim);
        update_state();
        status.wall_time       = tid::get_unscoped("t_tot").get_time();
        status.algo_time       = t_fes->get_time();
        auto truncation_errors = tensors.state->get_truncation_errors();
        print_status();
        check_convergence();

        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}, bond_lim {}, trnc_lim {:.2e}", status.step, status.iter, status.position,
                          status.direction, status.bond_lim, status.trnc_lim);

        reduce_bond_dimension_limit(settings::strategy::fes_rate, UpdateWhen::SATURATED, StorageEvent::FES_STATE);
        // It's important not to perform the last move, so we break now: that last state would not get optimized
        if(status.algo_stop != AlgorithmStop::NONE) break;

        move_center_point();

        // Retain the max truncation error, otherwise it is lost in the next pass
        tensors.state->keep_max_truncation_errors(truncation_errors);
    }
    tools::log->info("Finished {} finite entanglement scaling of state [{}] -- stop reason: {}", status.algo_type_sv(), tensors.state->get_name(),
                     status.algo_stop_sv());
    // Reset our logger
    tools::log            = tools::Logger::getLogger(status.algo_type_str());
    status.fes_is_running = false;
}

std::vector<xdmrg::OptMeta> xdmrg::get_opt_conf_list() {
    tools::log->trace("Configuring xDMRG optimization trial");
    std::vector<OptMeta> metas;

    /*
     *
     *  First trial
     *
     */
    OptMeta m1;
    m1.label = "m1";
    m1.retry = false;

    // The first decision is easy. Real or complex optimization
    if(tensors.is_real()) m1.optType = OptType::REAL;

    // Set the default svd limits
    m1.bond_lim = status.bond_lim;
    m1.trnc_lim = status.trnc_lim;

    // Set up a multiplier for number of iterations
    size_t iter_stuck_multiplier = status.algorithm_has_stuck_for > 0 ? settings::solver::iter_stuck_multiplier : 1;

    // Copy settings
    m1.max_sites     = std::min(2ul, settings::strategy::multisite_mps_site_def); // Default is 2-site dmrg, unless we specifically ask for 1-site
    m1.compress_otf  = settings::precision::use_compressed_mpo_squared_otf;
    m1.bfgs_max_iter = settings::solver::bfgs_max_iter * iter_stuck_multiplier;
    m1.bfgs_max_rank = status.algorithm_has_stuck_for == 0 ? 16 : 64; // Tested: around 8-32 seems to be a good compromise,but larger is more precise sometimes.
                                                                      // Overhead goes from 1.2x to 2x computation time at in 8 -> 64
    m1.eigs_iter_max = status.variance_mpo_converged_for > 0 or status.energy_variance_lowest < settings::precision::variance_convergence_threshold
                           ? std::min(settings::solver::eigs_iter_max, 10000ul)          // Avoid running too many iterations when already converged
                           : settings::solver::eigs_iter_max * iter_stuck_multiplier;    // Run as much as it takes before convergence

    m1.eigs_tol = std::clamp(status.energy_variance_lowest,                              // Increase precision as variance decreases
                             settings::solver::eigs_tol_min,                             // From min
                             settings::solver::eigs_tol_max);                            // to max
    if(status.algorithm_has_stuck_for > 0) m1.eigs_tol = settings::solver::eigs_tol_min; // Set to high precision when stuck
    if(status.fes_is_running) m1.eigs_tol = settings::solver::eigs_tol_max;              // No need for high precision during FES.

    m1.eigs_ncv = settings::solver::eigs_ncv;
    // Adjust the maximum number of sites to consider
    if(status.algorithm_has_succeeded)
        m1.max_sites = m1.min_sites; // No need to do expensive operations -- just finish
    else {
        using namespace settings::strategy;
        size_t has_stuck_for = status.algorithm_has_stuck_for;
        size_t saturated_for = status.algorithm_saturated_for * (status.algorithm_converged_for == 0); // Turn on only if non-converged
        switch(multisite_mps_when) {
            case MultisiteWhen::NEVER: break;
            case MultisiteWhen::STUCK: m1.max_sites = std::min(multisite_mps_site_def + has_stuck_for, multisite_mps_site_max); break;
            case MultisiteWhen::SATURATED: m1.max_sites = std::min(multisite_mps_site_def + saturated_for, multisite_mps_site_max); break;
            case MultisiteWhen::ALWAYS: m1.max_sites = std::min(multisite_mps_site_def + saturated_for + 1, multisite_mps_site_max); break;
        }
    }

    // Next we set up the mode at the early stages of the simulation
    // Note that we make stricter requirements as we go down the if-list
    bool prefer_eigs_never          = settings::solver::prefer_eigs_over_bfgs == OptEigs::NEVER;
    bool prefer_eigs_when_stuck     = settings::solver::prefer_eigs_over_bfgs == OptEigs::WHEN_STUCK and status.algorithm_has_stuck_for > 0;
    bool prefer_eigs_when_saturated = settings::solver::prefer_eigs_over_bfgs == OptEigs::WHEN_SATURATED and status.algorithm_saturated_for > 0;
    bool prefer_eigs_always         = settings::solver::prefer_eigs_over_bfgs == OptEigs::ALWAYS;

    if(m1.optSolver == OptSolver::BFGS and not prefer_eigs_never and
       (prefer_eigs_when_saturated or prefer_eigs_when_stuck or prefer_eigs_always or status.fes_is_running)) {
        m1.optMode   = OptMode::VARIANCE;
        m1.optSolver = OptSolver::EIGS;
    }

    if(status.iter < settings::xdmrg::opt_overlap_iters + settings::xdmrg::opt_subspace_iters and settings::xdmrg::opt_subspace_iters > 0) {
        // If early in the simulation, and the bond dimension is small enough we use shift-invert optimization
        m1.optMode   = OptMode::SUBSPACE;
        m1.optSolver = OptSolver::EIGS;
        m1.max_sites = settings::strategy::multisite_mps_site_max;
        if(settings::xdmrg::opt_subspace_bond_lim > 0 and m1.bond_lim > settings::xdmrg::opt_subspace_bond_lim) {
            tools::log->info("Will keep bond dimension back during variance|shift-invert optimization {} -> {}", m1.bond_lim,
                             settings::xdmrg::opt_subspace_bond_lim);
            m1.bond_lim = settings::xdmrg::opt_subspace_bond_lim;
        }
    }

    if(status.iter < settings::xdmrg::opt_overlap_iters) {
        // Very early in the simulation it is worth just following the overlap to get the overall structure of the final state
        m1.optMode   = OptMode::OVERLAP;
        m1.optSolver = OptSolver::EIGS;
        m1.max_sites = settings::strategy::multisite_mps_site_max;
        m1.retry     = false;
        if(settings::xdmrg::opt_overlap_bond_lim > 0 and m1.bond_lim > settings::xdmrg::opt_overlap_bond_lim) {
            tools::log->info("Will keep bond dimension back during overlap optimization {} -> {}", m1.bond_lim, settings::xdmrg::opt_subspace_bond_lim);
            m1.bond_lim = settings::xdmrg::opt_subspace_bond_lim;
        }
    }

    // Setup strong overrides to normal conditions, e.g. when the algorithm has already converged
    if(tensors.state->size_1site() > settings::solver::max_size_shift_invert) {
        // Make sure to avoid size-sensitive optimization modes if the 1-site problem size is huge
        // When this happens, we should use optimize VARIANCE using EIGS or BFGS instead.
        m1.optMode = OptMode::VARIANCE;
    }

    // Set up the maximum problem size here
    switch(m1.optMode) {
        case OptMode::OVERLAP:
        case OptMode::SUBSPACE: m1.max_problem_size = settings::solver::max_size_shift_invert; break;
        case OptMode::ENERGY:
        case OptMode::SIMPS:
        case OptMode::VARIANCE: m1.max_problem_size = settings::precision::max_size_multisite; break;
    }

    m1.chosen_sites = tools::finite::multisite::generate_site_list(*tensors.state, m1.max_problem_size, m1.max_sites, m1.min_sites, "meta 1");
    m1.problem_dims = tools::finite::multisite::get_dimensions(*tensors.state, m1.chosen_sites);
    m1.problem_size = tools::finite::multisite::get_problem_size(*tensors.state, m1.chosen_sites);

    // Do eigs (or eig) instead of bfgs when it's cheap
    if(m1.problem_size <= settings::solver::max_size_full_eigs) m1.optSolver = OptSolver::EIGS;

    // Let EIGS fix BFGS rnorm when enabled
    if(m1.optSolver == OptSolver::BFGS) m1.retry = settings::solver::bfgs_fix_rnorm_w_eigs;

    if(status.env_expansion_alpha > 0 and not status.fes_is_running) {
        // If we are doing 1-site dmrg, then we better use subspace expansion
        if(m1.chosen_sites.size() == 1) m1.alpha_expansion = status.env_expansion_alpha;
    }

    m1.validate();
    metas.emplace_back(m1);
    if(not m1.retry) return metas;
    /*
     *
     *  Second trial
     *
     *  NOTES
     *      1) OVERLAP does not get a second chance: It is supposed to pick best overlap, not improve variance
     *      2) By default, the result from m1 is used as a starting point for m2.
     *         However, if you change the number of sites, i.e. m2.max_sites,  then you need to start from scratch with the "current state"
     */

    OptMeta m2 = m1; // Start with m1 as a baseline
    m2.label   = "m2";
    m2.optInit = OptInit::LAST_RESULT;
    m2.retry   = false;
    if(m1.optMode == OptMode::SUBSPACE) {
        // I.e. if we did a SUBSPACE run that did not result in better variance, try VARIANCE optimization
        // This usually helps to fine-tune the result if the subspace had bad quality.
        m2.optWhen   = OptWhen::PREV_FAIL_WORSENED; // Don't worry about the gradient
        m2.optMode   = OptMode::VARIANCE;
        m2.optSolver = OptSolver::BFGS;
        metas.emplace_back(m2);
    } else if(m1.optSolver == OptSolver::EIGS and m1.optMode == OptMode::ENERGY) {
        // If we did a EIGS|ENERGY optimization that worsened the variance, run EIGS|VARIANCE with the last result as initial state
        m2.optWhen   = OptWhen::PREV_FAIL_WORSENED;
        m2.optSolver = OptSolver::EIGS;
        m2.optMode   = OptMode::VARIANCE;
        m2.optInit   = OptInit::LAST_RESULT;
        metas.emplace_back(m2);
    } else if(m1.optMode == OptMode::VARIANCE and
              (m1.optSolver == OptSolver::EIGS or (m1.optSolver == OptSolver::BFGS and settings::solver::bfgs_fix_rnorm_w_eigs))) {
        // If we did a VARIANCE optimization whose residual_norm did succeed, then run again for longer and more sites
        m2.optWhen   = OptWhen::PREV_FAIL_RESIDUAL | OptWhen::PREV_FAIL_OVERLAP | OptWhen::PREV_FAIL_GRADIENT | OptWhen::PREV_FAIL_WORSENED;
        m2.optSolver = OptSolver::EIGS;
        m2.optMode   = OptMode::VARIANCE;
        m2.optInit   = OptInit::LAST_RESULT;
        metas.emplace_back(m2);
    }
    for(const auto &config : metas) config.validate();
    return metas;
}

bool xdmrg::try_again(const std::vector<tools::finite::opt::opt_mps> &results, const xdmrg::OptMeta &meta) {
    if(results.empty()) return true;
    bool should = meta.should_proceed(results.back().get_optexit());
    if(should) tools::log->debug("Optimizer status: {} --> try again: {} | when: {}", flag2str(results.back().get_optexit()), should, flag2str(meta.optWhen));
    return should;
}

void xdmrg::update_state() {
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
        if(not try_again(results, meta)) break;

        // Try activating the sites asked for;
        tensors.activate_sites(meta.chosen_sites);
        if(tensors.active_sites.empty()) continue;
        tensors.rebuild_edges();

        // Hold the variance before the optimization step for comparison
        if(not variance_before_step) variance_before_step = measure::energy_variance(tensors); // Should just take value from cache

        // Use environment expansion if alpha_expansion is set
        // Note that this changes the mps and edges adjacent to "tensors.active_sites"
        if(meta.alpha_expansion) {
            auto pos_expanded = tensors.expand_environment(std::nullopt, EnvExpandMode::VAR); // nullopt implies a pos query
            if(not mps_original) mps_original = tensors.state->get_mps_sites(pos_expanded);
            tensors.expand_environment(meta.alpha_expansion, EnvExpandMode::VAR, svd::config(meta.bond_lim));
        }

        // Announce the current configuration for optimization
        tools::log->debug("Running meta {}/{}: {} | init {} | mode {} | space {} | type {} | sites {} | dims {} = {} | ε = {:.2e} | α = {:.3e}", idx_conf + 1,
                          confList.size(), meta.label, enum2sv(meta.optInit), enum2sv(meta.optMode), enum2sv(meta.optSolver), enum2sv(meta.optType),
                          meta.chosen_sites, tensors.state->active_dimensions(), tensors.state->active_problem_size(), meta.trnc_lim,
                          (meta.alpha_expansion ? meta.alpha_expansion.value() : std::numeric_limits<double>::quiet_NaN()));
        // Run the optimization
        switch(meta.optInit) {
            case OptInit::CURRENT_STATE: {
                auto initial_state = opt::get_opt_initial_mps(tensors);
                results.emplace_back(opt::find_excited_state(tensors, initial_state, status, meta));
                break;
            }
            case OptInit::LAST_RESULT: {
                if(results.empty()) throw except::logic_error("There are no previous results to select an initial state");
                results.emplace_back(opt::find_excited_state(tensors, results.back(), status, meta));
                break;
            }
        }

        // Save the neighboring mps that this result is compatible with,
        // so that we may use them if this result turns out to be the winner
        if(meta.alpha_expansion) {
            results.back().set_alpha(meta.alpha_expansion);
            auto pos_expanded         = tensors.expand_environment(std::nullopt, EnvExpandMode::VAR); // nullopt implies a pos query
            results.back().mps_backup = tensors.state->get_mps_sites(pos_expanded);                   // Backup the mps sites that this run was compatible with
        }

        // Reset the mps to the original if they were backed up earlier
        // This is so that the next meta starts with unchanged mps as neighbors
        if(mps_original) {
            tensors.state->set_mps_sites(mps_original.value());
            tensors.rebuild_edges(); // Make sure the edges are up-to-date
        }

        // We can now decide if we are happy with the result or not.
        results.back().set_relchange(results.back().get_variance() / variance_before_step.value());
        results.back().set_bond_limit(meta.bond_lim);
        results.back().set_trnc_limit(meta.trnc_lim);
        /* clang-format off */
        meta.optExit = OptExit::SUCCESS;
        if(results.back().get_grad_max()       > 1.000                         ) meta.optExit |= OptExit::FAIL_GRADIENT;
        if(results.back().get_rnorm()          > settings::solver::eigs_tol_max) meta.optExit |= OptExit::FAIL_RESIDUAL;
        if(results.back().get_eigs_nev()       == 0 and
                          meta.optSolver       == OptSolver::EIGS              ) meta.optExit |= OptExit::FAIL_RESIDUAL; // No convergence
        if(results.back().get_overlap()        < 0.010                         ) meta.optExit |= OptExit::FAIL_OVERLAP;
        if(results.back().get_relchange()      > 1.001                         ) meta.optExit |= OptExit::FAIL_WORSENED;
        else if(results.back().get_relchange() > 0.999                         ) meta.optExit |= OptExit::FAIL_NOCHANGE;

        results.back().set_optexit(meta.optExit);
        /* clang-format on */

        tools::log->trace(FMT_STRING("Optimization [{}|{}]: {}. Variance change {:8.2e} --> {:8.2e} ({:.3f} %)"), enum2sv(meta.optMode),
                          enum2sv(meta.optSolver), flag2str(meta.optExit), variance_before_step.value(), results.back().get_variance(),
                          results.back().get_relchange() * 100);
        if(results.back().get_relchange() > 1000) tools::log->error("Variance increase by over 1000x: Something is very wrong");
    }

    if(not results.empty()) {
        if(tools::log->level() <= spdlog::level::debug)
            for(const auto &r : results)
                tools::log->debug(FMT_STRING("Candidate: {:<24} | E {:<20.16f}| σ²H {:<8.2e} | rnorm {:8.2e} | overlap {:.16f} | "
                                             "alpha {:8.2e} | "
                                             "sites {} |"
                                             "{:20} | {} | time {:.2e} s"),
                                  r.get_name(), r.get_energy(), r.get_variance(), r.get_rnorm(), r.get_overlap(), r.get_alpha(), r.get_sites(),
                                  fmt::format("[{}][{}]", enum2sv(r.get_optmode()), enum2sv(r.get_optsolver())), flag2str(r.get_optexit()), r.get_time());

        // Sort the results in order of increasing variance
        auto comp_variance = [](const opt_mps &lhs, const opt_mps &rhs) {
            if(std::isnan(lhs.get_variance())) throw except::logic_error("Result error: state [{}] has variance NAN", lhs.get_name());
            if(std::isnan(rhs.get_variance())) throw except::logic_error("Result error: state [{}] has variance NAN", rhs.get_name());
            return lhs.get_variance() < rhs.get_variance();
        };
        std::sort(results.begin(), results.end(), comp_variance);

        // Take the best result
        const auto &winner = results.front();
        last_optspace      = winner.get_optsolver();
        last_optmode       = winner.get_optmode();
        tensors.activate_sites(winner.get_sites());
        tensors.state->tag_active_sites_normalized(false);
        if(not winner.mps_backup.empty()) {
            // Restore the neighboring mps that are compatible with this winner
            tensors.state->set_mps_sites(winner.mps_backup);
            tensors.rebuild_edges();
        }

        // Do the truncation with SVD
        tensors.merge_multisite_mps(winner.get_tensor(), svd::config(winner.get_bond_lim(), winner.get_trnc_lim()));
        tensors.rebuild_edges(); // This will only do work if edges were modified, which is the case in 1-site dmrg.
        if(tools::log->level() <= spdlog::level::trace)
            tools::log->trace("Truncation errors: {:8.2e}", fmt::join(tensors.state->get_truncation_errors_active(), ", "));

        if constexpr(settings::debug) {
            auto variance_before_svd = winner.get_variance();
            auto variance_after_svd  = tools::finite::measure::energy_variance(tensors);
            tools::log->debug("Variance check before SVD: {:8.2e}", variance_before_svd);
            tools::log->debug("Variance check after  SVD: {:8.2e}", variance_after_svd);
            tools::log->debug("Variance change from  SVD: {:.16f}%", 100 * variance_after_svd / variance_before_svd);
        }

        // Update current energy density ε
        status.energy_dens = (tools::finite::measure::energy_per_site(tensors) - status.energy_min) / (status.energy_max - status.energy_min);

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
    check_convergence_spin_parity_sector(settings::strategy::target_axis);
    if(std::max(status.variance_mpo_saturated_for, status.entanglement_saturated_for) > settings::strategy::max_saturation_iters or
       (status.variance_mpo_saturated_for > 0 and status.entanglement_saturated_for > 0))
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

    status.algorithm_has_succeeded = status.bond_limit_has_reached_max and status.algorithm_converged_for >= settings::strategy::min_converged_iters and
                                     status.algorithm_saturated_for >= settings::strategy::min_saturation_iters;
    status.algorithm_has_to_stop = status.bond_limit_has_reached_max and status.algorithm_has_stuck_for >= settings::strategy::max_stuck_iters;

    tools::log->info(
        "Algorithm report: converged {} (σ² {} Sₑ {} spin {}) | saturated {} (σ² {} Sₑ {}) | stuck {} | succeeded {} | has to stop {} | var prec limit {:8.2e}",
        status.algorithm_converged_for, status.variance_mpo_converged_for, status.entanglement_converged_for, status.spin_parity_has_converged,
        status.algorithm_saturated_for, status.variance_mpo_saturated_for, status.entanglement_saturated_for, status.algorithm_has_stuck_for,
        status.algorithm_has_succeeded, status.algorithm_has_to_stop, status.energy_variance_prec_limit);

    status.algo_stop = AlgorithmStop::NONE;
    if(status.iter >= settings::xdmrg::min_iters) {
        if(status.iter >= settings::xdmrg::max_iters) status.algo_stop = AlgorithmStop::MAX_ITERS;
        if(status.algorithm_has_succeeded) status.algo_stop = AlgorithmStop::SUCCESS;
        if(status.algorithm_has_to_stop) status.algo_stop = AlgorithmStop::SATURATED;
        if(status.num_resets > settings::strategy::max_resets) status.algo_stop = AlgorithmStop::MAX_RESET;
        if(status.entanglement_saturated_for > 0 and settings::xdmrg::finish_if_entanglm_saturated) status.algo_stop = AlgorithmStop::SATURATED;
        if(status.variance_mpo_saturated_for > 0 and settings::xdmrg::finish_if_variance_saturated) status.algo_stop = AlgorithmStop::SATURATED;
        if(settings::strategy::randomize_early and excited_state_number == 0 and tensors.state->find_largest_bond() >= 32 and
           tools::finite::measure::energy_variance(tensors) < 1e-4)
            status.algo_stop = AlgorithmStop::RANDOMIZE;
    }
}

void xdmrg::randomize_into_state_in_energy_window(ResetReason reason, StateInit state_type, const std::optional<std::string> &sector) {
    tools::log->info("Resetting to state in energy window -- reason: {}", enum2sv(reason));
    tools::log->info("Searching for state in normalized energy range: {} +- {}", status.energy_dens_target, status.energy_dens_window);

    status.num_resets++;
    if(reason == ResetReason::SATURATED and status.num_resets > settings::strategy::max_resets) {
        tools::log->info("Not allowed more resets due to saturation: num resets {} > max resets {}", status.num_resets, settings::strategy::max_resets);
        return;
    }
    auto t_rnd             = tid::tic_scope("rnd_state_ewin", tid::level::higher);
    int  counter           = 0;
    bool outside_of_window = true;
    tensors.activate_sites(settings::solver::max_size_full_eigs, 2);
    tensors.rebuild_edges();
    while(true) {
        randomize_state(ResetReason::FIND_WINDOW, state_type, std::nullopt, sector, -1); // Do not use the bitfield: set to -1
        status.energy_dens = tools::finite::measure::energy_normalized(tensors, status.energy_min, status.energy_max);
        outside_of_window  = std::abs(status.energy_dens - status.energy_dens_target) >= status.energy_dens_window;
        tools::log->info("New energy density: {:.16f} | window {} | outside of window: {}", status.energy_dens, status.energy_dens_window, outside_of_window);
        if(not outside_of_window) break;
        counter++;
        if(counter >= 200) throw except::runtime_error("Failed to find initial state in energy window after {}. retries: ", counter);
        if(counter % 10 == 0 and energy_window_growth_factor != 1.0) {
            double old_energy_dens_window = status.energy_dens_window;
            double new_energy_dens_window = std::min(energy_window_growth_factor * status.energy_dens_window, 0.5);

            tools::log->info("Can't find state in energy window.  Increasing energy window: {} --> {}", old_energy_dens_window, new_energy_dens_window);
            status.energy_dens_window = new_energy_dens_window;
        }
    }
    tools::log->info("Energy initial = {:.16f} | density = {:.8f} | retries = {}", tools::finite::measure::energy(tensors), status.energy_dens, counter);
    clear_convergence_status();
    init_energy_limits(std::nullopt, status.energy_dens_window);
    tools::log->info("Number of product state resets: {}", status.num_resets);
}

void xdmrg::find_energy_range() {
    tools::log->trace("Finding energy range");
    auto t_init = tid::tic_scope("init");
    // Here we define a set of tasks for fdmrg in order to produce the lowest and highest energy eigenstates,
    // We don't want it to randomize its own model, so we implant our current model before running the tasks.

    std::deque<fdmrg_task> gs_tasks = {fdmrg_task::INIT_CLEAR_STATUS, fdmrg_task::INIT_BOND_LIMITS, fdmrg_task::INIT_TRNC_LIMITS,
                                       fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE, fdmrg_task::FIND_GROUND_STATE};

    std::deque<fdmrg_task> hs_tasks = {fdmrg_task::INIT_CLEAR_STATUS, fdmrg_task::INIT_BOND_LIMITS, fdmrg_task::INIT_TRNC_LIMITS,
                                       fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE, fdmrg_task::FIND_HIGHEST_STATE};
    // Find the lowest energy state
    {
        auto  t_gs = tid::tic_scope("fDMRG");
        fdmrg fdmrg_gs(h5file);
        *fdmrg_gs.tensors.model = *tensors.model; // Copy the model
        fdmrg_gs.tensors.state->set_name("state_emin");
        tools::log = tools::Logger::setLogger(status.algo_type_str() + "-gs", settings::console::loglevel, settings::console::timestamp);
        fdmrg_gs.run_task_list(gs_tasks);
        status.energy_min = tools::finite::measure::energy(fdmrg_gs.tensors);
        write_to_file(*fdmrg_gs.tensors.state, *fdmrg_gs.tensors.model, *fdmrg_gs.tensors.edges, StorageEvent::EMIN_STATE, CopyPolicy::OFF);
    }

    // Find the highest energy state
    {
        auto  t_hs = tid::tic_scope("fDMRG");
        fdmrg fdmrg_hs(h5file);
        *fdmrg_hs.tensors.model = *tensors.model; // Copy the model
        fdmrg_hs.tensors.state->set_name("state_emax");
        tools::log = tools::Logger::setLogger(status.algo_type_str() + "-hs", settings::console::loglevel, settings::console::timestamp);
        fdmrg_hs.run_task_list(hs_tasks);
        status.energy_max = tools::finite::measure::energy(fdmrg_hs.tensors);
        write_to_file(*fdmrg_hs.tensors.state, *fdmrg_hs.tensors.model, *fdmrg_hs.tensors.edges, StorageEvent::EMAX_STATE, CopyPolicy::OFF);
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
    if(abs_t(status.delta_t) == 0) update_time_step();
    tools::log->info("Creating time evolution gates");
    time_gates_1site = qm::time::get_time_evolution_gates(status.delta_t, ham_gates_1body);
    time_gates_2site = qm::time::get_time_evolution_gates(status.delta_t, ham_gates_2body);
    time_gates_3site = qm::time::get_time_evolution_gates(status.delta_t, ham_gates_3body);
}
