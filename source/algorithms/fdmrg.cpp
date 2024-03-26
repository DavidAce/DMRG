#include "fdmrg.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "io/fmt_custom.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include "tools/common/prof.h"
#include "tools/finite/h5.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"
#include "tools/finite/print.h"

fdmrg::fdmrg() : AlgorithmFinite(AlgorithmType::fDMRG) { tools::log->trace("Constructing class_fdmrg (without a file)"); }
fdmrg::fdmrg(std::shared_ptr<h5pp::File> h5file_) : AlgorithmFinite(std::move(h5file_), AlgorithmType::fDMRG) { tools::log->trace("Constructing class_fdmrg"); }

std::string_view fdmrg::get_state_name() const {
    if(ritz == OptRitz::SR)
        return "state_emin";
    else
        return "state_emax";
}

void fdmrg::resume() {
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
        // The simplest is to inferr it from the state prefix itself
        auto name = tools::common::h5::resume::extract_state_name(state_prefix);

        // Initialize a custom task list
        std::deque<fdmrg_task> task_list;

        if(status.algorithm_has_succeeded)
            task_list = {fdmrg_task::POST_PRINT_RESULT};
        else {
            task_list.emplace_back(fdmrg_task::INIT_CLEAR_CONVERGENCE);
            // This could be a savepoint state
            // Simply "continue" the algorithm until convergence
            if(name.find("emax") != std::string::npos)
                task_list.emplace_back(fdmrg_task::FIND_HIGHEST_STATE);
            else if(name.find("emin") != std::string::npos)
                task_list.emplace_back(fdmrg_task::FIND_GROUND_STATE);
            else
                throw except::runtime_error("Unrecognized state name for fdmrg: [{}]", name);
            task_list.emplace_back(fdmrg_task::POST_DEFAULT);
        }
        run_task_list(task_list);
    }
}

void fdmrg::run_task_list(std::deque<fdmrg_task> &task_list) {
    while(not task_list.empty()) {
        auto task = task_list.front();
        switch(task) {
            case fdmrg_task::INIT_RANDOMIZE_MODEL: initialize_model(); break;
            case fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE: initialize_state(ResetReason::INIT, StateInit::RANDOM_PRODUCT_STATE); break;
            case fdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE: initialize_state(ResetReason::INIT, StateInit::RANDOM_ENTANGLED_STATE); break;
            case fdmrg_task::INIT_BOND_LIMITS: init_bond_dimension_limits(); break;
            case fdmrg_task::INIT_TRNC_LIMITS: init_truncation_error_limits(); break;
            case fdmrg_task::INIT_WRITE_MODEL: write_to_file(StorageEvent::MODEL); break;
            case fdmrg_task::INIT_CLEAR_STATUS: status.clear(); break;
            case fdmrg_task::INIT_CLEAR_CONVERGENCE: clear_convergence_status(); break;
            case fdmrg_task::INIT_DEFAULT: run_preprocessing(); break;
            case fdmrg_task::FIND_GROUND_STATE:
                ritz = OptRitz::SR;
                tensors.state->set_name("state_emin");
                run_algorithm();
                break;
            case fdmrg_task::FIND_HIGHEST_STATE:
                ritz = OptRitz::LR;
                tensors.state->set_name("state_emax");
                run_algorithm();
                break;
            case fdmrg_task::POST_WRITE_RESULT: write_to_file(StorageEvent::FINISHED); break;
            case fdmrg_task::POST_PRINT_RESULT: print_status_full(); break;
            case fdmrg_task::POST_PRINT_TIMERS: tools::common::timer::print_timers(); break;
            case fdmrg_task::POST_FES_ANALYSIS: run_fes_analysis(); break;
            case fdmrg_task::POST_DEFAULT: run_postprocessing(); break;
            case fdmrg_task::TIMER_RESET: tid::reset("fDMRG"); break;
        }
        task_list.pop_front();
    }
}

void fdmrg::run_default_task_list() {
    fdmrg_task fdmrg_task_find_state_ritz;
    switch(settings::fdmrg::ritz) {
        case OptRitz::SR: fdmrg_task_find_state_ritz = fdmrg_task::FIND_GROUND_STATE; break;
        case OptRitz::LR: fdmrg_task_find_state_ritz = fdmrg_task::FIND_HIGHEST_STATE; break;
        default: throw except::logic_error("fdmrg expects ritz SR or LR. Got: {}", enum2sv(settings::fdmrg::ritz));
    }

    std::deque<fdmrg_task> default_task_list = {
        fdmrg_task::INIT_DEFAULT,
        fdmrg_task_find_state_ritz,
        fdmrg_task::POST_DEFAULT,
    };

    run_task_list(default_task_list);
    if(not default_task_list.empty()) {
        for(auto &task : default_task_list) tools::log->critical("Unfinished task: {}", enum2sv(task));
        throw except::runtime_error("Simulation ended with unfinished tasks");
    }
}

void fdmrg::run_preprocessing() {
    tools::log->info("Running {} preprocessing", status.algo_type_sv());
    auto t_pre = tid::tic_scope("pre");
    status.clear();
    if(tensors.state->get_name().empty()) tensors.state->set_name(get_state_name());
    initialize_model(); // First use of random!
    tools::finite::print::model(*tensors.model);
    init_bond_dimension_limits();
    init_truncation_error_limits();
    initialize_state(ResetReason::INIT, settings::strategy::initial_state);
    set_parity_shift_mpo(); // This shifts the energy of the opposite spin parity sector, to resolve degeneracy/spectral pairing
    set_energy_shift_mpo();
    rebuild_tensors();
    tools::log->info("Finished {} preprocessing", status.algo_type_sv());
}

void fdmrg::run_algorithm() {
    if(tensors.state->get_name().empty()) tensors.state->set_name(get_state_name());
    tools::log->info("Starting {} algorithm with model [{}] for state [{}]", status.algo_type_sv(), enum2sv(settings::model::model_type),
                     tensors.state->get_name());
    auto t_run       = tid::tic_scope("run");
    status.algo_stop = AlgorithmStop::NONE;
    while(true) {
        update_state();
        print_status();
        check_convergence();
        write_to_file();

        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);

        // It's important not to perform the last move, so we break now: that last state would not get optimized
        if(status.algo_stop != AlgorithmStop::NONE) break;
        update_bond_dimension_limit();   // Will update bond dimension if the state precision is being limited by bond dimension
        update_truncation_error_limit(); // Will update truncation error limit if the state is being truncated
        update_expansion_factor_alpha(); // Will update the subspace expansion factor
        try_projection();
        move_center_point();
        status.wall_time = tid::get_unscoped("t_tot").get_time();
        status.algo_time = t_run->get_time();
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", status.algo_type_sv(), tensors.state->get_name(), status.algo_stop_sv());
    status.algorithm_has_finished = true;
    Eigen::Tensor<real, 1> vec    = tools::finite::measure::mps2tensor(*tensors.state).real();
    write_to_file(vec, "vec", StorageEvent::FINISHED);
}

void fdmrg::run_fes_analysis() {
    if(settings::strategy::fes_rate == 0) return;
    tools::log->warn("FES is not yet implemented for fdmrg");
}

void fdmrg::update_state() {
    auto    t_step = tid::tic_scope("step");
    OptMeta meta(ritz, OptFunc::ENERGY);
    if(tensors.is_real()) meta.optType = OptType::REAL; // Can do everything in real mode if the model is real
    std::optional<double> alpha_expansion = std::nullopt;
    tools::log->debug("Starting fDMRG iter {} | step {} | pos {} | dir {} | ritz {} | type {}", status.iter, status.step, status.position, status.direction,
                      enum2sv(ritz), enum2sv(meta.optType));
    tensors.activate_sites(settings::solver::eigs_max_size_shift_invert, settings::strategy::multisite_opt_site_def);
    if(not tensors.active_sites.empty()) {
        tensors.rebuild_edges();
        if(status.env_expansion_alpha > 0) {
            // If we are doing 1-site dmrg, then we better use subspace expansion
            if(tensors.active_sites.size() == 1) alpha_expansion = status.env_expansion_alpha;
            // Use subspace expansion if alpha_expansion was set
            if(alpha_expansion) tensors.expand_environment(alpha_expansion.value(), EnvExpandMode::ENE, svd::config(status.bond_lim));
        }
        auto initial_mps = tools::finite::opt::get_opt_initial_mps(tensors);
        auto result_mps  = tools::finite::opt::find_ground_state(tensors, initial_mps, status, meta);
        if constexpr(settings::debug) tools::log->debug("Variance after opt: {:8.2e} | norm {:.16f}", result_mps.get_variance(), result_mps.get_norm());
        tensors.merge_multisite_mps(result_mps.get_tensor(), svd::config(status.bond_lim, status.trnc_lim));
        tensors.rebuild_edges(); // This will only do work if edges were modified, which is the case in 1-site dmrg.
        if constexpr(settings::debug)
            tools::log->debug("Variance after svd: {:8.2e} | trunc: {}", tools::finite::measure::energy_variance(tensors),
                              tools::finite::measure::truncation_errors_active(*tensors.state));
        // Update record holder
        if(not tensors.active_sites.empty()) {
            auto var = tools::finite::measure::energy_variance(tensors);
            tools::log->trace("Updating variance record holder: var {:8.2e} | record {:8.2e}", var, status.energy_variance_lowest);
            if(var < status.energy_variance_lowest) status.energy_variance_lowest = var;
        }
        if constexpr(settings::debug) tensors.assert_validity();
    }
}

void fdmrg::check_convergence() {
    if(not tensors.position_is_inward_edge()) return;
    auto t_con = tid::tic_scope("convergence");
    update_variance_max_digits();
    check_convergence_variance();
    check_convergence_entg_entropy();
    check_convergence_spin_parity_sector(settings::strategy::target_axis);

    if(std::max(status.variance_mpo_saturated_for, status.entanglement_saturated_for) > settings::strategy::max_saturated_iters or
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
                                     status.algorithm_saturated_for >= settings::strategy::min_saturated_iters;
    status.algorithm_has_to_stop = status.bond_limit_has_reached_max and status.algorithm_has_stuck_for >= settings::strategy::max_stuck_iters;

    tools::log->info(
        "Sweep report: converged {} (σ² {} Sₑ {} spin {}) | saturated {} (σ² {} Sₑ {}) | stuck {} | succeeded: {} | has to stop: {} | σ²H target {:8.2e}",
        status.algorithm_converged_for, status.variance_mpo_converged_for, status.entanglement_converged_for, status.spin_parity_has_converged,
        status.algorithm_saturated_for, status.variance_mpo_saturated_for, status.entanglement_saturated_for, status.algorithm_has_stuck_for,
        status.algorithm_has_succeeded, status.algorithm_has_to_stop,
        std::max({settings::precision::variance_convergence_threshold, status.energy_variance_prec_limit}));
    status.algo_stop = AlgorithmStop::NONE;
    if(status.iter >= settings::fdmrg::min_iters) {
        if(status.iter >= settings::fdmrg::max_iters) status.algo_stop = AlgorithmStop::MAX_ITERS;
        if(status.algorithm_has_succeeded) status.algo_stop = AlgorithmStop::SUCCESS;
        if(status.algorithm_has_to_stop) status.algo_stop = AlgorithmStop::SATURATED;
        if(status.num_resets > settings::strategy::max_resets) status.algo_stop = AlgorithmStop::MAX_RESET;
    }
}
