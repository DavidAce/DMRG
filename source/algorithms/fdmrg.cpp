#include "fdmrg.h"
#include <config/settings.h>
#include <io/fmt.h>
#include <tensors/state/StateFinite.h>
#include <tid/tid.h>
#include <tools/common/h5.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/h5.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt.h>
#include <tools/finite/opt_mps.h>

fdmrg::fdmrg(std::shared_ptr<h5pp::File> h5pp_file_) : AlgorithmFinite(std::move(h5pp_file_), AlgorithmType::fDMRG) {
    tools::log->trace("Constructing class_fdmrg");
}

void fdmrg::resume() {
    // Resume can imply many things
    // 1) Resume a simulation which terminated prematurely
    // 2) Resume a previously successful simulation. This may be desireable if the config
    //    wants something that is not present in the file.
    //      a) A certain number of states
    //      b) A state inside of a particular energy window
    //      c) The ground or "roof" states
    // To guide the behavior, we check the setting ResumePolicy.

    auto state_prefix = tools::common::h5::resume::find_resumable_state(*h5pp_file, status.algo_type);
    if(state_prefix.empty()) throw std::runtime_error("Could not resume: no valid state candidates found for resume");
    tools::log->info("Resuming state [{}]", state_prefix);
    tools::finite::h5::load::simulation(*h5pp_file, state_prefix, tensors, status, status.algo_type);
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
            throw std::runtime_error(fmt::format("Unrecognized state name for fdmrg: [{}]", name));
        task_list.emplace_back(fdmrg_task::POST_DEFAULT);
    }
    run_task_list(task_list);
}

void fdmrg::run_task_list(std::deque<fdmrg_task> &task_list) {
    while(not task_list.empty()) {
        auto task = task_list.front();
        switch(task) {
            case fdmrg_task::INIT_RANDOMIZE_MODEL: randomize_model(); break;
            case fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE: randomize_state(ResetReason::INIT, StateInit::RANDOM_PRODUCT_STATE); break;
            case fdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE: randomize_state(ResetReason::INIT, StateInit::RANDOM_ENTANGLED_STATE); break;
            case fdmrg_task::INIT_BOND_DIM_LIMITS: init_bond_dimension_limits(); break;
            case fdmrg_task::INIT_WRITE_MODEL: write_to_file(StorageReason::MODEL); break;
            case fdmrg_task::INIT_CLEAR_STATUS: status.clear(); break;
            case fdmrg_task::INIT_CLEAR_CONVERGENCE: clear_convergence_status(); break;
            case fdmrg_task::INIT_DEFAULT: run_preprocessing(); break;
            case fdmrg_task::FIND_GROUND_STATE:
                ritz = StateRitz::SR;
                tensors.state->set_name("state_emin");
                run_algorithm();
                break;
            case fdmrg_task::FIND_HIGHEST_STATE:
                ritz = StateRitz::LR;
                tensors.state->set_name("state_emax");
                run_algorithm();
                break;
            case fdmrg_task::POST_WRITE_RESULT: write_to_file(StorageReason::FINISHED); break;
            case fdmrg_task::POST_PRINT_RESULT: print_status_full(); break;
            case fdmrg_task::POST_PRINT_PROFILING: tools::common::profile::print_profiling(status); break;
            case fdmrg_task::POST_DEFAULT: run_postprocessing(); break;
            case fdmrg_task::PROF_RESET: tid::reset("fDMRG"); break;
        }
        task_list.pop_front();
    }
}

void fdmrg::run_default_task_list() {
    std::deque<fdmrg_task> default_task_list = {
        fdmrg_task::INIT_DEFAULT,
        fdmrg_task::FIND_GROUND_STATE,
        fdmrg_task::POST_DEFAULT,
    };

    run_task_list(default_task_list);
    if(not default_task_list.empty()) {
        for(auto &task : default_task_list) tools::log->critical("Unfinished task: {}", enum2str(task));
        throw std::runtime_error("Simulation ended with unfinished tasks");
    }
}

void fdmrg::run_preprocessing() {
    tools::log->info("Running {} preprocessing", status.algo_type_sv());
    status.clear();
    randomize_model(); // First use of random!
    init_bond_dimension_limits();
    randomize_state(ResetReason::INIT, settings::strategy::initial_state);
    tools::log->info("Finished {} preprocessing", status.algo_type_sv());
}

void fdmrg::run_algorithm() {
    if(tensors.state->get_name().empty()) {
        if(ritz == StateRitz::SR)
            tensors.state->set_name("state_emin");
        else
            tensors.state->set_name("state_emax");
    }
    tools::log->info("Starting {} algorithm with model [{}] for state [{}]", status.algo_type_sv(), enum2str(settings::model::model_type),
                     tensors.state->get_name());
    auto t_run       = tid::tic_scope("run");
    status.algo_stop = AlgorithmStop::NONE;

    while(true) {
        single_fdmrg_step();
        check_convergence();
        write_to_file();
        print_status_update();

        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);

        // It's important not to perform the last move, so we break now: that last state would not get optimized
        if(status.algo_stop != AlgorithmStop::NONE) break;
        update_bond_dimension_limit(); // Will update bond dimension if the state precision is being limited by bond dimension
        try_projection();
        reduce_mpo_energy();
        move_center_point();
        status.wall_time = tid::get_unscoped("t_tot").get_time();
        status.algo_time = t_run->get_time();
        for(const auto &[key, u] : tid::internal::tid_db)
            for(const auto &t : tid::get_tree(u)) tools::log->info("{}", t.str());
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", status.algo_type_sv(), tensors.state->get_name(), status.algo_stop_sv());
    status.algorithm_has_finished = true;
}

void fdmrg::single_fdmrg_step() {
    /*!
     * \fn void single_DMRG_step(std::string ritz)
     */
    auto t_step = tid::tic_scope("step");
    tools::log->debug("Starting fDMRG iter {} | step {} | pos {} | dir {} | ritz {}", status.iter, status.step, status.position, status.direction,
                      enum2str(ritz));
    tensors.activate_sites(settings::precision::max_size_part_diag, 2);

    if(tensors.active_sites.empty())
        tensors.activate_sites({0}); // Activate a site so that edge checks can happen
    else {
        // Decide to use subspace expansion
        if(tensors.active_sites.size() == 1 or
           (status.algorithm_has_stuck_for > 0 and tools::finite::measure::energy_variance(tensors) > settings::precision::variance_convergence_threshold))
            sub_expansion_alpha = std::min(0.1, status.energy_variance_lowest);
        else
            sub_expansion_alpha = std::nullopt;

        // Use subspace expansion if alpha_expansion was set
        if(sub_expansion_alpha) tensors.expand_subspace(sub_expansion_alpha.value(), status.chi_lim);

        auto multisite_mps = tools::finite::opt::find_ground_state(tensors, status, ritz);
        if constexpr(settings::debug)
            tools::log->debug("Variance after opt: {:.8f} | norm {:.16f}", std::log10(multisite_mps.get_variance()), multisite_mps.get_norm());

        tensors.merge_multisite_tensor(multisite_mps.get_tensor(), status.chi_lim);
        if constexpr(settings::debug)
            tools::log->debug("Variance after svd: {:.8f} | trunc: {}", std::log10(tools::finite::measure::energy_variance(tensors)),
                              tools::finite::measure::truncation_errors_active(*tensors.state));
        // Update record holder
        if(not tensors.active_sites.empty()) {
            auto var = tools::finite::measure::energy_variance(tensors);
            tools::log->trace("Updating variance record holder: var {} | record {}", std::log10(var), std::log10(status.energy_variance_lowest));
            if(var < status.energy_variance_lowest) status.energy_variance_lowest = var;
        }
    }
}

void fdmrg::check_convergence() {
    if(not tensors.position_is_inward_edge()) return;
    auto t_con = tid::tic_scope("convergence");
    update_variance_max_digits();
    check_convergence_variance();
    check_convergence_entg_entropy();
    check_convergence_spin_parity_sector(settings::strategy::target_sector);

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

    status.algorithm_has_succeeded = status.algorithm_converged_for > min_converged_iters and status.algorithm_saturated_for > min_saturation_iters;
    status.algorithm_has_to_stop   = status.algorithm_has_stuck_for >= max_stuck_iters;

    if(status.algorithm_converged_for == 0 and status.variance_mpo_saturated_for * status.entanglement_saturated_for == 0 and
       status.algorithm_has_stuck_for != 0)
        throw std::logic_error("Should have zeroed");
    if(status.algorithm_converged_for == 0 and status.variance_mpo_saturated_for * status.entanglement_saturated_for > 0 and
       status.algorithm_has_stuck_for == 0)
        throw std::logic_error("Should not have zeroed");

    tools::log->info("Algorithm report: converged {} | saturated {} | stuck {} | succeeded {} | has to stop {} | var prec limit {:.6f}",
                     status.algorithm_converged_for, status.algorithm_saturated_for, status.algorithm_has_stuck_for, status.algorithm_has_succeeded,
                     status.algorithm_has_to_stop, std::log10(status.energy_variance_prec_limit));
    status.algo_stop = AlgorithmStop::NONE;
    if(status.iter >= settings::fdmrg::min_iters) {
        if(status.iter >= settings::fdmrg::max_iters) status.algo_stop = AlgorithmStop::MAX_ITERS;
        if(status.algorithm_has_succeeded) status.algo_stop = AlgorithmStop::SUCCEEDED;
        if(status.algorithm_has_to_stop) status.algo_stop = AlgorithmStop::SATURATED;
        if(status.num_resets > settings::strategy::max_resets) status.algo_stop = AlgorithmStop::MAX_RESET;
    }
}
