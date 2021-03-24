//
// Created by david on 2018-01-31.
//

#include "class_fdmrg.h"
#include <config/nmspc_settings.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/io.h>
#include <tools/finite/measure.h>
#include <tools/finite/ops.h>
#include <tools/finite/opt.h>
#include <unsupported/Eigen/CXX11/Tensor>

class_fdmrg::class_fdmrg(std::shared_ptr<h5pp::File> h5pp_file_) : class_algorithm_finite(std::move(h5pp_file_), AlgorithmType::fDMRG) {
    tools::log->trace("Constructing class_fdmrg");
}

void class_fdmrg::resume() {
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
    tools::finite::io::h5resume::load_simulation(*h5pp_file, state_prefix, tensors, status);
    // Our first task is to decide on a state name for the newly loaded state
    // The simplest is to inferr it from the state prefix itself
    auto name = tools::common::io::h5resume::extract_state_name(state_prefix);

    // Initialize a custom task list
    std::deque<fdmrg_task> task_list;

    if(status.algorithm_has_succeeded)
        task_list = {fdmrg_task::POST_PRINT_RESULT};
    else{
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

void class_fdmrg::run_task_list(std::deque<fdmrg_task> &task_list) {
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
            case fdmrg_task::POST_PRINT_PROFILING: tools::common::profile::print_profiling(algo_type); break;
            case fdmrg_task::POST_DEFAULT: run_postprocessing(); break;
            case fdmrg_task::PROF_RESET: tools::common::profile::reset_profiling(algo_type); break;
        }
        task_list.pop_front();
    }
}

void class_fdmrg::run_default_task_list() {
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

void class_fdmrg::run_preprocessing() {
    tools::log->info("Running {} preprocessing", algo_name);
    auto t_pre = tools::common::profile::prof[algo_type]["t_pre"]->tic_token();
    status.clear();
    randomize_model(); // First use of random!
    init_bond_dimension_limits();
    randomize_state(ResetReason::INIT, settings::strategy::initial_state);
    tools::log->info("Finished {} preprocessing", algo_name);
}

void class_fdmrg::run_algorithm() {
    if(tensors.state->get_name().empty()) {
        if(ritz == StateRitz::SR)
            tensors.state->set_name("state_emin");
        else
            tensors.state->set_name("state_emax");
    }
    tools::log->info("Starting {} algorithm with model [{}] for state [{}]", algo_name, enum2str(settings::model::model_type), tensors.state->get_name());
    auto t_sim = tools::common::profile::prof[algo_type]["t_sim"]->tic_token();
    stop_reason = StopReason::NONE;
    while(true) {
        single_fdmrg_step();
        check_convergence();
        write_to_file();
        print_status_update();
        print_profiling_lap();

        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);

        // It's important not to perform the last move, so we break now: that last state would not get optimized
        if(stop_reason != StopReason::NONE) break;

        update_bond_dimension_limit(); // Will update bond dimension if the state precision is being limited by bond dimension
        try_projection();
        reduce_mpo_energy();
        move_center_point();
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", algo_name, tensors.state->get_name(), enum2str(stop_reason));
    status.algorithm_has_finished = true;
}

void class_fdmrg::single_fdmrg_step() {
    /*!
     * \fn void single_DMRG_step(std::string ritz)
     */
    tools::log->trace("Starting single fdmrg step with ritz [{}]", enum2str(ritz));
    tensors.activate_sites(settings::precision::max_size_part_diag, 2);
    std::optional<double> alpha_expansion;
    if(tensors.active_sites.size() == 1 or status.algorithm_has_stuck_for == 1) alpha_expansion = std::min(0.1,status.energy_variance_lowest);

    if(tensors.active_sites.empty())
        tensors.activate_sites({0}); // Activate a site so that edge checks can happen
    else {
        // Use subspace expansion if alpha_expansion is set
        if(alpha_expansion)
            tensors.expand_subspace(alpha_expansion.value(), status.chi_lim);

        Eigen::Tensor<Scalar, 3> multisite_tensor = tools::finite::opt::find_ground_state(tensors, ritz);
        if constexpr(settings::debug)
            tools::log->debug("Variance after opt: {:.8f}", std::log10(tools::finite::measure::energy_variance(multisite_tensor, tensors)));

        tensors.merge_multisite_tensor(multisite_tensor, status.chi_lim);
        if constexpr(settings::debug)
            tools::log->debug("Variance after svd: {:.8f} | trunc: {}", std::log10(tools::finite::measure::energy_variance(tensors)),
                              tools::finite::measure::truncation_errors_active(*tensors.state));
        // Update record holder
        if(not tensors.active_sites.empty()) {
            tools::log->trace("Updating variance record holder");
            auto var = tools::finite::measure::energy_variance(tensors);
            if(var < status.energy_variance_lowest) status.energy_variance_lowest = var;
        }
    }
    status.wall_time = tools::common::profile::t_tot->get_measured_time();
    status.algo_time = tools::common::profile::prof[algo_type]["t_sim"]->get_measured_time();
}

void class_fdmrg::check_convergence() {
    if(not tensors.position_is_inward_edge()) return;

    auto t_con = tools::common::profile::prof[algo_type]["t_con"]->tic_token();

    check_convergence_variance();
    check_convergence_entg_entropy();
    check_convergence_spin_parity_sector(settings::strategy::target_sector);

    if(std::max(status.variance_mpo_saturated_for, status.entanglement_saturated_for) > max_saturation_iters or
       (status.variance_mpo_saturated_for > 0 and status.entanglement_saturated_for > 0)) status.algorithm_saturated_for++;
    else status.algorithm_saturated_for = 0;

    if(status.variance_mpo_converged_for > 0 and status.entanglement_converged_for > 0 and status.spin_parity_has_converged) status.algorithm_converged_for++;
    else status.algorithm_converged_for = 0;

    if(status.algorithm_saturated_for > 0 and status.algorithm_converged_for == 0) status.algorithm_has_stuck_for++;
    else status.algorithm_has_stuck_for = 0;

    status.algorithm_has_succeeded = status.algorithm_converged_for > min_converged_iters and status.algorithm_saturated_for > min_saturation_iters;
    status.algorithm_has_to_stop = status.algorithm_has_stuck_for >= max_stuck_iters;

    tools::log->info("Algorithm report: converged {} | saturated {} | stuck {} | succeeded {} | has to stop {}",
                      status.algorithm_converged_for, status.algorithm_saturated_for,  status.algorithm_has_stuck_for, status.algorithm_has_succeeded,
                      status.algorithm_has_to_stop);
    stop_reason = StopReason::NONE;
    if(status.iter >= settings::fdmrg::min_iters) {
        if(status.iter >= settings::fdmrg::max_iters) stop_reason = StopReason::MAX_ITERS;
        if(status.algorithm_has_succeeded) stop_reason = StopReason::SUCCEEDED;
        if(status.algorithm_has_to_stop) stop_reason = StopReason::SATURATED;
        if(status.num_resets > settings::strategy::max_resets) stop_reason = StopReason::MAX_RESET;
    }
}

bool   class_fdmrg::cfg_algorithm_is_on() { return settings::fdmrg::on; }
long   class_fdmrg::cfg_chi_lim_max() { return settings::fdmrg::chi_lim_max; }
size_t class_fdmrg::cfg_print_freq() { return settings::fdmrg::print_freq; }
bool   class_fdmrg::cfg_chi_lim_grow() { return settings::fdmrg::chi_lim_grow; }
long   class_fdmrg::cfg_chi_lim_init() { return settings::fdmrg::chi_lim_init; }
bool   class_fdmrg::cfg_store_wave_function() { return settings::fdmrg::store_wavefn; }
