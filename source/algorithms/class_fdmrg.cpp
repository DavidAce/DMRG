//
// Created by david on 2018-01-31.
//

#include "class_fdmrg.h"
#include <config/nmspc_settings.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/io.h>
#include <tools/finite/measure.h>
#include <tools/finite/mpo.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <tools/finite/opt.h>
#include <tools/finite/svd.h>

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
    tools::finite::io::h5resume::load_tensors(*h5pp_file, state_prefix, tensors, status);

    // Our first task is to decide on a state name for the newly loaded state
    // The simplest is to inferr it from the state prefix itself
    auto name   = tools::common::io::h5resume::extract_state_name(state_prefix);


    // Initialize a custom task list
    std::list<fdmrg_task> task_list;

    if(not status.algorithm_has_finished) {
        // This could be a checkpoint state
        // Simply "continue" the algorithm until convergence
        if(name.find("emax") != std::string::npos) task_list.emplace_back(fdmrg_task::FIND_HIGHEST_STATE);
        else if(name.find("emin") != std::string::npos) task_list.emplace_back(fdmrg_task::FIND_GROUND_STATE);
        else throw std::runtime_error(fmt::format("Unrecognized state name for fdmrg: [{}]",name));
        task_list.emplace_back(fdmrg_task::POST_DEFAULT);
        run_task_list(task_list);
    }
    // If we reached this point the current state has finished for one reason or another.
    // TODO: We may still have some more things to do, e.g. the config may be asking for more states
}

void class_fdmrg::run_task_list(std::list<fdmrg_task> &task_list) {
    while(not task_list.empty()) {
        auto task = task_list.front();
        switch(task) {
            case fdmrg_task::INIT_RANDOMIZE_MODEL: randomize_model(); break;
            case fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE: randomize_into_product_state(ResetReason::INIT); break;
            case fdmrg_task::INIT_BOND_DIM_LIMITS: init_bond_dimension_limits(); break;
            case fdmrg_task::INIT_WRITE_MODEL: write_to_file(StorageReason::MODEL); break;
            case fdmrg_task::INIT_CLEAR_STATUS: status.clear(); break;
            case fdmrg_task::INIT_DEFAULT: run_preprocessing(); break;
            case fdmrg_task::FIND_GROUND_STATE:
                ritz = StateRitz::SR;
                run_algorithm();
                break;
            case fdmrg_task::FIND_HIGHEST_STATE:
                ritz = StateRitz::LR;
                run_algorithm();
                break;
            case fdmrg_task::POST_WRITE_RESULT: write_to_file(StorageReason::FINISHED); break;
            case fdmrg_task::POST_PRINT_RESULT: print_status_full(); break;
            case fdmrg_task::POST_DEFAULT: run_postprocessing(); break;
        }
        task_list.pop_front();
    }
}

void class_fdmrg::run_default_task_list() {
    std::list<fdmrg_task> default_task_list = {
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
    tools::common::profile::t_pre->tic();
    status.clear();
    randomize_model(); // First use of random!
    init_bond_dimension_limits();
    randomize_into_product_state(ResetReason::INIT);
    auto spin_components = tools::finite::measure::spin_components(*tensors.state);
    tools::log->info("Initial spin components: {}", spin_components);
    tools::common::profile::t_pre->toc();
    tools::log->info("Finished {} preprocessing", algo_name);
}

void class_fdmrg::run_algorithm() {
    if(state_name.empty()) state_name = ritz == StateRitz::SR ? "state_emin" : "state_emax";
    tools::log->info("Starting {} algorithm with model [{}] for state [{}]", algo_name, enum2str(settings::model::model_type), state_name);
    while(true) {
        single_fdmrg_step();
        print_status_update();
        write_to_file();
        check_convergence();
        update_truncation_limit();     // Will update SVD threshold iff the state precision is being limited by truncation error
        update_bond_dimension_limit(); // Will update bond dimension if the state precision is being limited by bond dimension
        try_projection();

        // It's important not to perform the last move.
        // That last state would not get optimized
        if(tensors.position_is_any_edge()) {
            if(status.iter >= settings::fdmrg::max_iters) {
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
        }
        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);
        move_center_point();
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", algo_name, state_name, enum2str(stop_reason));
    status.algorithm_has_finished = true;
}

void class_fdmrg::single_fdmrg_step() {
    /*!
     * \fn void single_DMRG_step(std::string ritz)
     */
    tools::log->trace("Starting single fdmrg step with ritz [{}]", enum2str(ritz));
    tools::common::profile::t_sim->tic();
    tensors.activate_sites(settings::precision::max_size_part_diag, settings::strategy::multisite_max_sites);
    Eigen::Tensor<Scalar, 3> multisite_tensor = tools::finite::opt::find_ground_state(tensors, ritz);
    tensors.merge_multisite_tensor(multisite_tensor, status.chi_lim);

    // Update record holder
    auto var = tools::finite::measure::energy_variance_per_site(tensors);
    if(var < status.lowest_recorded_variance) status.lowest_recorded_variance = var;

    tools::common::profile::t_sim->toc();
    status.wall_time = tools::common::profile::t_tot->get_age();
    status.simu_time = tools::common::profile::t_sim->get_measured_time();
}

void class_fdmrg::check_convergence() {
    tools::common::profile::t_con->tic();

    if(tensors.position_is_any_edge()) {
        check_convergence_variance();
        check_convergence_entg_entropy();
    }

    status.algorithm_has_converged = status.variance_mpo_has_converged and status.entanglement_has_converged;

    status.algorithm_has_saturated =
        ((status.variance_mpo_saturated_for >= min_saturation_iters and status.entanglement_saturated_for >= min_saturation_iters) or
         (tensors.state->get_iteration() > settings::fdmrg::min_iters and not tensors.state->any_sites_updated()));

    status.algorithm_has_succeeded = status.algorithm_has_converged and status.algorithm_has_saturated;

    status.algorithm_has_got_stuck = status.algorithm_has_saturated and not status.algorithm_has_succeeded;

    if(tensors.state->position_is_any_edge()) { status.algorithm_has_stuck_for = status.algorithm_has_got_stuck ? status.algorithm_has_stuck_for + 1 : 0; }

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

bool   class_fdmrg::cfg_algorithm_is_on() { return settings::fdmrg::on; }
long   class_fdmrg::cfg_chi_lim_max() { return settings::fdmrg::chi_lim_max; }
size_t class_fdmrg::cfg_print_freq() { return settings::fdmrg::print_freq; }
bool   class_fdmrg::cfg_chi_lim_grow() { return settings::fdmrg::chi_lim_grow; }
long   class_fdmrg::cfg_chi_lim_init() { return settings::fdmrg::chi_lim_init; }
bool   class_fdmrg::cfg_store_wave_function() { return settings::fdmrg::store_wavefn; }
