//
// Created by david on 2018-01-31.
//

#include "class_fdmrg.h"
#include <config/nmspc_settings.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/mpo.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <tools/finite/opt.h>
#include <tools/finite/svd.h>

class_fdmrg::class_fdmrg(std::shared_ptr<h5pp::File> h5pp_file_) : class_algorithm_finite(std::move(h5pp_file_), AlgorithmType::fDMRG) {
    tools::log->trace("Constructing class_fdmrg");
}

void class_fdmrg::run_task_list(std::list<fdmrg_task> &task_list) {
    while(not task_list.empty()) {
        auto task = task_list.front();
        switch(task) {
            case fdmrg_task::INIT_RANDOMIZE_MODEL: randomize_model(); break;
            case fdmrg_task::INIT_RANDOM_PRODUCT_STATE: reset_to_random_product_state(ResetReason::INIT); break;
            case fdmrg_task::INIT_BOND_DIM_LIMITS: reset_bond_dimension_limits(); break;
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

void class_fdmrg::run_algorithm() {
    if(state_name.empty())
        state_name = ritz == StateRitz::SR ? "state_emin" : "state_emax";
    tools::log->info("Starting {} algorithm with model [{}] for state [{}]", algo_name, enum2str(settings::model::model_type), state_name);
    while(true) {
        single_fDMRG_step();
        print_status_update();
        write_to_file();
        copy_from_tmp();
        check_convergence();
        update_truncation_limit();     // Will update SVD threshold iff the state precision is being limited by truncation error
        update_bond_dimension_limit(); // Will update bond dimension if the state precision is being limited by bond dimension
        try_projection(*tensors.state);

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
    tools::log->info("Finished {} simulation -- reason: {}", algo_name, enum2str(stop_reason));
    status.algorithm_has_finished = true;
}

void class_fdmrg::single_fDMRG_step() {
    /*!
     * \fn void single_DMRG_step(std::string ritz)
     */
    tools::log->trace("Starting single fDMRG step with ritz [{}]", enum2str(ritz));
    tools::common::profile::t_sim->tic();
    tensors.activate_sites(settings::precision::max_size_part_diag, settings::strategy::multisite_max_sites);
    Eigen::Tensor<Scalar, 3> multisite_tensor = tools::finite::opt::find_ground_state(tensors, ritz);
    tensors.merge_multisite_tensor(multisite_tensor, status.chi_lim);

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

    if(status.iter <= settings::fdmrg::min_iters) { clear_convergence_status(); }

    if(status.variance_mpo_has_converged and status.entanglement_has_converged
       //           and status.variance_mpo_saturated_for >= min_saturation_iters
       //           and status.entanglement_saturated_for >= min_saturation_iters
    ) {
        tools::log->debug("Simulation has converged");
        status.algorithm_has_converged = true;
    }

    if(status.chi_lim_has_reached_chi_max and
       (status.variance_mpo_saturated_for >= max_saturation_iters or status.entanglement_saturated_for >= max_saturation_iters)) {
        tools::log->debug("Simulation has to stop");
        status.algorithm_has_to_stop = true;
    }

    if(tensors.position_is_any_edge() and status.variance_mpo_has_saturated and not status.variance_mpo_has_converged and not status.algorithm_has_converged and
       not has_projected) {
        tools::log->info("Projecting to {} due to saturation", settings::strategy::target_sector);
        *tensors.state = tools::finite::ops::get_projection_to_nearest_sector(*tensors.state, settings::strategy::target_sector);
        has_projected  = true;
    }

    tools::common::profile::t_con->toc();
}

bool class_fdmrg::algo_on() { return settings::fdmrg::on; }
long class_fdmrg::chi_lim_max() { return settings::fdmrg::chi_lim_max; }
// size_t class_fDMRG::write_freq() { return settings::fdmrg::write_freq; }
size_t class_fdmrg::print_freq() { return settings::fdmrg::print_freq; }
bool   class_fdmrg::chi_lim_grow() { return settings::fdmrg::chi_lim_grow; }
long   class_fdmrg::chi_lim_init() { return settings::fdmrg::chi_lim_init; }
bool   class_fdmrg::store_wave_function() { return settings::fdmrg::store_wavefn; }
