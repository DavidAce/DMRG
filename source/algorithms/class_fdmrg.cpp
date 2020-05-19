//
// Created by david on 2018-01-31.
//

#include "class_fdmrg.h"
#include <config/nmspc_settings.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/ops.h>
#include <tools/finite/opt.h>
#include <tools/finite/svd.h>

class_fdmrg::class_fdmrg(std::shared_ptr<h5pp::File> h5pp_file_) : class_algorithm_finite(std::move(h5pp_file_), AlgorithmType::fDMRG) {
    tools::log->trace("Constructing class {}", algo_name);
}

void class_fdmrg::run_simulation() {
    if(ritz == StateRitz::SR) state_name = "state_emin";
    else state_name = "state_emax";
    tools::log->info("Starting {} simulation of model [{}] for state [{}]", algo_name, enum2str(settings::model::model_type), state_name);
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
            if(status.num_resets > settings::precision::max_resets) {
                stop_reason = StopReason::MAX_RESET;
                break;
            }
        }
        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);
        move_center_point();
    }
    tools::log->info("Finished {} simulation -- reason: {}", algo_name, enum2str(stop_reason));
}

void class_fdmrg::single_fDMRG_step() {
    /*!
     * \fn void single_DMRG_step(std::string ritz)
     */
    tools::log->trace("Starting single fDMRG step with ritz: [{}]", enum2str(ritz));
    tools::common::profile::t_sim->tic();
    Eigen::Tensor<Scalar, 3> multisite_tensor = tools::finite::opt::find_ground_state(tensors, ritz);
    tensors.merge_multisite_tensor(multisite_tensor);
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

    if(status.iter <= settings::fdmrg::min_iters) {
        clear_saturation_status();
    }

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

    if(tensors.position_is_any_edge() and status.variance_mpo_has_saturated and not status.variance_mpo_has_converged and
       not status.algorithm_has_converged and not has_projected) {
        tools::log->info("Projecting to {} due to saturation", settings::strategy::target_parity_sector);
        *tensors.state        = tools::finite::ops::get_projection_to_closest_parity_sector(*tensors.state, settings::strategy::target_parity_sector);
        has_projected = true;
    }

    tools::common::profile::t_con->toc();
}

bool   class_fdmrg::algo_on() { return settings::fdmrg::on; }
long   class_fdmrg::chi_max() { return settings::fdmrg::chi_max; }
//size_t class_fDMRG::write_freq() { return settings::fdmrg::write_freq; }
size_t class_fdmrg::print_freq() { return settings::fdmrg::print_freq; }
bool   class_fdmrg::chi_grow() { return settings::fdmrg::chi_grow; }
long   class_fdmrg::chi_init() { return settings::fdmrg::chi_init; }
bool   class_fdmrg::store_wave_function() { return settings::fdmrg::store_wavefn; }
