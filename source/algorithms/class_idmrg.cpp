//
// Created by david on 2018-01-18.
//

#include "class_idmrg.h"
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <tensors/class_tensors_infinite.h>
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/infinite/opt.h>
using namespace std;
using namespace Textra;

class_idmrg::class_idmrg(std::shared_ptr<h5pp::File> h5ppFile_) : class_algorithm_infinite(std::move(h5ppFile_), AlgorithmType::iDMRG) {
    tools::log->trace("Constructing class_idmrg");
    tensors.initialize(settings::model::model_type);
}

void class_idmrg::run_simulation() {
    if(ritz == StateRitz::SR)
        tensors.state->set_name("state_emin");
    else
        tensors.state->set_name("state_emax");
    tools::log->info("Starting {} simulation of model [{}] for state [{}]", algo_name, enum2str(settings::model::model_type), tensors.state->get_name());
    auto t_sim = tools::common::profile::prof[algo_type]["t_sim"]->tic_token();
    while(true) {
        single_iDMRG_step();
        print_status_update();
        write_to_file();
        copy_from_tmp();
        check_convergence();
        update_bond_dimension_limit(); // Will update bond dimension if the state precision is being limited by bond dimension

        // It's important not to perform the last move.
        // That last state would not get optimized
        if(status.iter >= settings::idmrg::max_iters) {
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

        if(status.iter >= settings::idmrg::max_iters) {
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

        update_bond_dimension_limit();
        tensors.enlarge();
        status.iter++;
        status.step++;
    }
    tools::log->info("Finished {} simulation -- reason: {}", algo_name, enum2str(stop_reason));
}

void class_idmrg::single_iDMRG_step() {
    /*!
     * \fn void single_DMRG_step(class_superblock &state)
     */
    tools::log->trace("Starting single iDMRG step with ritz: [{}]", enum2str(ritz));
    Eigen::Tensor<Scalar, 3> twosite_tensor = tools::infinite::opt::find_ground_state(tensors, ritz);
    tensors.merge_twosite_tensor(twosite_tensor, status.chi_lim);
    status.wall_time = tools::common::profile::t_tot->get_measured_time();
    status.algo_time = tools::common::profile::prof[algo_type]["t_sim"]->get_measured_time();
}

void class_idmrg::check_convergence() {
    tools::log->trace("Checking convergence");
    auto t_con = tools::common::profile::prof[algo_type]["t_con"]->tic_token();
    check_convergence_entg_entropy();
    check_convergence_variance_mpo();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    if(status.entanglement_has_converged and status.variance_mpo_has_converged and status.variance_ham_has_converged and status.variance_mom_has_converged and
       status.chi_lim_has_reached_chi_max) {
        status.algorithm_has_converged = true;
    }
}

bool   class_idmrg::cfg_algorithm_is_on() { return settings::idmrg::on; }
long   class_idmrg::cfg_chi_lim_max() { return settings::idmrg::chi_lim_max; }
size_t class_idmrg::cfg_print_freq() { return settings::idmrg::print_freq; }
bool   class_idmrg::cfg_chi_lim_grow() { return settings::idmrg::chi_lim_grow; }
long   class_idmrg::cfg_chi_lim_init() { return settings::idmrg::chi_lim_init; }
