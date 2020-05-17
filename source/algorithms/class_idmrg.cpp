//
// Created by david on 2018-01-18.
//

#include "class_idmrg.h"
#include <h5pp/h5pp.h>
#include <config/nmspc_settings.h>
#include <tensors/class_tensors_infinite.h>
#include <tensors/state/class_state_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/edges/class_edges_infinite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/infinite/opt.h>

using namespace std;
using namespace Textra;

class_idmrg::class_idmrg(std::shared_ptr<h5pp::File> h5ppFile_)
    : class_algorithm_infinite(std::move(h5ppFile_), AlgorithmType::iDMRG) {
    tools::log->trace("Constructing class {}", algo_name);
    tensors.initialize(settings::model::model_type);
}



void class_idmrg::run_simulation() {
    if(ritz == StateRitz::SR) state_name = "state_emin";
    else state_name = "state_emax";
    tools::log->info("Starting {} simulation of model [{}] for state [{}]", algo_name, enum2str(settings::model::model_type), state_name);
    while(true){
        single_iDMRG_step();
        print_status_update();
        write_to_file();
        copy_from_tmp();
        check_convergence();
        update_truncation_limit();     // Will update SVD threshold iff the state precision is being limited by truncation error
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
        if(status.num_resets > settings::precision::max_resets) {
            stop_reason = StopReason::MAX_RESET;
            break;
        }

        if (status.iter >= settings::idmrg::max_iters)  {stop_reason = StopReason::MAX_ITERS; break;}
        if (status.algorithm_has_succeeded)             {stop_reason = StopReason::SUCCEEDED; break;}
        if (status.algorithm_has_to_stop)               {stop_reason = StopReason::SATURATED; break;}

        enlarge_environment();
        update_bond_dimension_limit();
        swap();
        status.iter++;
        status.step++;
    }
    tools::log->info("Finished {} simulation -- reason: {}", algo_name,enum2str(stop_reason));
}


void class_idmrg::single_iDMRG_step(){
/*!
 * \fn void single_DMRG_step(class_superblock &state)
 */
    tools::log->trace("Starting single iDMRG step with ritz: [{}]", enum2str(ritz));
    tools::common::profile::t_sim->tic();
    Eigen::Tensor<Scalar,3> theta = tools::infinite::opt::find_ground_state(tensors,ritz);
    tools::infinite::opt::truncate_theta(theta, *tensors.state);
    tensors.state->clear_measurements();
    tools::common::profile::t_sim->toc();
    status.wall_time = tools::common::profile::t_tot->get_age();
    status.simu_time = tools::common::profile::t_sim->get_measured_time();
}



void class_idmrg::check_convergence(){
    tools::log->trace("Checking convergence");
    tools::common::profile::t_con->tic();
    check_convergence_entg_entropy();
    check_convergence_variance_mpo();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    if(status.entanglement_has_converged and status.variance_mpo_has_converged and status.variance_ham_has_converged and status.variance_mom_has_converged and
       status.chi_lim_has_reached_chi_max)
    {
        status.algorithm_has_converged = true;
    }
    tools::common::profile::t_con->toc();
}




bool   class_idmrg::algo_on()   {return settings::idmrg::on;}
long   class_idmrg::chi_max()   {return settings::idmrg::chi_max;}
size_t class_idmrg::print_freq(){return settings::idmrg::print_freq;}
bool   class_idmrg::chi_grow()  {return settings::idmrg::chi_grow;}
long   class_idmrg::chi_init()  {return settings::idmrg::chi_init;}



