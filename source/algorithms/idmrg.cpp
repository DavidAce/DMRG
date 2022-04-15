#include "idmrg.h"
#include "config/settings.h"
#include "math/tenx.h"
#include "tensors/edges/EdgesInfinite.h"
#include "tensors/model/ModelInfinite.h"
#include "tensors/state/StateInfinite.h"
#include "tensors/TensorsInfinite.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/common/prof.h"
#include "tools/infinite/opt.h"

idmrg::idmrg(std::shared_ptr<h5pp::File> h5ppFile_) : AlgorithmInfinite(std::move(h5ppFile_), AlgorithmType::iDMRG) {
    tools::log->trace("Constructing class_idmrg");
    tensors.initialize(settings::model::model_type);
}

void idmrg::run_simulation() {
    if(ritz == OptRitz::SR)
        tensors.state->set_name("state_emin");
    else
        tensors.state->set_name("state_emax");
    tools::log->info("Starting {} simulation of model [{}] for state [{}]", status.algo_type_sv(), enum2sv(settings::model::model_type),
                     tensors.state->get_name());
    auto t_algo = tid::tic_scope(status.algo_type_sv());
    while(true) {
        single_iDMRG_step();
        check_convergence();
        print_status_update();
        write_to_file();

        // It's important not to perform the last move.
        // That last state would not get optimized
        if(status.iter >= settings::idmrg::max_iters) {
            status.algo_stop = AlgorithmStop::MAX_ITERS;
            break;
        }
        if(status.algorithm_has_succeeded) {
            status.algo_stop = AlgorithmStop::SUCCESS;
            break;
        }
        if(status.algorithm_has_to_stop) {
            status.algo_stop = AlgorithmStop::SATURATED;
            break;
        }
        if(status.num_resets > settings::strategy::max_resets) {
            status.algo_stop = AlgorithmStop::MAX_RESET;
            break;
        }

        if(status.iter >= settings::idmrg::max_iters) {
            status.algo_stop = AlgorithmStop::MAX_ITERS;
            break;
        }
        if(status.algorithm_has_succeeded) {
            status.algo_stop = AlgorithmStop::SUCCESS;
            break;
        }
        if(status.algorithm_has_to_stop) {
            status.algo_stop = AlgorithmStop::SATURATED;
            break;
        }

        update_bond_dimension_limit(); // Will update bond dimension if the state precision is being limited by bond dimension
        tensors.enlarge();
        status.iter++;
        status.step++;
        status.wall_time = tid::get_unscoped("t_tot").get_time();
        status.algo_time = t_algo->get_time();
    }
    tools::log->info("Finished {} simulation -- reason: {}", status.algo_type_sv(), status.algo_stop_sv());
}

void idmrg::single_iDMRG_step() {
    /*!
     * \fn void single_DMRG_step(class_superblock &state)
     */
    tools::log->trace("Starting single iDMRG step with ritz: [{}]", enum2sv(ritz));
    Eigen::Tensor<Scalar, 3> twosite_tensor = tools::infinite::opt::find_ground_state(tensors, ritz);
    tensors.merge_twosite_tensor(twosite_tensor, status.bond_lim);
}

void idmrg::check_convergence() {
    tools::log->trace("Checking convergence");
    auto t_con = tid::tic_scope("conv");
    check_convergence_entg_entropy();
    check_convergence_variance_mpo();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    if(status.entanglement_converged_for > 0 and status.variance_mpo_converged_for > 0 and status.variance_ham_converged_for > 0 and
       status.variance_mom_converged_for > 0 and status.bond_limit_has_reached_max) {
        status.algorithm_converged_for++;
    } else
        status.algorithm_converged_for = 0;
}
