//
// Created by david on 2018-01-18.
//

#include "class_iDMRG.h"
#include <h5pp/h5pp.h>
#include <tools/infinite/opt.h>
#include <tools/common/prof.h>
#include <tools/common/log.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_infinite.h>
using namespace std;
using namespace Textra;

class_iDMRG::class_iDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
    : class_algorithm_infinite(std::move(h5ppFile_), SimulationType::iDMRG) {
    tools::log->trace("Constructing class {}", sim_name);

}



void class_iDMRG::run_simulation() {
    if(ritz == StateRitz::SR) state_name = "state_emin";
    else state_name = "state_emax";
    tools::log->info("Starting {} simulation of model [{}] for state [{}]", sim_name, enum2str(settings::model::model_type), state_name);
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
        if(sim_status.iter >= settings::idmrg::max_iters) {
            stop_reason = StopReason::MAX_ITERS;
            break;
        }
        if(sim_status.simulation_has_succeeded) {
            stop_reason = StopReason::SUCCEEDED;
            break;
        }
        if(sim_status.simulation_has_to_stop) {
            stop_reason = StopReason::SATURATED;
            break;
        }
        if(sim_status.num_resets > settings::precision::max_resets) {
            stop_reason = StopReason::MAX_RESET;
            break;
        }

        if (sim_status.iter >= settings::idmrg::max_iters)  {stop_reason = StopReason::MAX_ITERS; break;}
        if (sim_status.simulation_has_succeeded)                 {stop_reason = StopReason::SUCCEEDED; break;}
        if (sim_status.simulation_has_to_stop)                   {stop_reason = StopReason::SATURATED; break;}


        enlarge_environment();
        update_bond_dimension_limit();
        swap();
        sim_status.iter++;
        sim_status.step++;
    }
    tools::log->info("Finished {} simulation -- reason: {}",sim_name,enum2str(stop_reason));
}


void class_iDMRG::single_iDMRG_step(){
/*!
 * \fn void single_DMRG_step(class_superblock &state)
 */
    tools::log->trace("Starting single iDMRG step with ritz: [{}]", enum2str(ritz));
    tools::common::profile::t_sim->tic();
    Eigen::Tensor<Scalar,4> theta = tools::infinite::opt::find_ground_state(*state,ritz);
    tools::infinite::opt::truncate_theta(theta, *state);
    state->unset_measurements();
    tools::common::profile::t_sim->toc();
    sim_status.wall_time = tools::common::profile::t_tot->get_age();
    sim_status.simu_time = tools::common::profile::t_sim->get_measured_time();
}



void class_iDMRG::check_convergence(){
    tools::log->trace("Checking convergence");
    tools::common::profile::t_con->tic();
    check_convergence_entg_entropy();
    check_convergence_variance_mpo();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    if(sim_status.entanglement_has_converged and
       sim_status.variance_mpo_has_converged and
       sim_status.variance_ham_has_converged and
       sim_status.variance_mom_has_converged and
       sim_status.chi_lim_has_reached_chi_max)
    {
        sim_status.simulation_has_converged = true;
    }
    tools::common::profile::t_con->toc();
}




bool   class_iDMRG::sim_on()   {return settings::idmrg::on;}
long   class_iDMRG::chi_max()   {return settings::idmrg::chi_max;}
//size_t class_iDMRG::write_freq(){return settings::idmrg::write_freq;}
size_t class_iDMRG::print_freq(){return settings::idmrg::print_freq;}
bool   class_iDMRG::chi_grow()  {return settings::idmrg::chi_grow;}
long   class_iDMRG::chi_init()  {return settings::idmrg::chi_init;}



