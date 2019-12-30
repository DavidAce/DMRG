//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <h5pp/h5pp.h>
#include <io/class_h5table_buffer.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_infinite.h>
#include <tools/nmspc_tools.h>
#include <math/nmspc_math.h>
#include "class_iDMRG.h"
using namespace std;
using namespace Textra;

class_iDMRG::class_iDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
    : class_algorithm_infinite(std::move(h5ppFile_),"iDMRG", SimulationType::iDMRG) {
    log->trace("Constructing class_iDMRG");

}



void class_iDMRG::run_simulation() {
    if (not settings::idmrg::on) { return; }
    log->info("Starting {} simulation", sim_name);
    while(true){
        single_DMRG_step("SR");
        write_state();
        write_measurements();
        write_sim_status();
        write_profiling();
        copy_from_tmp();
        print_status_update();
        check_convergence();

        // It's important not to perform the last swap.
        // That last state would not get optimized

        if (sim_status.iteration >= settings::idmrg::max_steps)  {stop_reason = StopReason::MAX_ITERS; break;}
        if (sim_status.simulation_has_succeeded)                 {stop_reason = StopReason::SUCCEEDED; break;}
        if (sim_status.simulation_has_to_stop)                   {stop_reason = StopReason::SATURATED; break;}


        enlarge_environment();
        update_bond_dimension_limit();
        swap();
        sim_status.iteration++;
        sim_status.step++;
    }
    switch(stop_reason){
        case StopReason::MAX_ITERS : log->info("Finished {} simulation -- reason: MAX ITERS",sim_name) ;break;
        case StopReason::SUCCEEDED : log->info("Finished {} simulation -- reason: SUCCEEDED",sim_name) ;break;
        case StopReason::SATURATED : log->info("Finished {} simulation -- reason: SATURATED",sim_name) ;break;
        default: log->info("Finished {} simulation -- reason: NONE GIVEN",sim_name);
    }

}


void class_iDMRG::single_DMRG_step(std::string ritz){
/*!
 * \fn void single_DMRG_step(class_superblock &state)
 */
    log->trace("Starting infinite DMRG step");
    tools::common::profile::t_sim.tic();
    Eigen::Tensor<Scalar,4> theta = tools::infinite::opt::find_ground_state(*state,ritz);
    tools::infinite::opt::truncate_theta(theta, *state);
    state->unset_measurements();
    tools::common::profile::t_sim.toc();
    sim_status.wall_time = tools::common::profile::t_tot.get_age();
    sim_status.simu_time = tools::common::profile::t_sim.get_measured_time();
}



void class_iDMRG::check_convergence(){
    log->trace("Checking convergence");
    tools::common::profile::t_con.tic();
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
    tools::common::profile::t_con.toc();
}




bool   class_iDMRG::sim_on()   {return settings::idmrg::on;}
long   class_iDMRG::chi_max()   {return settings::idmrg::chi_max;}
size_t class_iDMRG::num_sites() {return 2u;}
size_t class_iDMRG::write_freq(){return settings::idmrg::write_freq;}
size_t class_iDMRG::print_freq(){return settings::idmrg::print_freq;}
bool   class_iDMRG::chi_grow()  {return settings::idmrg::chi_grow;}
long   class_iDMRG::chi_init()  {return settings::idmrg::chi_init;}



