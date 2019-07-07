//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <h5pp/h5pp.h>
#include <io/class_hdf5_table_buffer2.h>
#include <simulation/nmspc_settings.h>
#include <state/class_infinite_state.h>
#include <state/tools/nmspc_tools.h>
#include <math/nmspc_math.h>
#include <spdlog/spdlog.h>
#include "class_iDMRG.h"
using namespace std;
using namespace Textra;

class_iDMRG::class_iDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
    : class_algorithm_infinite(std::move(h5ppFile_),"iDMRG", SimulationType::iDMRG) {

//    initialize_superblock(settings::model::initial_state);

}



void class_iDMRG::run_simulation() {
    if (not settings::idmrg::on) { return; }
    log->info("Starting {} simulation", sim_name);
    t_tot.tic();
    while(true){
        single_DMRG_step("SR");
        print_status_update();
//        store_table_entry_progress();
        store_profiling_deltas();
        check_convergence();

        // It's important not to perform the last swap.
        // That last state would not get optimized

        if (sim_status.iteration >= settings::idmrg::max_steps)  {stop_reason = StopReason::MAX_STEPS; break;}
        if (sim_status.simulation_has_converged)                 {stop_reason = StopReason::CONVERGED; break;}
        if (sim_status.simulation_has_to_stop)                   {stop_reason = StopReason::SATURATED; break;}


        enlarge_environment();
        swap();
        sim_status.iteration++;
    }
    t_tot.toc();
    switch(stop_reason){
        case StopReason::MAX_STEPS : log->info("Finished {} simulation -- reason: MAX_STEPS",sim_name) ;break;
        case StopReason::CONVERGED : log->info("Finished {} simulation -- reason: CONVERGED",sim_name) ;break;
        case StopReason::SATURATED : log->info("Finished {} simulation -- reason: SATURATED",sim_name) ;break;
        default: log->info("Finished {} simulation -- reason: NONE GIVEN",sim_name);
    }


}


void class_iDMRG::single_DMRG_step(std::string ritz){
/*!
 * \fn void single_DMRG_step(class_superblock &state)
 */
    log->trace("Starting infinite DMRG step");
    t_sim.tic();
    Eigen::Tensor<Scalar,4> theta = tools::infinite::opt::find_ground_state(*state,ritz);
    tools::infinite::opt::truncate_theta(theta, *state, sim_status.chi_temp, settings::precision::SVDThreshold);
    state->unset_measurements();
    t_sim.toc();
    sim_status.wall_time = t_tot.get_age();
    sim_status.simu_time = t_sim.get_age();
}



void class_iDMRG::check_convergence(){
    log->trace("Checking convergence");
    t_con.tic();
    check_convergence_entg_entropy();
    check_convergence_variance_mpo();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    update_bond_dimension();
    if(sim_status.entanglement_has_converged and
       sim_status.variance_mpo_has_converged and
       sim_status.variance_ham_has_converged and
       sim_status.variance_mom_has_converged and
       sim_status.bond_dimension_has_reached_max)
    {
        sim_status.simulation_has_converged = true;
    }
    t_con.toc();
}






bool   class_iDMRG::sim_on()   {return settings::idmrg::on;}
long   class_iDMRG::chi_max()   {return settings::idmrg::chi_max;}
size_t class_iDMRG::num_sites() {return 2u;}
size_t class_iDMRG::store_freq(){return settings::idmrg::store_freq;}
size_t class_iDMRG::print_freq(){return settings::idmrg::print_freq;}
bool   class_iDMRG::chi_grow()  {return settings::idmrg::chi_grow;}



