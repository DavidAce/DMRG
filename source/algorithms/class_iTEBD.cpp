//
// Created by david on 2018-01-18.
//
#include "class_iTEBD.h"
#include <tools/common/prof.h>
#include <tools/common/log.h>
#include <tools/infinite/opt.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_infinite.h>
#include <state/class_mps_2site.h>
#include <model/class_model_base.h>
#include <general/nmspc_quantum_mechanics.h>
#include <h5pp/h5pp.h>
using namespace std;
using namespace Textra;

class_iTEBD::class_iTEBD(std::shared_ptr<h5pp::File> h5ppFile_)
        : class_algorithm_infinite(std::move(h5ppFile_),"iTEBD", SimulationType::iTEBD) {
    sim_status.delta_t      = settings::itebd::delta_t0;
    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    h_evn = state->HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    h_odd = state->HB->single_site_hamiltonian(1,2,SX,SY, SZ);
}



void class_iTEBD::run_preprocessing() {
    tools::common::profile::t_pre->tic();
    sim_status.delta_t = settings::itebd::delta_t0;
    unitary_time_evolving_operators = qm::timeEvolution::get_2site_evolution_gates(sim_status.delta_t, settings::itebd::suzuki_order, h_evn, h_odd);
    tools::common::profile::t_pre->toc();
}


void class_iTEBD::run_simulation()    {
    tools::log->info("Starting {} simulation", sim_name);
    while(sim_status.iter < settings::itebd::max_steps and not sim_status.simulation_has_converged) {
        single_TEBD_step();
        sim_status.phys_time += sim_status.delta_t;
        write_to_file();
//        write_measurements();
//        write_sim_status();
//        write_profiling();
        copy_from_tmp();
        print_status_update();
        check_convergence();
        sim_status.iter++;
    }
}


void class_iTEBD::run_postprocessing(){
    tools::common::profile::t_pos->tic();
    print_status_full();
    tools::common::profile::t_pos->toc();
    tools::common::profile::print_profiling();
}

void class_iTEBD::single_TEBD_step(){
/*!
 * \fn single_TEBD_step(class_superblock &state)
 * \brief infinite Time evolving block decimation.
 */
    tools::common::profile::t_sim->tic();
    for (auto &U: unitary_time_evolving_operators){
        Eigen::Tensor<Scalar,4> theta = tools::infinite::opt::time_evolve_theta(*state ,U);
        tools::infinite::opt::truncate_theta(theta, *state);
        if (&U != &unitary_time_evolving_operators.back()) {
            state->swap_AB();        }
    }
    state->unset_measurements();
    tools::common::profile::t_sim->toc();
    sim_status.wall_time = tools::common::profile::t_tot->get_age();
    sim_status.simu_time = tools::common::profile::t_sim->get_measured_time();
}



void class_iTEBD::check_convergence(){
    tools::common::profile::t_con->tic();
    check_convergence_entg_entropy();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    update_bond_dimension_limit();
    check_convergence_time_step();
    if(sim_status.entanglement_has_converged and
       sim_status.variance_ham_has_converged and
       sim_status.variance_mom_has_converged and
       sim_status.chi_lim_has_reached_chi_max and
       sim_status.time_step_has_converged)
    {
        sim_status.simulation_has_converged = true;
    }
    tools::common::profile::t_con->toc();
}

void class_iTEBD::check_convergence_time_step(){
    if(sim_status.delta_t <= settings::itebd::delta_tmin){
        sim_status.time_step_has_converged = true;
    }else if (sim_status.chi_lim_has_reached_chi_max and sim_status.entanglement_has_converged) {
        sim_status.delta_t = std::max(settings::itebd::delta_tmin, sim_status.delta_t * 0.5);
        unitary_time_evolving_operators = qm::timeEvolution::get_2site_evolution_gates(-sim_status.delta_t, settings::itebd::suzuki_order, h_evn, h_odd);
//        state->H->update_evolution_step_size(-sim_status.delta_t, settings::itebd::suzuki_order);
        clear_saturation_status();
    }
}

//void class_iTEBD::store_log_entry_progress(bool force){
//    if (not force){
//        if (math::mod(sim_status.iteration, settings::itebd::write_freq) != 0) {return;}
//    }
//    compute_observables();
//    t_sto->tic();
//    log_itebd->append_record(
//            sim_status.iteration,
//            state->measurements.bond_dimension.value(),
//            settings::itebd::chi_lim,
//            sim_status.delta_t,
//            state->measurements.energy_per_site.value(),
//            state->measurements.energy_per_site_ham.value(),
//            state->measurements.energy_per_site_mom.value(),
//            state->measurements.energy_variance_per_site.value(),
//            state->measurements.energy_variance_per_site_ham.value(),
//            state->measurements.energy_variance_per_site_mom.value(),
//            state->measurements.entanglement_entropy.value(),
//            state->measurements.truncation_error.value(),
//            sim_status.phys_time,
//            t_tot.get_age());
//
//    t_sto->toc();
//}




bool   class_iTEBD::sim_on()    {return settings::itebd::on;}
long   class_iTEBD::chi_max()   {return settings::itebd::chi_max;}
size_t class_iTEBD::write_freq(){return settings::itebd::write_freq;}
size_t class_iTEBD::print_freq(){return settings::itebd::print_freq;}
bool   class_iTEBD::chi_grow()  {return settings::itebd::chi_grow;}
long   class_iTEBD::chi_init()  {return settings::itebd::chi_init;}

