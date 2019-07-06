//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <io/class_hdf5_table_buffer2.h>
#include <simulation/nmspc_settings.h>
#include <state/class_infinite_state.h>
#include <state/class_mps_2site.h>
#include <state/tools/nmspc_tools.h>
#include <model/class_model_base.h>
#include <math/nmspc_math.h>
#include <general/nmspc_quantum_mechanics.h>
#include <spdlog/spdlog.h>
#include <h5pp/h5pp.h>
#include "class_iTEBD.h"
using namespace std;
using namespace Textra;
//using namespace std::complex_literals;

class_iTEBD::class_iTEBD(std::shared_ptr<h5pp::File> h5ppFile_)
        : class_algorithm_infinite(std::move(h5ppFile_),"iTEBD", SimulationType::iTEBD) {
//    initialize_constants();
//    table_itebd = std::make_unique<class_hdf5_table<class_table_tebd>>(h5pp_file, sim_name + "/measurements", "simulation_progress", sim_name);
    sim_status.delta_t      = settings::itebd::delta_t0;
//    initialize_superblock(settings::model::initial_state);
    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    h_evn = state->HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    h_odd = state->HB->single_site_hamiltonian(1,2,SX,SY, SZ);
}



void class_iTEBD::run_preprocessing() {
    t_tot.tic();
    sim_status.delta_t = settings::itebd::delta_t0;
    unitary_time_evolving_operators = qm::timeEvolution::get_2site_evolution_gates(sim_status.delta_t, settings::itebd::suzuki_order, h_evn, h_odd);
    t_tot.toc();
}


void class_iTEBD::run_simulation()    {
    log->info("Starting {} simulation", sim_name);
    t_tot.tic();
    while(sim_status.iteration < settings::itebd::max_steps and not sim_status.simulation_has_converged) {
        single_TEBD_step(sim_status.chi_temp);
        sim_status.phys_time += sim_status.delta_t;
//        store_table_entry_progress();
        store_profiling_deltas();
        print_status_update();
        check_convergence();
        sim_status.iteration++;
    }
    t_tot.toc();
}


void class_iTEBD::run_postprocessing(){
    print_status_full();
    print_profiling();
}

void class_iTEBD::single_TEBD_step(long chi){
/*!
 * \fn single_iTEBD_step(class_superblock &state)
 * \brief infinite Time evolving block decimation.
 * \param state A class containing MPS, environment and Hamiltonian MPO objects.
 */
    t_sim.tic();
    for (auto &U: unitary_time_evolving_operators){
        t_evo.tic();
        Eigen::Tensor<Scalar,4> theta = tools::infinite::opt::time_evolve_theta(*state ,U);
        t_evo.toc();

        t_svd.tic();
        tools::infinite::opt::truncate_theta(theta, *state, chi, settings::precision::SVDThreshold);
        t_svd.toc();

        if (&U != &unitary_time_evolving_operators.back()) {
            state->swap_AB();        }
    }
    state->unset_measurements();
    t_sim.toc();
    sim_status.wall_time = t_tot.get_age();
    sim_status.simu_time = t_sim.get_age();
}


//void class_iTEBD::initialize_constants(){
//    using namespace settings;
//    max_steps    = itebd::max_steps   ;
//    settings::itebd::delta_t0     = itebd::settings::itebd::delta_t0    ;
//    settings::itebd::delta_tmin   = itebd::settings::itebd::delta_tmin  ;
//    settings::itebd::suzuki_order = itebd::settings::itebd::suzuki_order;
//    settings::itebd::chi_max      = itebd::settings::itebd::chi_max     ;
//    chi_grow     = itebd::chi_grow    ;
//    print_freq   = itebd::print_freq  ;
//    settings::itebd::store_freq   = itebd::settings::itebd::store_freq  ;
//}



void class_iTEBD::check_convergence_time_step(){
    if(sim_status.delta_t <= settings::itebd::delta_tmin){
        sim_status.time_step_has_converged = true;
    }else if (sim_status.bond_dimension_has_reached_max and sim_status.entanglement_has_converged) {
        sim_status.delta_t = std::max(settings::itebd::delta_tmin, sim_status.delta_t * 0.5);
        unitary_time_evolving_operators = qm::timeEvolution::get_2site_evolution_gates(-sim_status.delta_t, settings::itebd::suzuki_order, h_evn, h_odd);
//        state->H->update_evolution_step_size(-sim_status.delta_t, settings::itebd::suzuki_order);
        clear_saturation_status();
    }
}

void class_iTEBD::check_convergence(){
    t_con.tic();
    check_convergence_entg_entropy();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    update_bond_dimension();
    check_convergence_time_step();
    if(sim_status.entanglement_has_converged and
       sim_status.variance_ham_has_converged and
       sim_status.variance_mom_has_converged and
       sim_status.bond_dimension_has_reached_max and
       sim_status.time_step_has_converged)
    {
        sim_status.simulation_has_converged = true;
    }
    t_con.toc();
}


//void class_iTEBD::store_table_entry_progress(bool force){
//    if (not force){
//        if (math::mod(sim_status.iteration, settings::itebd::store_freq) != 0) {return;}
//    }
//    compute_observables();
//    t_sto.tic();
//    table_itebd->append_record(
//            sim_status.iteration,
//            state->measurements.bond_dimension.value(),
//            settings::itebd::chi_max,
//            sim_status.delta_t,
//            state->measurements.energy_per_site.value(),
//            state->measurements.energy_per_site_ham.value(),
//            state->measurements.energy_per_site_mom.value(),
//            state->measurements.energy_variance_per_site.value(),
//            state->measurements.energy_variance_per_site_ham.value(),
//            state->measurements.energy_variance_per_site_mom.value(),
//            state->measurements.current_entanglement_entropy.value(),
//            state->measurements.truncation_error.value(),
//            sim_status.phys_time,
//            t_tot.get_age());
//
//    t_sto.toc();
//}

bool   class_iTEBD::sim_on()    {return settings::itebd::on;}
long   class_iTEBD::chi_max()   {return settings::itebd::chi_max;}
size_t class_iTEBD::num_sites() {return 2u;}
size_t class_iTEBD::store_freq(){return settings::itebd::store_freq;}
size_t class_iTEBD::print_freq(){return settings::itebd::print_freq;}
bool   class_iTEBD::chi_grow()  {return settings::itebd::chi_grow;}

