//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <io/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_state/class_superblock.h>
#include <mps_state/class_mps_2site.h>
#include <mps_tools/nmspc_mps_tools.h>
#include <model/class_hamiltonian_base.h>
#include <general/nmspc_math.h>
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
    sim_state.delta_t      = settings::itebd::delta_t0;
//    initialize_superblock(settings::model::initial_state);
    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    h_evn = superblock->HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    h_odd = superblock->HB->single_site_hamiltonian(1,2,SX,SY, SZ);
}



void class_iTEBD::run_preprocessing() {
    t_tot.tic();
    sim_state.delta_t = settings::itebd::delta_t0;
    unitary_time_evolving_operators = qm::timeEvolution::get_2site_evolution_gates(sim_state.delta_t, settings::itebd::suzuki_order, h_evn, h_odd);
    t_tot.toc();
}


void class_iTEBD::run_simulation()    {
    log->info("Starting {} simulation", sim_name);
    t_tot.tic();
    while(sim_state.iteration < settings::itebd::max_steps and not sim_state.simulation_has_converged) {
        single_TEBD_step(sim_state.chi_temp);
        sim_state.phys_time += sim_state.delta_t;
        store_table_entry_progress();
        store_profiling_deltas();
        print_status_update();
        check_convergence();
        sim_state.iteration++;
    }
    t_tot.toc();



}
void class_iTEBD::run_postprocessing(){
    print_status_full();
    print_profiling();
}

void class_iTEBD::single_TEBD_step(long chi){
/*!
 * \fn single_iTEBD_step(class_superblock &superblock)
 * \brief infinite Time evolving block decimation.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 */
    t_sim.tic();
    for (auto &U: unitary_time_evolving_operators){
        t_evo.tic();
        Eigen::Tensor<Scalar,4> theta = superblock->evolve_MPS(superblock->MPS->get_theta() ,U);
        t_evo.toc();

        t_svd.tic();
        superblock->truncate_MPS(theta, chi, settings::precision::SVDThreshold);
        t_svd.toc();

        if (&U != &unitary_time_evolving_operators.back()) {
            superblock->swap_AB();        }
    }
    superblock->unset_measurements();
    t_sim.toc();
    sim_state.wall_time = t_tot.get_age();
    sim_state.simu_time = t_sim.get_age();
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
    if(sim_state.delta_t <= settings::itebd::delta_tmin){
        sim_state.time_step_has_converged = true;
    }else if (sim_state.bond_dimension_has_reached_max and sim_state.entanglement_has_converged) {
        sim_state.delta_t = std::max(settings::itebd::delta_tmin, sim_state.delta_t * 0.5);
        unitary_time_evolving_operators = qm::timeEvolution::get_2site_evolution_gates(-sim_state.delta_t, settings::itebd::suzuki_order, h_evn, h_odd);
//        superblock->H->update_evolution_step_size(-sim_state.delta_t, settings::itebd::suzuki_order);
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
    if(sim_state.entanglement_has_converged and
       sim_state.variance_ham_has_converged and
       sim_state.variance_mom_has_converged and
       sim_state.bond_dimension_has_reached_max and
       sim_state.time_step_has_converged)
    {
        sim_state.simulation_has_converged = true;
    }
    t_con.toc();
}


//void class_iTEBD::store_table_entry_progress(bool force){
//    if (not force){
//        if (Math::mod(sim_state.iteration, settings::itebd::store_freq) != 0) {return;}
//    }
//    compute_observables();
//    t_sto.tic();
//    table_itebd->append_record(
//            sim_state.iteration,
//            superblock->measurements.bond_dimension.value(),
//            settings::itebd::chi_max,
//            sim_state.delta_t,
//            superblock->measurements.energy_per_site.value(),
//            superblock->measurements.energy_per_site_ham.value(),
//            superblock->measurements.energy_per_site_mom.value(),
//            superblock->measurements.energy_variance_per_site.value(),
//            superblock->measurements.energy_variance_per_site_ham.value(),
//            superblock->measurements.energy_variance_per_site_mom.value(),
//            superblock->measurements.current_entanglement_entropy.value(),
//            superblock->measurements.truncation_error.value(),
//            sim_state.phys_time,
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

