//
// Created by david on 2018-02-09.
//


#include <iomanip>
#include <sim_parameters/nmspc_sim_settings.h>
#include <IO/class_hdf5_table_buffer.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mps.h>
#include <mps_routines/class_environment_storage.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include "class_excited_state_DMRG.h"
using namespace std;
using namespace Textra;

class_excited_state_DMRG::class_excited_state_DMRG(std::shared_ptr<class_hdf5_file> hdf5_)
        : class_algorithm_base() {

    simtype = SimulationType::xDMRG;
    simulation_name = "xDMRG";
    group_name = simulation_name;
    table_name = simulation_name;
    hdf5         = std::move(hdf5_);
    table_buffer = std::make_shared<class_hdf5_table_buffer>(hdf5, group_name, table_name);
    superblock   = std::make_shared<class_superblock>();
    measurement  = std::make_shared<class_measurement>(superblock, simtype);
    env_storage  = std::make_shared<class_environment_storage>(max_length, superblock, hdf5);
}



void class_excited_state_DMRG::run() {
    if (!settings::xdmrg::on) { return; }
    ccout(0) << "\nStarting " << simulation_name << " simulation" << std::endl;
    t_tot.tic();
    initialize_random_chain();
    std::cout << "starting sweeps" << std::endl;
    while(sweep < max_sweeps) {
        single_DMRG_step(chi_max);
        env_storage_overwrite_MPS();         //Needs to occurr after update_MPS...
        store_table_entry();
        print_status_update();

        position = enlarge_environment(direction);
        position = env_storage_move();
    }
    t_tot.toc();
    print_status_full();
    print_profiling();

}

int class_excited_state_DMRG::initialize_random_chain() {

    while (true) {
        auto r1 = rn::uniform_complex_1();
        auto r2 = rn::uniform_complex_1();
        superblock->MPS->GA(1, 0, 0) = r1.real();
        superblock->MPS->GA(0, 0, 0) = r1.imag();
        superblock->MPS->GB(1, 0, 0) = r2.real();
        superblock->MPS->GB(0, 0, 0) = r2.imag();
        position = env_storage_insert();
        if (superblock->chain_length < max_length) {
            enlarge_environment();
            position++;
        } else {
            break;
        }

    }
    return position;
}



int class_excited_state_DMRG::env_storage_insert() {
    t_ste.tic();
    int position = env_storage->insert();
    t_ste.toc();
    return position;
}

int class_excited_state_DMRG::env_storage_load(){
    t_ste.tic();
    int position = env_storage->load();
    t_ste.toc();
    return position;
}
void class_excited_state_DMRG::env_storage_overwrite_MPS(){
    t_ste.tic();
    env_storage->overwrite_MPS();
    t_ste.toc();
}
int class_excited_state_DMRG::env_storage_move(){
    t_ste.tic();
    int position = env_storage->move(direction, sweep);
    t_ste.toc();
    return position;
}


void class_excited_state_DMRG::store_table_entry(){
    t_sto.tic();
    table_buffer->emplace_back(superblock->chi,
                               chi_max,
                               measurement->get_energy(),
                               measurement->get_entropy(),
                               measurement->get_variance1(),
                               measurement->get_variance2(),
                               measurement->get_truncation_error(),
                               position,
                               superblock->chain_length,
                               0,
                               t_tot.get_age(),
                               0);
    t_sto.toc();
}

void class_excited_state_DMRG::print_profiling(){
    if (settings::profiling::on) {
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_ste.print_time_w_percent(t_tot);
        t_env.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        t_eig.print_time_w_percent(t_sim);
        t_svd.print_time_w_percent(t_sim);
        t_mps.print_time_w_percent(t_sim);
    }
}

void class_excited_state_DMRG::print_status_update() {
    if (position != middle_of_chain) {return;}
    if (print_freq == 0) {return;}

    t_obs.tic();
    std::cout << setprecision(16) << fixed << left;
    ccout(1) << left  << simulation_name << " ";
    ccout(1) << left  << "Pos: "                    << setw(6) << position;
    ccout(1) << left  << "Dir: "                    << setw(3) << direction;
    ccout(1) << left  << "Sweep: "                  << setw(4) << sweep;
    ccout(1) << left  << "E: "                      << setw(21) << setprecision(16)    << fixed   << measurement->get_energy();
    ccout(1) << left  << "S: "                      << setw(21) << setprecision(16)    << fixed   << measurement->get_entropy();
    ccout(1) << left  << "\u03C7_max: "             << setw(4)  << setprecision(3)     << fixed   << chi_max;
    ccout(1) << left  << "\u03C7: "                 << setw(4)  << setprecision(3)     << fixed   << measurement->get_chi();
    ccout(1) << left  << "Var(E): "                 << setw(21) << setprecision(16)    << fixed   << measurement->get_variance();
    ccout(1) << left  << "SVD truncation: "         << setw(21) << setprecision(16)    << fixed   << measurement->get_truncation_error();
    ccout(1) << left  << "Chain length: "           << setw(25) << setprecision(12)    << fixed   << superblock->chain_length;
    ccout(1) << std::endl;
    t_obs.toc();

}

void class_excited_state_DMRG::print_status_full(){
    std::cout << std::endl;
    std::cout << " -- Final results -- " << simulation_name << std::endl;
    ccout(0)  << setw(20) << "Energy               = " << setprecision(16) << fixed      << measurement->get_energy()        << std::endl;
    ccout(0)  << setw(20) << "Entanglement Entropy = " << setprecision(16) << fixed      << measurement->get_entropy()       << std::endl;
    ccout(0)  << setw(20) << "chi_max              = " << setprecision(4)  << fixed      << chi_max                          << std::endl;
    ccout(0)  << setw(20) << "chi                  = " << setprecision(4)  << fixed      << measurement->get_chi()           << std::endl;
    ccout(0)  << setw(20) << "Variance             = " << setprecision(16) << scientific << measurement->get_variance()      << std::endl;
    ccout(0)  << setw(20) << "SVD truncation       = " << setprecision(4)  << scientific << measurement->get_truncation_error() << std::endl;
    ccout(0) << setw(20)  << "Chain length         = " << setprecision(12) << fixed      << superblock->chain_length << std::endl;
    std::cout << std::endl;
}