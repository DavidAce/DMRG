//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <sim_parameters/n_sim_settings.h>
#include <IO/class_hdf5_table_buffer.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <general/n_math.h>
#include "class_imaginary_TEBD.h"

using namespace std;
using namespace Textra;

class_imaginary_TEBD::class_imaginary_TEBD(std::shared_ptr<class_hdf5_file> hdf5_)
        : class_base() {
    simtype = SimulationType::iTEBD;
    simulation_name = "iTEBD";
    std::string group_name = simulation_name;
    std::string table_name = simulation_name;
    hdf5         = std::move(hdf5_);
    table_buffer = std::make_shared<class_hdf5_table_buffer>(hdf5, group_name, table_name);
    superblock   = std::make_shared<class_superblock>();
    observables  = std::make_shared<class_measurement>(superblock, simtype);
}


void class_imaginary_TEBD::run() {
    if (!settings::itebd::on) { return; }
    ccout(0) << "\nStarting " << simulation_name << " simulation" << std::endl;
    superblock->H.update_timestep(settings::itebd::delta_t, 1);
    t_tot.tic();
    while(iteration < max_steps) {
        store_table_entry();
        single_TEBD_step(chi_max);
        swap();
        iteration++;
        phys_time += superblock->H.timestep;
        print_status_update();
    }
    t_tot.toc();
    print_status_full();
    print_profiling();
}


void class_imaginary_TEBD::store_table_entry(){
    t_sto.tic();
    table_buffer->emplace_back(superblock->chi,
                               chi_max,
                               observables->get_energy(),
                               observables->get_entropy(),
                               observables->get_variance1(),
                               observables->get_variance2(),
                               observables->get_truncation_error(),
                               iteration,
                               superblock->chain_length,
                               superblock->H.timestep,
                               t_tot.get_age(),
                               phys_time);
    t_sto.toc();
}

void class_imaginary_TEBD::print_profiling(){
    if (settings::profiling::on) {
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_udt.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        t_evo.print_time_w_percent(t_sim);
        t_svd.print_time_w_percent(t_sim);
        t_mps.print_time_w_percent(t_sim);
    }
}



void class_imaginary_TEBD::print_status_update() {
    if (Math::mod(iteration, print_freq) != 0) {return;}
    if (print_freq == 0) {return;}

    t_obs.tic();
    std::cout << setprecision(16) << fixed << left;
    ccout(1) << left  << simulation_name;
    ccout(1) << left  << "Step: "                   << setw(10) << iteration;
    ccout(1) << left  << "E: "                      << setw(21) << setprecision(16)    << fixed   << observables->get_energy();
    ccout(1) << left  << "S: "                      << setw(21) << setprecision(16)    << fixed   << observables->get_entropy();
    ccout(1) << left  << "\u03C7_max: "             << setw(4)  << setprecision(3)     << fixed   << chi_max;
    ccout(1) << left  << "\u03C7: "                 << setw(4)  << setprecision(3)     << fixed   << observables->get_chi();
    ccout(1) << left  << "Var(E): "                 << setw(21) << setprecision(16)    << fixed   << observables->get_variance();
    ccout(1) << left  << "SVD truncation: "         << setw(21) << setprecision(16)    << fixed   << observables->get_truncation_error();
    ccout(1) << left  << "Timestep: "               << setw(25) << setprecision(12)    << fixed   << superblock->H.timestep;
    ccout(1) << std::endl;
    t_obs.toc();
}

void class_imaginary_TEBD::print_status_full(){
    std::cout << std::endl;
    std::cout << " -- Final results -- " << simulation_name << std::endl;
    ccout(0)  << setw(20) << "Energy               = " << setprecision(16) << fixed      << observables->get_energy()        << std::endl;
    ccout(0)  << setw(20) << "Entanglement Entropy = " << setprecision(16) << fixed      << observables->get_entropy()       << std::endl;
    ccout(0)  << setw(20) << "chi_max              = " << setprecision(4)  << fixed      << chi_max                          << std::endl;
    ccout(0)  << setw(20) << "chi                  = " << setprecision(4)  << fixed      << observables->get_chi()           << std::endl;
    ccout(0)  << setw(20) << "Variance             = " << setprecision(16) << scientific << observables->get_variance()      << std::endl;
    ccout(0)  << setw(20) << "SVD truncation       = " << setprecision(4)  << scientific << observables->get_truncation_error() << std::endl;
    ccout(0)  << setw(20) << "Timestep             = " << setprecision(12) << fixed      << superblock->H.timestep << std::endl;
    std::cout << std::endl;
}
