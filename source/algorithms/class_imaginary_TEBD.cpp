//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <sim_parameters/nmspc_sim_settings.h>
#include <IO/class_hdf5_table_buffer.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mpo.h>
#include <general/nmspc_math.h>
#include "class_imaginary_TEBD.h"

using namespace std;
using namespace Textra;

class_imaginary_TEBD::class_imaginary_TEBD(std::shared_ptr<class_hdf5_file> hdf5_)
        : class_algorithm_base() {
    simtype = SimulationType::iTEBD;
    simulation_name = "iTEBD";
    std::string group_name = simulation_name;
    std::string table_name = simulation_name;
    hdf5         = std::move(hdf5_);
    table_buffer = std::make_shared<class_hdf5_table_buffer>(hdf5, group_name, table_name);
    superblock   = std::make_shared<class_superblock>();
    measurement  = std::make_shared<class_measurement>(superblock, simtype);
}


void class_imaginary_TEBD::run() {
    if (!settings::itebd::on) { return; }
    ccout(0) << "\nStarting " << simulation_name << " simulation" << std::endl;
    t_tot.tic();
    superblock->H->update_timestep(settings::itebd::delta_t0, 1);
    while(delta_t > delta_tmin  and iteration < max_steps) {
        reduce_timestep();
        store_table_entry();
        single_TEBD_step(chi_max);
        swap();
        iteration++;
        phys_time += superblock->H->timestep;
        print_status_update();
    }
    t_tot.toc();
    print_status_full();
    print_profiling();
}
void class_imaginary_TEBD::reduce_timestep(){
    t_udt.tic();
    old_entropy = new_entropy;
    new_entropy = measurement->get_entropy();
    double rel_diff = std::fabs((new_entropy-old_entropy)/new_entropy);
    if (rel_diff < 1e-6 ){
        delta_t *=0.5;
        superblock->H->update_timestep(delta_t, 1);
    }
    t_udt.toc();
}


void class_imaginary_TEBD::store_table_entry(){
    t_sto.tic();
    table_buffer->emplace_back(superblock->chi,
                               chi_max,
                               measurement->get_energy(),
                               measurement->get_entropy(),
                               measurement->get_variance1(),
                               measurement->get_variance2(),
                               measurement->get_truncation_error(),
                               iteration,
                               superblock->chain_length,
                               superblock->H->timestep,
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
    ccout(1) << left  << simulation_name << " ";
    ccout(1) << left  << "Step: "                   << setw(10) << iteration;
    ccout(1) << left  << "E: "                      << setw(21) << setprecision(16)    << fixed   << measurement->get_energy();
    ccout(1) << left  << "S: "                      << setw(21) << setprecision(16)    << fixed   << measurement->get_entropy();
    ccout(1) << left  << "\u03C7_max: "             << setw(4)  << setprecision(3)     << fixed   << chi_max;
    ccout(1) << left  << "\u03C7: "                 << setw(4)  << setprecision(3)     << fixed   << measurement->get_chi();
    ccout(1) << left  << "Var(E): "                 << setw(21) << setprecision(16)    << fixed   << measurement->get_variance();
    ccout(1) << left  << "SVD truncation: "         << setw(21) << setprecision(16)    << fixed   << measurement->get_truncation_error();
    ccout(1) << left  << "Timestep: "               << setw(25) << setprecision(12)    << fixed   << superblock->H->timestep;
    ccout(1) << std::endl;
    t_obs.toc();
}

void class_imaginary_TEBD::print_status_full(){
    std::cout << std::endl;
    std::cout << " -- Final results -- " << simulation_name << std::endl;
    ccout(0)  << setw(20) << "Energy               = " << setprecision(16) << fixed      << measurement->get_energy()        << std::endl;
    ccout(0)  << setw(20) << "Entanglement Entropy = " << setprecision(16) << fixed      << measurement->get_entropy()       << std::endl;
    ccout(0)  << setw(20) << "chi_max              = " << setprecision(4)  << fixed      << chi_max                          << std::endl;
    ccout(0)  << setw(20) << "chi                  = " << setprecision(4)  << fixed      << measurement->get_chi()           << std::endl;
    ccout(0)  << setw(20) << "Variance             = " << setprecision(16) << scientific << measurement->get_variance()      << std::endl;
    ccout(0)  << setw(20) << "SVD truncation       = " << setprecision(4)  << scientific << measurement->get_truncation_error() << std::endl;
    ccout(0)  << setw(20) << "Timestep             = " << setprecision(12) << fixed      << superblock->H->timestep << std::endl;
    std::cout << std::endl;
}
