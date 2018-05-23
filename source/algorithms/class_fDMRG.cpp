//
// Created by david on 2018-01-31.
//


#include <iomanip>
#include <sim_parameters/nmspc_sim_settings.h>
#include <IO/class_hdf5_table_buffer.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_finite_chain_storage.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include "class_fDMRG.h"

using namespace std;
using namespace Textra;

class_fDMRG::class_fDMRG(std::shared_ptr<class_hdf5_file> hdf5_)
        : class_base_algorithm(std::move(hdf5_),"fDMRG", SimulationType::fDMRG) {
}



void class_fDMRG::run() {
    if (!settings::fdmrg::on) { return; }
    ccout(0) << "\nStarting " << sim_name << " simulation" << std::endl;
    t_tot.tic();
//    chi_temp = chi_max;
    initialize_chain();
    while(true) {
        single_DMRG_step(chi_temp);
        env_storage_overwrite_MPS();         //Needs to occurr after update_MPS...
        store_table_entry();
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if(sweeps >= max_sweeps) {break;}
        position = enlarge_environment(direction);
        position = env_storage_move();
        update_chi();
    }
    t_tot.toc();
    print_status_full();
    print_profiling();
    measurement->compute_finite_norm();
    measurement->compute_finite_energy();
    measurement->compute_finite_variance();

}

int class_fDMRG::initialize_chain() {
    while(true){
        single_DMRG_step(chi_max);
        print_status_update();
        position = env_storage_insert();
        if (superblock->environment_size + 2ul < (unsigned long) max_length) {
            position = enlarge_environment();
            swap();
        } else {
            break;
        }
    }
    return position;
}



void class_fDMRG::print_profiling(){
    if (settings::profiling::on) {
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_ste.print_time_w_percent(t_tot);
        t_prt.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        std::cout << std::endl;
        print_profiling_sim(t_sim);
        measurement->print_profiling(t_obs);
        std::cout << std::endl;
    }
}

void class_fDMRG::print_profiling_sim(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "Simulation breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_opt.print_time_w_percent(t_parent);
        t_svd.print_time_w_percent(t_parent);
        t_env.print_time_w_percent(t_parent);
        t_mps.print_time_w_percent(t_parent);
        t_chi.print_time_w_percent(t_parent);
        std::cout << std::endl;
    }
}