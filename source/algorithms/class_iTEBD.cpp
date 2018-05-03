//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <complex>
#include <sim_parameters/nmspc_sim_settings.h>
#include <IO/class_hdf5_table_buffer.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mpo.h>
#include <general/nmspc_math.h>
#include "class_iTEBD.h"
using namespace std;
using namespace Textra;
using namespace std::complex_literals;

class_iTEBD::class_iTEBD(std::shared_ptr<class_hdf5_file> hdf5_)
        : class_base_algorithm(std::move(hdf5_),"iTEBD","iTEBD", SimulationType::iTEBD) {
    delta_t      = delta_t0;
}


void class_iTEBD::run() {
    if (!settings::itebd::on) { return; }
    ccout(0) << "\nStarting " << table_name << " simulation" << std::endl;
    t_tot.tic();
    delta_t = delta_t0;
    superblock->H->update_evolution_step_size(-delta_t0, suzuki_order);
    while(iteration < max_steps ) {
        single_TEBD_step(chi_temp);
        iteration++;
        phys_time += superblock->H->step_size;
        position = enlarge_environment();
        store_table_entry();
        print_status_update();
        update_chi();
    }
    t_tot.toc();
    print_status_full();
    print_profiling();
}



void class_iTEBD::print_profiling(){
    if (settings::profiling::on) {
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_prt.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        print_profiling_sim(t_sim);
        measurement->print_profiling(t_obs);
    }
}

void class_iTEBD::print_profiling_sim(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "\n Simulation breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_evo.print_time_w_percent(t_parent);
        t_svd.print_time_w_percent(t_parent);
        t_env.print_time_w_percent(t_parent);
        t_chi.print_time_w_percent(t_parent);
        t_udt.print_time_w_percent(t_parent);
    }
}