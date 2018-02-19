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
#include "class_iTEBD.h"

using namespace std;
using namespace Textra;

class_iTEBD::class_iTEBD(std::shared_ptr<class_hdf5_file> hdf5_)
        : class_algorithm_base(std::move(hdf5_),"iTEBD","iTEBD", SimulationType::iTEBD) {
    delta_t      = delta_t0;
}


void class_iTEBD::run() {
    if (!settings::itebd::on) { return; }
    ccout(0) << "\nStarting " << table_name << " simulation" << std::endl;
    t_tot.tic();
    delta_t = delta_t0;
    superblock->H->update_timestep(delta_t0, 1);
    while(iteration < max_steps and not simulation_has_converged ) {
        single_TEBD_step(chi_temp);
        print_status_update();
        store_table_entry();

        iteration++;
        position ++;
        superblock->chain_length  += 2;
        phys_time += superblock->H->timestep;
        update_chi();
        swap();
    }
    t_tot.toc();
    print_status_full();
    print_profiling();
}



void class_iTEBD::print_profiling(){
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
