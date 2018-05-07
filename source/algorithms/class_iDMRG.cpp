//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <sim_parameters/nmspc_sim_settings.h>
#include <IO/class_hdf5_table_buffer.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <general/nmspc_math.h>
#include "class_iDMRG.h"
using namespace std;
using namespace Textra;

class_iDMRG::class_iDMRG(std::shared_ptr<class_hdf5_file> hdf5_)
    : class_base_algorithm(std::move(hdf5_),"iDMRG","iDMRG", SimulationType::iDMRG) {
}



void class_iDMRG::run() {
    if (!settings::idmrg::on) { return; }
    ccout(0) << "\nStarting " << table_name << " simulation" << std::endl;
        t_tot.tic();
        while(iteration < max_steps){
            single_DMRG_step(chi_temp);
            print_status_update();
            store_table_entry();
            position = enlarge_environment();
            iteration++;
            update_chi();
            swap();
        }
        t_tot.toc();
        print_status_full();
        print_profiling();
        superblock->t_eig.print_time();

}

void class_iDMRG::print_profiling(){
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

void class_iDMRG::print_profiling_sim(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "\n Simulation breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_opt.print_time_w_percent(t_parent);
        t_svd.print_time_w_percent(t_parent);
        t_env.print_time_w_percent(t_parent);
        t_mps.print_time_w_percent(t_parent);
        t_chi.print_time_w_percent(t_parent);
    }
}