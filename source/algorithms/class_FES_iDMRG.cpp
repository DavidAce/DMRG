
//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <sim_parameters/nmspc_sim_settings.h>
#include <IO/class_hdf5_table_buffer.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <general/nmspc_math.h>
#include "class_FES_iDMRG.h"
#include "class_iDMRG.h"
#include <thread>
#include <chrono>
using namespace std;
using namespace Textra;

class_FES_iDMRG::class_FES_iDMRG(std::shared_ptr<class_hdf5_file> hdf5_) {
    hdf5        = std::move(hdf5_);
    sim_type    = SimulationType ::FES_iDMRG;
    sim_name    = "FES_iDMRG";
    table_name  = "FES_iDMRG";
    initialize_constants();
}


void class_FES_iDMRG::run() {
    if (!settings::fes_idmrg::on) { return; }
    t_tot.tic();
    auto chi_max_list = Math::LinSpaced(chi_num, chi_min, chi_max);
    for(auto &chi_max : chi_max_list ) {
        std::string table_name_chi =  table_name +"_" + to_string(chi_max);
        class_iDMRG iDMRG(hdf5, sim_name, table_name_chi, SimulationType::FES_iDMRG);
        iDMRG.chi_max           = chi_max;
        iDMRG.run();
    }
}



void class_FES_iDMRG::print_profiling(){
    if (settings::profiling::on) {
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_env.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        t_eig.print_time_w_percent(t_sim);
        t_svd.print_time_w_percent(t_sim);
        t_mps.print_time_w_percent(t_sim);
    }
}
