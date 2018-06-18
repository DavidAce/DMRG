//
// Created by david on 2018-01-31.
//


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
#include <algorithms/class_iTEBD.h>
#include <algorithms/class_FES_iTEBD.h>

using namespace std;
using namespace Textra;

class_FES_iTEBD::class_FES_iTEBD(std::shared_ptr<class_hdf5_file> hdf5_) {
    hdf5 = std::move(hdf5_);
    sim_type    = SimulationType ::FES_iTEBD;
    sim_name    = "FES_iTEBD";
    table_name  = "FES_iTEBD";
    initialize_constants();
}

void class_FES_iTEBD::run() {
    if (!settings::fes_itebd::on) { return; }
    t_tot.tic();
    auto chi_max_list = Math::LinSpaced(chi_num, chi_min, chi_max);
    for(auto &chi_max : chi_max_list ) {
        std::string table_name_chi =  table_name +"_" + to_string(chi_max);
        class_iTEBD iTEBD(hdf5, sim_name, table_name_chi, SimulationType::FES_iTEBD);
        iTEBD.chi_max           = chi_max;
        iTEBD.run();
    }
}



void class_FES_iTEBD::print_profiling(){
    if (settings::profiling::on) {
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_env.print_time_w_percent(t_tot);
        t_prt.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        print_profiling_sim(t_sim);
        measurement->print_profiling(t_obs);
    }
}

void class_FES_iTEBD::print_profiling_sim(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "\n Simulation breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_opt.print_time_w_percent(t_parent);
        t_svd.print_time_w_percent(t_parent);
        t_mps.print_time_w_percent(t_parent);
        t_chi.print_time_w_percent(t_parent);
    }
}