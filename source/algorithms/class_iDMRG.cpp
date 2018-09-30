//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <IO/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <general/nmspc_math.h>
#include "class_iDMRG.h"
using namespace std;
using namespace Textra;

class_iDMRG::class_iDMRG(std::shared_ptr<class_hdf5_file> hdf5_)
    : class_algorithm_base(std::move(hdf5_),"iDMRG", SimulationType::iDMRG) {
    initialize_constants();
    table_idmrg = std::make_unique<class_hdf5_table<class_table_dmrg>>(hdf5, sim_name,sim_name);
    measurement  = std::make_shared<class_measurement>(superblock, sim_type);
    initialize_state(settings::model::initial_state);

}



void class_iDMRG::run() {
    if (!settings::idmrg::on) { return; }
    ccout(0) << "\nStarting " << sim_name << " simulation" << std::endl;
    t_tot.tic();
    while(iteration < max_steps){// and not simulation_has_converged){
        single_DMRG_step(chi_temp);
        print_status_update();
        store_table_entry_to_file();
        store_profiling_to_file_delta();
        enlarge_environment();
        check_convergence_all();
        swap();
        iteration++;
    }
    t_tot.toc();
    print_status_full();
    print_profiling();
    superblock->t_eig.print_time();
}

void class_iDMRG::initialize_constants(){
    using namespace settings;
    max_steps    = idmrg::max_steps;
    chi_max      = idmrg::chi_max   ;
    chi_grow     = idmrg::chi_grow  ;
    print_freq   = idmrg::print_freq;
    store_freq   = idmrg::store_freq;

}

void class_iDMRG::store_table_entry_to_file(){
    if (Math::mod(iteration, store_freq) != 0) {return;}
    compute_observables();
    t_sto.tic();
    table_idmrg->append_record(
                             iteration,
                             measurement->get_chain_length(),
                             iteration,
                             measurement->get_chi(),
                             chi_max,
                             measurement->get_energy_mpo(),
                             measurement->get_energy_ham(),
                             measurement->get_energy_mom(),
                             std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN(),
                             measurement->get_variance_mpo(),
                             measurement->get_variance_ham(),
                             measurement->get_variance_mom(),
                             measurement->get_entanglement_entropy(),
                             measurement->get_truncation_error(),
                             t_tot.get_age());

    t_sto.toc();
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
        t_con.print_time_w_percent(t_parent);
    }
}