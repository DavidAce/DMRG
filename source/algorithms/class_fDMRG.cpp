//
// Created by david on 2018-01-31.
//


#include <iomanip>
#include <IO/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_finite_chain_sweeper.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include "class_fDMRG.h"

using namespace std;
using namespace Textra;

class_fDMRG::class_fDMRG(std::shared_ptr<class_hdf5_file> hdf5_)
        : class_base_algorithm(std::move(hdf5_),"fDMRG", SimulationType::fDMRG) {
    initialize_constants();
    table_fdmrg = std::make_unique<class_hdf5_table<class_table_dmrg>>(hdf5, sim_name,sim_name);
    env_storage = std::make_shared<class_finite_chain_sweeper>(max_length, superblock, hdf5 );
    measurement  = std::make_shared<class_measurement>(superblock, env_storage, sim_type);
    initialize_state(settings::model::initial_state);
}



void class_fDMRG::run() {
    if (!settings::fdmrg::on) { return; }
    ccout(0) << "\nStarting " << sim_name << " simulation" << std::endl;
    t_tot.tic();
    initialize_chain();
    while(true) {
        single_DMRG_step(chi_temp);
        env_storage_overwrite_MPS();         //Needs to occurr after update_MPS...
        store_table_entry_to_file();
        store_profiling_to_file();
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if(iteration >= max_sweeps) {break;}
        enlarge_environment(env_storage->get_direction());
        env_storage_move();
        update_chi();
        iteration = env_storage->get_sweeps();
    }
    t_tot.toc();
    print_status_full();
    print_profiling();
    measurement->compute_finite_norm();
    measurement->compute_finite_energy();
    measurement->compute_energy_variance_finite_chain();

}

void class_fDMRG::initialize_chain() {
    while(true){
        single_DMRG_step(chi_max);
        print_status_update();
        env_storage_insert();
        if (superblock->environment_size + 2ul < (unsigned long) max_length) {
            enlarge_environment();
            swap();
        } else {
            break;
        }
    }
}


void class_fDMRG::initialize_constants(){
    using namespace settings;
    max_length   = fdmrg::max_length;
    max_sweeps   = fdmrg::max_sweeps;
    chi_max      = fdmrg::chi_max;
    chi_grow     = fdmrg::chi_grow;
    print_freq   = fdmrg::print_freq;
    store_freq   = fdmrg::store_freq;
}


void class_fDMRG::update_chi(){
    t_chi.tic();
    if(entropy_has_converged()) {
        simulation_has_converged = chi_temp == chi_max;
        chi_temp = chi_grow ? std::min(chi_max, chi_temp + 4) : chi_max;
    }
    if(not chi_grow){
        chi_temp = chi_max;
    }
    t_chi.toc();
}


void class_fDMRG::store_table_entry_to_file(){
    if (Math::mod(iteration, store_freq) != 0) {return;}
    if (not env_storage->position_is_the_middle()) {return;}
    if (store_freq == 0){return;}
    compute_observables();
    t_sto.tic();
    table_fdmrg->append_record(
            iteration,
            measurement->get_chain_length(),
            env_storage->get_position(),
            measurement->get_chi(),
            chi_max,
            measurement->get_energy1(),
            measurement->get_energy2(),
            measurement->get_energy3(),
            measurement->get_energy4(),
            measurement->get_energy5(),
            measurement->get_energy6(),
            measurement->get_variance1(),
            measurement->get_variance2(),
            measurement->get_variance3(),
            measurement->get_variance4(),
            measurement->get_variance5(),
            measurement->get_variance6(),
            measurement->get_entanglement_entropy(),
            measurement->get_truncation_error(),
            t_tot.get_age());
    t_sto.toc();
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