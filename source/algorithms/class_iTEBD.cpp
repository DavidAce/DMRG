//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <IO/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mpo.h>
#include <general/nmspc_math.h>
#include "class_iTEBD.h"
using namespace std;
using namespace Textra;
using namespace std::complex_literals;

class_iTEBD::class_iTEBD(std::shared_ptr<class_hdf5_file> hdf5_)
        : class_algorithm_base(std::move(hdf5_),"iTEBD", SimulationType::iTEBD) {
    initialize_constants();
    table_itebd = std::make_unique<class_hdf5_table<class_table_tebd>>(hdf5, sim_name,sim_name);
    delta_t      = delta_t0;
    measurement  = std::make_shared<class_measurement>(superblock, sim_type);
    initialize_state(settings::model::initial_state);

}


void class_iTEBD::run() {
    if (!settings::itebd::on) { return; }
    ccout(0) << "\nStarting " << sim_name << " simulation" << std::endl;
    t_tot.tic();
    delta_t = delta_t0;
    superblock->H->update_evolution_step_size(-delta_t0, suzuki_order);
    while(iteration < max_steps and not simulation_has_converged) {
        single_TEBD_step(chi_max_temp);
        phys_time += superblock->H->step_size;
        store_table_entry_to_file();
        store_profiling_to_file();
        print_status_update();
        check_convergence_overall();
        iteration++;
    }
    t_tot.toc();
    print_status_full();
    print_profiling();
}

void class_iTEBD::initialize_constants(){
    using namespace settings;
    max_steps    = itebd::max_steps   ;
    delta_t0     = itebd::delta_t0    ;
    delta_tmin   = itebd::delta_tmin  ;
    suzuki_order = itebd::suzuki_order;
    chi_max      = itebd::chi_max     ;
    chi_grow     = itebd::chi_grow    ;
    print_freq   = itebd::print_freq  ;
    store_freq   = itebd::store_freq  ;
}



void class_iTEBD::check_convergence_time_step(){
    if(delta_t <= delta_tmin){
        time_step_has_converged = true;
    }else if (bond_dimension_has_converged and entanglement_has_converged) {
        delta_t = std::max(delta_tmin, delta_t * 0.5);
        superblock->H->update_evolution_step_size(-delta_t, suzuki_order);
        clear_convergence_checks();
    }
}

void class_iTEBD::check_convergence_overall(){
    t_con.tic();
    check_convergence_entanglement();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    check_convergence_bond_dimension();
    check_convergence_time_step();
    if(entanglement_has_converged and
       variance_ham_has_converged and
       variance_mom_has_converged and
       bond_dimension_has_converged and
       time_step_has_converged)
    {
        simulation_has_converged = true;
    }
    t_con.toc();
}


void class_iTEBD::store_table_entry_to_file(){
    if (Math::mod(iteration, store_freq) != 0) {return;}
    compute_observables();
    t_sto.tic();
    table_itebd->append_record(
            iteration,
            measurement->get_chi(),
            chi_max,
            delta_t,
            std::numeric_limits<double>::quiet_NaN(),
            measurement->get_energy_ham(),
            measurement->get_energy_mom(),
            std::numeric_limits<double>::quiet_NaN(),
            measurement->get_variance_ham(),
            measurement->get_variance_mom(),
            measurement->get_entanglement_entropy(),
            measurement->get_truncation_error(),
            phys_time,
            t_tot.get_age());

    t_sto.toc();
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
        t_con.print_time_w_percent(t_parent);
        t_udt.print_time_w_percent(t_parent);
    }
}