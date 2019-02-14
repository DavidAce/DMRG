//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <IO/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/nmspc_mps_tools.h>
#include <general/nmspc_math.h>
#include <spdlog/spdlog.h>
#include "class_iDMRG.h"
using namespace std;
using namespace Textra;

class_iDMRG::class_iDMRG(std::shared_ptr<class_hdf5_file> hdf5_)
    : class_algorithm_base(std::move(hdf5_),"iDMRG", SimulationType::iDMRG) {
    table_idmrg = std::make_unique<class_hdf5_table<class_table_dmrg>>(hdf5, sim_name,sim_name);
    initialize_state(settings::model::initial_state);

}



void class_iDMRG::run() {
    if (!settings::idmrg::on) { return; }
    ccout(0) << "\nStarting " << sim_name << " simulation" << std::endl;
    t_tot.tic();
    while(sim_state.iteration < settings::idmrg::max_steps){// and not simulation_has_converged){
        single_DMRG_step();
        print_status_update();
        store_progress_to_file();
        store_profiling_to_file_delta();
        check_convergence();
        enlarge_environment();
        swap();
        sim_state.iteration++;
    }
    t_tot.toc();
    print_status_full();
    print_profiling();
    superblock->t_eig.print_time();
}

void class_iDMRG::run_simulation()    {}
void class_iDMRG::run_preprocessing() {}
void class_iDMRG::run_postprocessing(){}

void class_iDMRG::store_state_to_file(bool force){
    if(not force){
        if (Math::mod(sim_state.iteration, settings::idmrg::store_freq) != 0) {return;}
        if (settings::fdmrg::store_freq == 0){return;}
    }
    spdlog::trace("Storing storing mps to file");
    t_sto.tic();
    MPS_Tools::Infinite::Hdf5::write_superblock_state(*superblock,*hdf5,sim_name);
    t_sto.toc();
}

void class_iDMRG::store_progress_to_file(bool force){
    if (not force){
        if (Math::mod(sim_state.iteration, settings::idmrg::store_freq) != 0) {return;}
    }
    compute_observables();
    using namespace MPS_Tools::Common::Measure;
    t_sto.tic();
    table_idmrg->append_record(
            sim_state.iteration,
            superblock->measurements.length,
            sim_state.iteration,
            superblock->measurements.bond_dimension,
            settings::idmrg::chi_max,
            superblock->measurements.energy_per_site_mpo,
            superblock->measurements.energy_per_site_ham,
            superblock->measurements.energy_per_site_mom,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            superblock->measurements.energy_variance_per_site_mpo ,
            superblock->measurements.energy_variance_per_site_ham,
            superblock->measurements.energy_variance_per_site_mom,
            superblock->measurements.current_entanglement_entropy,
            superblock->measurements.truncation_error,
            t_tot.get_age());


    t_sto.toc();
}


long   class_iDMRG::chi_max()   {return settings::idmrg::chi_max;}
int    class_iDMRG::num_sites() {return 2;}
int    class_iDMRG::store_freq(){return settings::idmrg::store_freq;}
int    class_iDMRG::print_freq(){return settings::idmrg::print_freq;}
bool   class_iDMRG::chi_grow()  {return settings::idmrg::chi_grow;}



void class_iDMRG::print_profiling(){
    if (settings::profiling::on) {
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_prt.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        print_profiling_sim(t_sim);
        superblock->print_profiling(t_obs);
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