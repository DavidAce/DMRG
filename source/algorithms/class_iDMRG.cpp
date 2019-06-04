//
// Created by david on 2018-01-18.
//

#include <iomanip>
#include <h5pp/h5pp.h>
#include <io/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/nmspc_mps_tools.h>
#include <general/nmspc_math.h>
#include <spdlog/spdlog.h>
#include "class_iDMRG.h"
using namespace std;
using namespace Textra;

class_iDMRG::class_iDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
    : class_algorithm_base(std::move(h5ppFile_),"iDMRG", SimulationType::iDMRG) {
    table_idmrg       = std::make_unique<class_hdf5_table<class_table_dmrg>>(h5ppFile, sim_name + "/measurements", "simulation_progress",sim_name);
//    initialize_superblock(settings::model::initial_state);

}



void class_iDMRG::run() {
    if (!settings::idmrg::on) { return; }
    log->info("Starting {} simulation", sim_name);
    t_tot.tic();
    while(true){
        single_DMRG_step();
        print_status_update();
        store_table_entry_progress();
        store_profiling_to_file_delta();
        check_convergence();

        // It's important not to perform the last swap.
        // That last state would not get optimized

        if (sim_state.iteration >= settings::idmrg::max_steps)  {stop_reason = StopReason::MAX_STEPS; break;}
        if (sim_state.simulation_has_converged)                 {stop_reason = StopReason::CONVERGED; break;}
        if (sim_state.simulation_has_to_stop)                   {stop_reason = StopReason::SATURATED; break;}


        enlarge_environment();
        swap();
        sim_state.iteration++;
    }
    t_tot.toc();
    switch(stop_reason){
        case StopReason::MAX_STEPS : log->info("Finished {} simulation -- reason: MAX_STEPS",sim_name) ;break;
        case StopReason::CONVERGED : log->info("Finished {} simulation -- reason: CONVERGED",sim_name) ;break;
        case StopReason::SATURATED : log->info("Finished {} simulation -- reason: SATURATED",sim_name) ;break;
        default: log->info("Finished {} simulation -- reason: NONE GIVEN",sim_name);
    }


    print_status_full();
    print_profiling();
    superblock->t_eig.print_time();
    h5ppFile->writeDataset(true, sim_name + "/simOK");
}

void class_iDMRG::run_simulation()    {}
void class_iDMRG::run_preprocessing() {}
void class_iDMRG::run_postprocessing(){}

void class_iDMRG::store_state_and_measurements_to_file(bool force){
    if(not force){
        if (Math::mod(sim_state.iteration, settings::idmrg::store_freq) != 0) {return;}
        if (settings::fdmrg::store_freq == 0){return;}
    }
    log->trace("Storing storing mps to file");
    t_sto.tic();
    MPS_Tools::Infinite::H5pp::write_all_superblock(*superblock, *h5ppFile, sim_name);
    t_sto.toc();
}

void class_iDMRG::store_table_entry_progress(bool force){
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