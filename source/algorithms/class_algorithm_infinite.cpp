//
// Created by david on 2019-06-24.
//

#include "class_algorithm_infinite.h"
#include <mps_state/class_superblock.h>
#include <mps_tools/nmspc_mps_tools.h>
#include <io/class_hdf5_table_buffer2.h>
#include <general/nmspc_math.h>
#include <h5pp/h5pp.h>


class_algorithm_infinite::class_algorithm_infinite(
        std::shared_ptr<h5pp::File> h5ppFile_,
        std::string sim_name,
        SimulationType sim_type
        )
    : class_algorithm_base(std::move(h5ppFile_),sim_name, sim_type)
{
    table_dmrg     = std::make_unique<class_hdf5_table<class_table_dmrg>>(h5pp_file, sim_name + "/measurements", "simulation_progress", sim_name);
    superblock     = std::make_shared<class_superblock>(sim_type,sim_name);
}




void class_algorithm_infinite::run() {
    if (not sim_on()) { return; }
    run_preprocessing();
    run_simulation();
    run_postprocessing();
}

void class_algorithm_infinite::run_preprocessing() {

}

void class_algorithm_infinite::run_postprocessing(){
    print_status_full();
    print_profiling();
    h5pp_file->writeDataset(true, sim_name + "/simOK");
}

void class_algorithm_infinite::compute_observables(){
    log->trace("Starting all measurements on current superblock");
    t_sim.tic();
    t_obs.tic();
    superblock->do_all_measurements();
    t_obs.toc();
    t_sim.toc();
}


void class_algorithm_infinite::reset_to_random_state(const std::string parity) {
    log->trace("Resetting MPS to random product state");
    sim_state.iteration = 0;

    // Randomize state
    *superblock = mpstools::infinite::mps::set_random_state(*superblock,parity);
    clear_saturation_status();
}


void class_algorithm_infinite::enlarge_environment(){
    log->trace("Enlarging environment" );
    t_sim.tic();
    t_env.tic();
    superblock->enlarge_environment(0);
    t_env.toc();
    t_sim.toc();
}

void class_algorithm_infinite::swap(){
    log->trace("Swap AB sites on superblock");
    t_sim.tic();
    superblock->swap_AB();
    t_sim.toc();
}



void class_algorithm_infinite::check_convergence_variance_ham(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    log->trace("Checking convergence of variance ham");

    threshold       = std::isnan(threshold)       ? settings::precision::VarConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::VarSaturationThreshold  : slope_threshold;
    check_saturation_using_slope(B_ham_vec,
                                 V_ham_vec,
                                 X_ham_vec,
                                 mpstools::common::measure::energy_variance_per_site_ham(*superblock),
                                 sim_state.iteration,
                                 1,
                                 slope_threshold,
                                 V_ham_slope,
                                 sim_state.variance_ham_has_saturated);
    sim_state.variance_ham_has_converged = mpstools::common::measure::energy_variance_per_site_ham(*superblock) < threshold;
}

void class_algorithm_infinite::check_convergence_variance_mom(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    log->trace("Checking convergence of variance mom");

    threshold       = std::isnan(threshold)       ? settings::precision::VarConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::VarSaturationThreshold  : slope_threshold;
    check_saturation_using_slope(B_mom_vec,
                                 V_mom_vec,
                                 X_mom_vec,
                                 mpstools::common::measure::energy_variance_per_site_mom(*superblock),
                                 sim_state.iteration,
                                 1,
                                 slope_threshold,
                                 V_mom_slope,
                                 sim_state.variance_mom_has_saturated);
    sim_state.variance_mom_has_converged = mpstools::common::measure::energy_variance_per_site_mom(*superblock) < threshold;
}

void class_algorithm_infinite::check_convergence_entg_entropy(double slope_threshold) {
    //Based on the the slope of entanglement middle_entanglement_entropy
    // This one is cheap to compute.
    log->debug("Checking convergence of entanglement");

    slope_threshold = std::isnan(slope_threshold) ? settings::precision::EntEntrSaturationThreshold  : slope_threshold;
    check_saturation_using_slope(BS_vec,
                                 S_vec,
                                 XS_vec,
                                 mpstools::common::measure::current_entanglement_entropy(*superblock),
                                 sim_state.iteration,
                                 1,
                                 slope_threshold,
                                 S_slope,
                                 sim_state.entanglement_has_saturated);
    sim_state.entanglement_has_converged = sim_state.entanglement_has_saturated;
}



void class_algorithm_infinite::store_state_and_measurements_to_file(bool force){
    if(not force){
        if (Math::mod(sim_state.iteration, store_freq()) != 0) {return;}
        if (store_freq() == 0){return;}
    }
    log->trace("Storing storing mps to file");
    t_sto.tic();
    mpstools::infinite::io::write_all_superblock(*superblock, *h5pp_file, sim_name);
    t_sto.toc();
}




//void class_algorithm_infinite::store_table_entry_progress(bool force){
//    if (not force){
//        if (Math::mod(sim_state.iteration, settings::idmrg::store_freq) != 0) {return;}
//    }
//    compute_observables();
//    using namespace mpstools::common::measure;
//    t_sto.tic();
//    table_dmrg->append_record(
//            sim_state.iteration,
//            superblock->measurements.length.value(),
//            sim_state.iteration,
//            superblock->measurements.bond_dimension.value(),
//            settings::idmrg::chi_max,
//            superblock->measurements.energy_per_site.value(),
//            superblock->measurements.energy_per_site_ham.value(),
//            superblock->measurements.energy_per_site_mom.value(),
//            std::numeric_limits<double>::quiet_NaN(),
//            std::numeric_limits<double>::quiet_NaN(),
//            std::numeric_limits<double>::quiet_NaN(),
//            superblock->measurements.energy_variance_per_site.value(),
//            superblock->measurements.energy_variance_per_site_ham.value(),
//            superblock->measurements.energy_variance_per_site_mom.value(),
//            superblock->measurements.current_entanglement_entropy.value(),
//            superblock->measurements.truncation_error.value(),
//            t_tot.get_age());
//
//
//    t_sto.toc();
//}


void class_algorithm_infinite::clear_saturation_status(){
    log->trace("Clearing saturation status");

    BS_vec.clear();
    S_vec.clear();
    XS_vec.clear();

    B_mpo_vec.clear();
    V_mpo_vec.clear();
    X_mpo_vec.clear();
    B_ham_vec.clear();
    V_ham_vec.clear();
    X_ham_vec.clear();
    B_mom_vec.clear();
    V_mom_vec.clear();
    X_mom_vec.clear();

    sim_state.entanglement_has_saturated      = false;
    sim_state.variance_mpo_has_saturated      = false;
    sim_state.variance_ham_has_saturated      = false;
    sim_state.variance_mom_has_saturated      = false;

    sim_state.variance_mpo_saturated_for = 0;
    sim_state.variance_ham_saturated_for = 0;
    sim_state.variance_mom_saturated_for = 0;



    sim_state.entanglement_has_converged = false;
    sim_state.variance_mpo_has_converged = false;
    sim_state.variance_ham_has_converged = false;
    sim_state.variance_mom_has_converged = false;

    sim_state.bond_dimension_has_reached_max = false;
    sim_state.simulation_has_to_stop         = false;
}



void class_algorithm_infinite::print_profiling(){
    if (settings::profiling::on) {
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_prt.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        print_profiling_sim(t_sim);
   }
}

void class_algorithm_infinite::print_profiling_sim(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "\n Simulation breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_opt.print_time_w_percent(t_parent);
        t_evo.print_time_w_percent(t_parent);
        t_svd.print_time_w_percent(t_parent);
        t_env.print_time_w_percent(t_parent);
        t_mps.print_time_w_percent(t_parent);
        t_con.print_time_w_percent(t_parent);
        t_udt.print_time_w_percent(t_parent);
    }
}





void class_algorithm_infinite::print_status_update() {
    if (Math::mod(sim_state.iteration, print_freq()) != 0) {return;}
//    if (not state->position_is_the_middle()) {return;}
    if (print_freq() == 0) {return;}
    compute_observables();
    using namespace std;
    t_prt.tic();
    std::stringstream report;
    report << setprecision(16) << fixed << left;
    report << left  << sim_name << " ";
    report << left  << "Iter: "                       << setw(6) << sim_state.iteration;
    report << left  << "E: ";

    switch(sim_type) {
        case SimulationType::iDMRG:
            report << setw(21) << setprecision(16)    << fixed   << superblock->measurements.energy_per_site_mpo.value();
            report << setw(21) << setprecision(16)    << fixed   << superblock->measurements.energy_per_site_ham.value();
            report << setw(21) << setprecision(16)    << fixed   << superblock->measurements.energy_per_site_mom.value();
            break;
        case SimulationType::iTEBD:
            report << setw(21) << setprecision(16)    << fixed   << superblock->measurements.energy_per_site_ham.value();
            report << setw(21) << setprecision(16)    << fixed   << superblock->measurements.energy_per_site_mom.value();
            break;
        default: throw std::runtime_error("Wrong simulation type");

    }

    report << left  << "log₁₀ σ²(E): ";
    switch(sim_type) {
        case SimulationType::iDMRG:
            report << setw(12) << setprecision(4)    << fixed   << std::log10(superblock->measurements.energy_variance_per_site_mpo.value());
            report << setw(12) << setprecision(4)    << fixed   << std::log10(superblock->measurements.energy_variance_per_site_ham.value());
            report << setw(12) << setprecision(4)    << fixed   << std::log10(superblock->measurements.energy_variance_per_site_mom.value());
            break;
        case SimulationType::iTEBD:
            report << setw(12) << setprecision(4)    << fixed   << std::log10(superblock->measurements.energy_variance_per_site_ham.value());
            report << setw(12) << setprecision(4)    << fixed   << std::log10(superblock->measurements.energy_variance_per_site_mom.value());
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }


    report << left  << "S: "                          << setw(21) << setprecision(16)    << fixed   << superblock->measurements.current_entanglement_entropy.value();
    report << left  << "χmax: "                       << setw(4)  << setprecision(3)     << fixed   << chi_max();
    report << left  << "χ: "                          << setw(4)  << setprecision(3)     << fixed   << superblock->measurements.bond_dimension.value();
    report << left  << "log₁₀ trunc: "                << setw(10) << setprecision(4)     << fixed   << std::log10(superblock->measurements.truncation_error.value());
    report << left  << "Sites: "                      << setw(6)  << setprecision(1)     << fixed   << superblock->measurements.length.value();
    switch(sim_type){
        case SimulationType::iTEBD:
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    report << left  << " Convergence [";
    switch(sim_type){
        case SimulationType::iDMRG:
            report << left  << " S-"   << std::boolalpha << setw(6) << sim_state.entanglement_has_converged;
            report << left  << " σ²-"  << std::boolalpha << setw(6) << sim_state.variance_mpo_has_converged;
            break;
        case SimulationType::iTEBD:
            report << left  << " S-"  << std::boolalpha << setw(6) << sim_state.entanglement_has_converged;
            break;
        default: throw std::runtime_error("Wrong simulation type");

    }
    report << left  << "]";
    report << left  << " Saturation [";
    switch(sim_type){
        case SimulationType::iDMRG:
            report << left  << " σ²- " << setw(2) << sim_state.variance_mpo_saturated_for << " steps";
            report << left  << " S-"   << std::boolalpha << setw(6) << sim_state.entanglement_has_saturated;
            break;
        case SimulationType::iTEBD:
            report << left  << " S-"   << std::boolalpha << setw(6) << sim_state.entanglement_has_saturated;
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    report << left  << "]";
    report << left  << " Time: "                          << setw(10) << setprecision(2)    << fixed   << t_tot.get_age() ;
    report << left << " Memory [";
    report << left << "Rss: "     << process_memory_in_mb("VmRSS")<< " MB ";
    report << left << "RssPeak: "  << process_memory_in_mb("VmHWM")<< " MB ";
    report << left << "VmPeak: "  << process_memory_in_mb("VmPeak")<< " MB";
    report << left << "]";
    log->info(report.str());
    t_prt.toc();
}

void class_algorithm_infinite::print_status_full(){
    compute_observables();
    using namespace std;
    using namespace mpstools::common::measure;
    t_prt.tic();
    log->info("--- Final results  --- {} ---", sim_name);
    log->info("Iterations            = {:<16d}"    , sim_state.iteration);
    switch(sim_type){
        case SimulationType::iDMRG:
            log->info("Energy MPO            = {:<16.16f}" , superblock->measurements.energy_per_site_mpo.value());
            log->info("Energy HAM            = {:<16.16f}" , superblock->measurements.energy_per_site_ham.value());
            log->info("Energy MOM            = {:<16.16f}" , superblock->measurements.energy_per_site_mom.value());
            break;
        case SimulationType::iTEBD:
            log->info("Energy HAM            = {:<16.16f}" , superblock->measurements.energy_per_site_ham.value());
            log->info("Energy MOM            = {:<16.16f}" , superblock->measurements.energy_per_site_mom.value());
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    switch(sim_type){
        case SimulationType::iDMRG:
            log->info("log₁₀ σ²(E) MPO       = {:<16.16f}" , superblock->measurements.energy_per_site_mpo.value());
            log->info("log₁₀ σ²(E) HAM       = {:<16.16f}" , superblock->measurements.energy_per_site_ham.value());
            log->info("log₁₀ σ²(E) MOM       = {:<16.16f}" , superblock->measurements.energy_per_site_mom.value());
            break;
        case SimulationType::iTEBD:
            log->info("log₁₀ σ²(E) HAM       = {:<16.16f}" , superblock->measurements.energy_per_site_ham.value());
            log->info("log₁₀ σ²(E) MOM       = {:<16.16f}" , superblock->measurements.energy_per_site_mom.value());
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }

    log->info("Entanglement Entropy  = {:<16.16f}" , superblock->measurements.current_entanglement_entropy.value());
    log->info("χmax                  = {:<16d}"    , chi_max()                                            );
    log->info("χ                     = {:<16d}"    , superblock->measurements.bond_dimension.value()      );
    log->info("log₁₀ truncation:     = {:<16.16f}" , log10(superblock->measurements.truncation_error.value()));

    switch(sim_type){
        case SimulationType::iTEBD:
            log->info("δt                    = {:<16.16f}" , sim_state.delta_t);
            break;

        default: throw std::runtime_error("Wrong simulation type");
    }

    log->info("Simulation converged  = {:<}"    , sim_state.simulation_has_converged);

    switch(sim_type){
        case SimulationType::iDMRG:
            log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}" , S_slope,sim_state.entanglement_has_converged, sim_state.entanglement_has_saturated);
            log->info("σ² MPO slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_mpo_slope ,sim_state.variance_mpo_has_converged, sim_state.variance_mpo_has_saturated);
            log->info("σ² HAM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_ham_slope ,sim_state.variance_ham_has_converged, sim_state.variance_ham_has_saturated);
            log->info("σ² MOM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_mom_slope ,sim_state.variance_mom_has_converged, sim_state.variance_mom_has_saturated);
            break;
        case SimulationType::iTEBD:
            log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}" , S_slope,sim_state.entanglement_has_converged, sim_state.entanglement_has_saturated);
            log->info("σ² HAM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_ham_slope ,sim_state.variance_ham_has_converged, sim_state.variance_ham_has_saturated);
            log->info("σ² MOM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_mom_slope ,sim_state.variance_mom_has_converged, sim_state.variance_mom_has_saturated);
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}" , S_slope,sim_state.entanglement_has_converged, sim_state.entanglement_has_saturated);
    log->info("Time                  = {:<16.16f}" , t_tot.get_age());
    log->info("Peak memory           = {:<6.1f} MB" , process_memory_in_mb("VmPeak"));
    t_prt.toc();
}

