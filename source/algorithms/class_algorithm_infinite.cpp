//
// Created by david on 2019-06-24.
//

#include "class_algorithm_infinite.h"
#include <state/class_infinite_state.h>
#include <tools/nmspc_tools.h>
#include <io/class_hdf5_log_buffer.h>
#include <math/nmspc_math.h>
#include <h5pp/h5pp.h>


class_algorithm_infinite::class_algorithm_infinite(
        std::shared_ptr<h5pp::File> h5ppFile_,
        std::string sim_name,
        SimulationType sim_type
        )
    : class_algorithm_base(std::move(h5ppFile_),sim_name, sim_type)
{
    state      = std::make_unique<class_infinite_state>(sim_type,sim_name);
}




void class_algorithm_infinite::run() {
    if (not sim_on()) { return; }
    t_tot.tic();
    run_preprocessing();
    run_simulation();
    run_postprocessing();
    t_tot.toc();
}

void class_algorithm_infinite::run_preprocessing() {
    t_pre.tic();
    t_pre.toc();
}

void class_algorithm_infinite::run_postprocessing(){
    t_pos.tic();
    print_status_full();
    print_profiling();
    h5pp_file->writeDataset(true, sim_name + "/simOK");
    t_pos.toc();
}

void class_algorithm_infinite::compute_observables(){
    log->trace("Starting all measurements on current state");
    state->do_all_measurements();
}

void class_algorithm_infinite::update_bond_dimension_limit(std::optional<long> max_bond_dim){
    if(not max_bond_dim.has_value()) {
        log->debug("No max bond dim given, setting {}", chi_max());
        max_bond_dim = chi_max();
    }

    sim_status.chi_lim_has_reached_chi_max = state->get_chi_lim() == max_bond_dim;
    if(not sim_status.chi_lim_has_reached_chi_max){
        if(chi_grow()){
            // Here the settings specify to grow the bond dimension limit progressively during the simulation
            // Only do this if the simulation is stuck.
            if(sim_status.simulation_has_got_stuck){
                long chi_new_limit = std::min(max_bond_dim.value(), state->get_chi_lim() * 2);
                log->debug("Updating bond dimension limit {} -> {}", state->get_chi_lim(), chi_new_limit);
                state->set_chi_lim(chi_new_limit);
                clear_saturation_status();
            }else{
                log->debug("chi_grow is ON but sim is not stuck -> Kept current bond dimension limit {}", state->get_chi_lim());
            }
        }else{
            // Here the settings specify to just set the limit to maximum chi directly
            log->debug("Setting bond dimension limit to maximum = {}", chi_max());
            state->set_chi_lim(max_bond_dim.value());
        }
    }else{
        log->debug("Chi limit has reached max: {} -> Kept current bond dimension limit {}", chi_max(),state->get_chi_lim());
    }
    sim_status.chi_max = max_bond_dim.value();
    sim_status.chi_lim = state->get_chi_lim();
}



void class_algorithm_infinite::reset_to_random_state(const std::string parity, int seed_state) {
    log->trace("Resetting MPS to random product state");
    sim_status.iteration = 0;

    // Randomize state
    *state = tools::infinite::mps::set_random_state(*state,parity, seed_state);
    clear_saturation_status();
}


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

    sim_status.entanglement_has_saturated      = false;
    sim_status.variance_mpo_has_saturated      = false;
    sim_status.variance_ham_has_saturated      = false;
    sim_status.variance_mom_has_saturated      = false;

    sim_status.variance_mpo_saturated_for = 0;
    sim_status.variance_ham_saturated_for = 0;
    sim_status.variance_mom_saturated_for = 0;


    sim_status.entanglement_has_converged = false;
    sim_status.variance_mpo_has_converged = false;
    sim_status.variance_ham_has_converged = false;
    sim_status.variance_mom_has_converged = false;

    sim_status.chi_lim_has_reached_chi_max = false;
    sim_status.simulation_has_to_stop      = false;
}



void class_algorithm_infinite::enlarge_environment(){
    log->trace("Enlarging environment" );
    state->enlarge_environment(0);
}

void class_algorithm_infinite::swap(){
    log->trace("Swap AB sites on state");
    state->swap_AB();
}

void class_algorithm_infinite::check_convergence_variance_mpo(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    log->debug("Checking convergence of variance mpo");
    threshold       = std::isnan(threshold) ? settings::precision::varianceConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::varianceSlopeThreshold : slope_threshold;
    compute_observables();

    auto report = check_saturation_using_slope(
                    B_mpo_vec,
                    V_mpo_vec,
                    X_mpo_vec,
                    tools::infinite::measure::energy_variance_per_site_mpo(*state),
                    sim_status.iteration,
                    1,
                    slope_threshold);
    if(report.has_computed) V_mpo_slope  = report.slope;
    sim_status.variance_mpo_has_saturated = report.has_saturated;
    sim_status.variance_mpo_saturated_for = (int) count(B_mpo_vec.begin(), B_mpo_vec.end(), true);
    sim_status.variance_mpo_has_converged =  state->measurements.energy_variance_per_site_mpo.value() < threshold;

}

void class_algorithm_infinite::check_convergence_variance_ham(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    log->trace("Checking convergence of variance ham");

    threshold       = std::isnan(threshold) ? settings::precision::varianceConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::varianceSlopeThreshold : slope_threshold;
    auto report  = check_saturation_using_slope(
            B_ham_vec,
            V_ham_vec,
            X_ham_vec,
            tools::infinite::measure::energy_variance_per_site_ham(*state),
            sim_status.iteration,
            1,
            slope_threshold);
    if(report.has_computed) V_ham_slope  = report.slope;
    sim_status.variance_ham_has_saturated = report.has_saturated;
    sim_status.variance_ham_has_converged = tools::infinite::measure::energy_variance_per_site_ham(*state) < threshold;
}

void class_algorithm_infinite::check_convergence_variance_mom(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    log->trace("Checking convergence of variance mom");

    threshold       = std::isnan(threshold) ? settings::precision::varianceConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::varianceSlopeThreshold : slope_threshold;
    auto report = check_saturation_using_slope(B_mom_vec,
            V_mom_vec,
            X_mom_vec,
            tools::infinite::measure::energy_variance_per_site_mom(*state),
            sim_status.iteration,
            1,
            slope_threshold);
    if(report.has_computed) V_mom_slope  = report.slope;
    sim_status.variance_mom_has_saturated = report.has_saturated;
    sim_status.variance_mom_has_converged = tools::infinite::measure::energy_variance_per_site_mom(*state) < threshold;
}

void class_algorithm_infinite::check_convergence_entg_entropy(double slope_threshold) {
    //Based on the the slope of entanglement entanglement_entropy_midchain
    // This one is cheap to compute.
    log->debug("Checking convergence of entanglement");

    slope_threshold = std::isnan(slope_threshold) ? settings::precision::entropySlopeThreshold : slope_threshold;
    auto report = check_saturation_using_slope(
            BS_vec,
            S_vec,
            XS_vec,
            tools::infinite::measure::current_entanglement_entropy(*state),
            sim_status.iteration,
            1,
            slope_threshold);
    if(report.has_computed) S_slope       = report.slope;
    sim_status.entanglement_has_saturated = report.has_saturated;
    sim_status.entanglement_has_converged = sim_status.entanglement_has_saturated;
}

void class_algorithm_infinite::write_measurements(bool force){
    if(not force){
        if (math::mod(sim_status.iteration, write_freq()) != 0) {return;}
        if (write_freq() == 0){return;}
    }
    log->trace("Writing all measurements to file");
    state->unset_measurements();
    compute_observables();
    h5pp_file->writeDataset(false, sim_name + "/simOK");
    tools::infinite::io::write_all_measurements(*state, *h5pp_file, sim_name);
    h5pp_file->writeDataset(true, sim_name + "/simOK");
}

void class_algorithm_infinite::write_state(bool force){
    if(not force){
        if (math::mod(sim_status.iteration, write_freq()) != 0) {return;}
        if (write_freq() == 0){return;}
        if (settings::output::storage_level <= StorageLevel::NONE){return;}
    }
    log->trace("Writing state to file");
    h5pp_file->writeDataset(false, sim_name + "/simOK");
    tools::infinite::io::write_all_state(*state, *h5pp_file, sim_name);
    h5pp_file->writeDataset(true, sim_name + "/simOK");
}


void class_algorithm_infinite::write_status(bool force){
    if (not force){
        if (math::mod(sim_status.iteration, write_freq()) != 0) {return;}
        if (write_freq() == 0){return;}
        if (settings::output::storage_level <= StorageLevel::NONE){return;}
    }
    log->trace("Writing simulation status to file");
    h5pp_file->writeDataset(false, sim_name + "/simOK");
    tools::common::io::write_simulation_status(sim_status, *h5pp_file, sim_name);
    h5pp_file->writeDataset(true, sim_name + "/simOK");
}


//void class_algorithm_infinite::store_log_entry_progress(bool force){
//    if (not force){
//        if (math::mod(sim_status.iteration, settings::idmrg::write_freq) != 0) {return;}
//    }
//    compute_observables();
//    using namespace tools::infinite::measure;
//    t_sto.tic();
//    log_dmrg->append_record(
//            sim_status.iteration,
//            state->measurements.length.value(),
//            sim_status.iteration,
//            state->measurements.bond_dimension.value(),
//            settings::idmrg::chi_lim,
//            state->measurements.energy_per_site.value(),
//            state->measurements.energy_per_site_ham.value(),
//            state->measurements.energy_per_site_mom.value(),
//            std::numeric_limits<double>::quiet_NaN(),
//            std::numeric_limits<double>::quiet_NaN(),
//            std::numeric_limits<double>::quiet_NaN(),
//            state->measurements.energy_variance_per_site.value(),
//            state->measurements.energy_variance_per_site_ham.value(),
//            state->measurements.energy_variance_per_site_mom.value(),
//            state->measurements.current_entanglement_entropy.value(),
//            state->measurements.truncation_error.value(),
//            t_tot.get_age());
//
//
//    t_sto.toc();
//}



void class_algorithm_infinite::print_status_update() {
    if (math::mod(sim_status.iteration, print_freq()) != 0) {return;}
//    if (not state->position_is_the_middle()) {return;}
    if (print_freq() == 0) {return;}
    compute_observables();
    using namespace std;
    t_prt.tic();
    std::stringstream report;
    report << setprecision(16) << fixed << left;
    report << left  << sim_name << " ";
    report << left  << "Iter: "                       << setw(6) << sim_status.iteration;
    report << left  << "E: ";

    switch(sim_type) {
        case SimulationType::iDMRG:
            report << setw(21) << setprecision(16)    << fixed   << state->measurements.energy_per_site_mpo.value();
            report << setw(21) << setprecision(16)    << fixed   << state->measurements.energy_per_site_ham.value();
            report << setw(21) << setprecision(16)    << fixed   << state->measurements.energy_per_site_mom.value();
            break;
        case SimulationType::iTEBD:
            report << setw(21) << setprecision(16)    << fixed   << state->measurements.energy_per_site_ham.value();
            report << setw(21) << setprecision(16)    << fixed   << state->measurements.energy_per_site_mom.value();
            break;
        default: throw std::runtime_error("Wrong simulation type");

    }

    report << left  << "log₁₀ σ²(E): ";
    switch(sim_type) {
        case SimulationType::iDMRG:
            report << setw(12) << setprecision(4)    << fixed   << std::log10(state->measurements.energy_variance_per_site_mpo.value());
            report << setw(12) << setprecision(4)    << fixed   << std::log10(state->measurements.energy_variance_per_site_ham.value());
            report << setw(12) << setprecision(4)    << fixed   << std::log10(state->measurements.energy_variance_per_site_mom.value());
            break;
        case SimulationType::iTEBD:
            report << setw(12) << setprecision(4)    << fixed   << std::log10(state->measurements.energy_variance_per_site_ham.value());
            report << setw(12) << setprecision(4)    << fixed   << std::log10(state->measurements.energy_variance_per_site_mom.value());
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }


    report << left  << "S: "                          << setw(21) << setprecision(16)    << fixed   << state->measurements.current_entanglement_entropy.value();
    report << left  << "χmax: "                       << setw(4)  << setprecision(3)     << fixed   << chi_max();
    report << left  << "χ: "                          << setw(4)  << setprecision(3)     << fixed   << state->measurements.bond_dimension.value();
    report << left  << "log₁₀ trunc: "                << setw(10) << setprecision(4)     << fixed   << std::log10(state->measurements.truncation_error.value());
    report << left  << "Sites: "                      << setw(6)  << setprecision(1)     << fixed   << state->measurements.length.value();
    switch(sim_type){
        case SimulationType::iTEBD:
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    report << left  << " Convergence [";
    switch(sim_type){
        case SimulationType::iDMRG:
            report << left  << " S-"   << std::boolalpha << setw(6) << sim_status.entanglement_has_converged;
            report << left  << " σ²-"  << std::boolalpha << setw(6) << sim_status.variance_mpo_has_converged;
            break;
        case SimulationType::iTEBD:
            report << left  << " S-"  << std::boolalpha << setw(6) << sim_status.entanglement_has_converged;
            break;
        default: throw std::runtime_error("Wrong simulation type");

    }
    report << left  << "]";
    report << left  << " Saturation [";
    switch(sim_type){
        case SimulationType::iDMRG:
            report << left  << " σ²- " << setw(2) << sim_status.variance_mpo_saturated_for << " steps";
            report << left  << " S-"   << std::boolalpha << setw(6) << sim_status.entanglement_has_saturated;
            break;
        case SimulationType::iTEBD:
            report << left  << " S-"   << std::boolalpha << setw(6) << sim_status.entanglement_has_saturated;
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
    using namespace tools::infinite::measure;
    t_prt.tic();
    log->info("--- Final results  --- {} ---", sim_name);
    log->info("Iterations            = {:<16d}"    , sim_status.iteration);
    switch(sim_type){
        case SimulationType::iDMRG:
            log->info("Energy MPO            = {:<16.16f}" , state->measurements.energy_per_site_mpo.value());
            log->info("Energy HAM            = {:<16.16f}" , state->measurements.energy_per_site_ham.value());
            log->info("Energy MOM            = {:<16.16f}" , state->measurements.energy_per_site_mom.value());
            break;
        case SimulationType::iTEBD:
            log->info("Energy HAM            = {:<16.16f}" , state->measurements.energy_per_site_ham.value());
            log->info("Energy MOM            = {:<16.16f}" , state->measurements.energy_per_site_mom.value());
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    switch(sim_type){
        case SimulationType::iDMRG:
            log->info("log₁₀ σ²(E) MPO       = {:<16.16f}" , state->measurements.energy_per_site_mpo.value());
            log->info("log₁₀ σ²(E) HAM       = {:<16.16f}" , state->measurements.energy_per_site_ham.value());
            log->info("log₁₀ σ²(E) MOM       = {:<16.16f}" , state->measurements.energy_per_site_mom.value());
            break;
        case SimulationType::iTEBD:
            log->info("log₁₀ σ²(E) HAM       = {:<16.16f}" , state->measurements.energy_per_site_ham.value());
            log->info("log₁₀ σ²(E) MOM       = {:<16.16f}" , state->measurements.energy_per_site_mom.value());
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }

    log->info("Entanglement Entropy  = {:<16.16f}" , state->measurements.current_entanglement_entropy.value());
    log->info("χmax                  = {:<16d}"    , chi_max()                                            );
    log->info("χ                     = {:<16d}"    , state->measurements.bond_dimension.value()      );
    log->info("log₁₀ truncation:     = {:<16.16f}" , log10(state->measurements.truncation_error.value()));

    switch(sim_type){
        case SimulationType::iTEBD:
            log->info("δt                    = {:<16.16f}" , sim_status.delta_t);
            break;

        default: throw std::runtime_error("Wrong simulation type");
    }

    log->info("Simulation saturated  = {:<}"    , sim_status.simulation_has_saturated);
    log->info("Simulation converged  = {:<}"    , sim_status.simulation_has_converged);
    log->info("Simulation succeeded  = {:<}"    , sim_status.simulation_has_succeeded);
    log->info("Simulation got stuck  = {:<}"    , sim_status.simulation_has_got_stuck);
    switch(sim_type){
        case SimulationType::iDMRG:
            log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}" , S_slope,sim_status.entanglement_has_converged, sim_status.entanglement_has_saturated);
            log->info("σ² MPO slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_mpo_slope ,sim_status.variance_mpo_has_converged, sim_status.variance_mpo_has_saturated);
            log->info("σ² HAM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_ham_slope ,sim_status.variance_ham_has_converged, sim_status.variance_ham_has_saturated);
            log->info("σ² MOM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_mom_slope ,sim_status.variance_mom_has_converged, sim_status.variance_mom_has_saturated);
            break;
        case SimulationType::iTEBD:
            log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}" , S_slope,sim_status.entanglement_has_converged, sim_status.entanglement_has_saturated);
            log->info("σ² HAM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_ham_slope ,sim_status.variance_ham_has_converged, sim_status.variance_ham_has_saturated);
            log->info("σ² MOM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_mom_slope ,sim_status.variance_mom_has_converged, sim_status.variance_mom_has_saturated);
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}" , S_slope,sim_status.entanglement_has_converged, sim_status.entanglement_has_saturated);
    log->info("Time                  = {:<16.16f}" , t_tot.get_age());
    log->info("Peak memory           = {:<6.1f} MB" , process_memory_in_mb("VmPeak"));
    t_prt.toc();
}

