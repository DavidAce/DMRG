//
// Created by david on 2019-06-24.
//


#include "class_algorithm_finite.h"
#include <h5pp/h5pp.h>
#include <mps_tools/nmspc_mps_tools.h>
#include <mps_state/class_finite_chain_state.h>
#include <io/class_hdf5_table_buffer2.h>
#include <general/nmspc_math.h>
#include <iomanip>

class_algorithm_finite::class_algorithm_finite(std::shared_ptr<h5pp::File> h5ppFile_, std::string sim_name, SimulationType sim_type, size_t num_sites)
    : class_algorithm_base(std::move(h5ppFile_), sim_name,sim_type)

{
    log->trace("Constructing class_algorithm_finite");
    state = std::make_unique<class_finite_chain_state>();
    mpstools::finite::mpo::initialize(*state, num_sites, settings::model::model_type);
    mpstools::finite::mps::initialize(*state, num_sites);
    mpstools::finite::mpo::randomize(*state);
    mpstools::finite::mps::randomize(*state);
    mpstools::finite::debug::check_integrity(*state);
    mpstools::finite::mps::project_to_closest_parity(*state, settings::model::symmetry);
    mpstools::finite::debug::check_integrity(*state);

    table_dmrg       = std::make_unique<class_hdf5_table<class_table_dmrg>>        (h5pp_file, sim_name + "/measurements", "simulation_progress", sim_name);
    table_dmrg_chain = std::make_unique<class_hdf5_table<class_table_finite_chain>>(h5pp_file, sim_name + "/measurements", "simulation_progress_full_chain", sim_name);

    min_saturation_length = 1;
    max_saturation_length = 4;
    mpstools::finite::print::print_hamiltonians(*state);
    mpstools::finite::print::print_state(*state);

}


void class_algorithm_finite::run()
/*!
 * \brief Dispatches finite DMRG stages.
 * This function manages the stages of simulation differently depending on whether
 * the data already existed in hdf5 storage or not.
 *
 * There can be two main scenarios that split into cases:
 * 1) The hdf5 file existed already and contains
 *      a) nothing recognizeable (previous crash?)       -- run full simulation from scratch.
 *      b) a converged simulation but no MPS             -- run full simulation from scratch.
 *      c) a not-yet-converged MPS                       -- resume simulation, reset the number of sweeps first.
 *      d) a converged MPS                               -- not much to do... run postprocessing
 * 2) The hdf5 file did not exist                        -- run full simulation from scratch.

 *
 */
{
    if (not sim_on()) { return; }
    log->info("Starting {}",sim_name);
    t_tot.tic();
    if (h5pp_file->getCreateMode() == h5pp::CreateMode::OPEN){
        // This is case 1
        bool finOK_exists = h5pp_file->linkExists("common/finOK");
        bool simOK_exists = h5pp_file->linkExists(sim_name + "/simOK");
        bool mps_exists   = h5pp_file->linkExists(sim_name + "/state/mps");
        bool finOK = false;
        bool simOK = false;
        if(finOK_exists) finOK = h5pp_file->readDataset<bool>("common/finOK");
        if(simOK_exists) simOK = h5pp_file->readDataset<bool>(sim_name + "/simOK");


        if (not simOK or not finOK){
            //Case 1 a -- run full simulation from scratch.
            log->trace("Case 1a");
            run_preprocessing();
            run_simulation();
        }else if(simOK and not mps_exists){
            // Case 1 b
            log->trace("Case 1b");
            run_preprocessing();
            run_simulation();
        }else if(simOK and mps_exists){
            // We can go ahead and load the state from hdf5
            log->trace("Loading MPS from file");
            try{
                mpstools::finite::io::load_from_hdf5(*h5pp_file, *state, sim_state, sim_name);
            }
            catch(std::exception &ex){
                log->error("Failed to load from hdf5: {}", ex.what());
                throw std::runtime_error("Failed to resume from file: " + std::string(ex.what()));
            }
            catch(...){log->error("Unknown error when trying to resume from file.");}

            bool convergence_was_reached  = h5pp_file->readDataset<bool>(sim_name + "/sim_state/simulation_has_converged");
            if(not convergence_was_reached){
                // Case 1 c -- resume simulation, reset the number of sweeps first.
                log->trace("Case 1c");
                settings::xdmrg::max_sweeps += state->get_sweeps();
                run_simulation();

            }else {
                // Case 1 d -- not much else to do.. redo postprocessing for good measure.
                log->trace("Case 1d");
            }
        }
    }else {
        // This is case 2
        log->trace("Case 2");
        run_preprocessing();
        run_simulation();
    }
    t_tot.toc();
    run_postprocessing();
}


void class_algorithm_finite::run_preprocessing(){


}

void class_algorithm_finite::run_postprocessing(){
    log->info("Running {} postprocessing",sim_name);
    mpstools::finite::debug::check_integrity(*state);
    state->unset_measurements();
    state->do_all_measurements();
    mpstools::finite::io::write_all_measurements(*state, *h5pp_file, sim_name);
    mpstools::finite::io::write_all_state(*state,*h5pp_file, sim_name);
    mpstools::finite::io::write_closest_parity_projection(*state, *h5pp_file, sim_name, settings::model::symmetry);

    //  Write the wavefunction (this is only defined for short enough state ( L < 14 say)
    if(store_wave_function()){
        h5pp_file->writeDataset(mpstools::finite::measure::mps_wavefn(*state), sim_name + "/state/psi");
    }
    print_status_full();
    print_profiling();
    h5pp_file->writeDataset(true, sim_name + "/simOK");
    log->info("Finished {} postprocessing",sim_name);
}



void class_algorithm_finite::single_DMRG_step(std::string ritz){
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */
    log->trace("Starting single DMRG step");
    t_sim.tic();
    t_opt.tic();
    Eigen::Tensor<Scalar,4> theta;
    std::tie(theta, sim_state.energy_now) = mpstools::finite::opt::find_optimal_ground_state(*state,sim_state, ritz);
    t_opt.toc();
    t_svd.tic();
    mpstools::finite::opt::truncate_theta(theta, *state, sim_state.chi_temp, settings::precision::SVDThreshold);
    t_svd.toc();
    state->unset_measurements();
    t_sim.toc();
    sim_state.wall_time = t_tot.get_age();
    sim_state.simu_time = t_sim.get_age();
}


void class_algorithm_finite::move_center_point(){
    log->trace("Moving center point ");
    t_sim.tic();
    t_ste.tic();
    size_t move_steps = state->active_sites.empty() ? 1 : std::max(1ul,state->active_sites.size()-2ul);
    try{
        for(size_t i = 0; i < move_steps;i++){
            mpstools::finite::mps::move_center_point(*state);
            sim_state.step += move_steps;
        }
    }catch(std::exception & e){
        mpstools::finite::print::print_state(*state);
        throw std::runtime_error("Failed to move center point: " + std::string(e.what()));
    }
    t_ste.toc();
    t_sim.toc();
}

void class_algorithm_finite::reset_to_random_state(const std::string parity) {
    log->trace("Resetting MPS to random product state");
    if (state->get_length() != (size_t)num_sites()) throw std::range_error("System size mismatch");
    // Randomize state
    mpstools::finite::mps::randomize(*state);
    mpstools::finite::mps::project_to_closest_parity(*state,parity);
    clear_saturation_status();
    sim_state.iteration = state->reset_sweeps();
}


void class_algorithm_finite::check_convergence_variance(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    log->debug("Checking convergence of variance mpo");
    threshold       = std::isnan(threshold)       ? settings::precision::VarConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::VarSaturationThreshold  : slope_threshold;
    compute_observables();
    check_saturation_using_slope(B_mpo_vec,
                                 V_mpo_vec,
                                 X_mpo_vec,
                                 state->measurements.energy_variance_per_site.value(),
                                 sim_state.iteration,
                                 1,
                                 slope_threshold,
                                 V_mpo_slope,
                                 sim_state.variance_mpo_has_saturated);
//    int rewind = std::min((int)B_mpo_vec.size(), (int)measurement->get_chain_length()/4 );
//    auto b_it = B_mpo_vec.begin();
//    std::advance(b_it, -rewind);
    sim_state.variance_mpo_saturated_for = (int) count(B_mpo_vec.begin(), B_mpo_vec.end(), true);
    sim_state.variance_mpo_has_converged =  state->measurements.energy_variance_per_site.value() < threshold;

}


void class_algorithm_finite::check_convergence_entg_entropy(double slope_threshold) {
    //Based on the the slope of entanglement middle_entanglement_entropy
    // This one is cheap to compute.
    log->debug("Checking convergence of entanglement");

    slope_threshold = std::isnan(slope_threshold) ? settings::precision::EntEntrSaturationThreshold  : slope_threshold;
    double entropysum = 0;
    for(auto &entropy : mpstools::finite::measure::entanglement_entropies(*state)){entropysum += entropy;}
    entropysum /= state->get_length();
    check_saturation_using_slope(BS_vec,
                                 S_vec,
                                 XS_vec,
                                 entropysum,
                                 sim_state.step,
                                 1,
                                 slope_threshold,
                                 S_slope,
                                 sim_state.entanglement_has_saturated);
    sim_state.entanglement_has_converged = sim_state.entanglement_has_saturated;
}


void class_algorithm_finite::clear_saturation_status(){
    log->trace("Clearing saturation status");

    BS_vec.clear();
    S_vec.clear();
    XS_vec.clear();

    B_mpo_vec.clear();
    V_mpo_vec.clear();
    X_mpo_vec.clear();

    sim_state.entanglement_has_saturated      = false;
    sim_state.variance_mpo_has_saturated      = false;

    sim_state.variance_mpo_saturated_for = 0;

    sim_state.entanglement_has_converged = false;
    sim_state.variance_mpo_has_converged = false;

    sim_state.bond_dimension_has_reached_max = false;
    sim_state.simulation_has_to_stop         = false;
}


void class_algorithm_finite::compute_observables(){
    log->trace("Starting all measurements on current mps");
    t_sim.tic();
    t_obs.tic();
    state->do_all_measurements();
    t_obs.toc();
    t_sim.toc();
}


void class_algorithm_finite::store_state_and_measurements_to_file(bool force){
    if(not force){
        if (not settings::hdf5::save_progress){return;}
        if (Math::mod(sim_state.iteration, store_freq()) != 0) {return;}
        if (not state->position_is_any_edge()) {return;}
        if (settings::xdmrg::store_freq == 0){return;}
        if (settings::hdf5::storage_level <= StorageLevel::NONE){return;}
    }
    state->unset_measurements();
    compute_observables();
    log->trace("Writing all measurements to file");
    t_sto.tic();
    h5pp_file->writeDataset(false, sim_name + "/simOK");
    mpstools::finite::io::write_all_measurements(*state, *h5pp_file, sim_name);
    mpstools::finite::io::write_closest_parity_projection(*state, *h5pp_file, sim_name, settings::model::symmetry);
    //  Write the wavefunction (this is only defined for short enough state ( L < 14 say)
    if(store_wave_function()){
        h5pp_file->writeDataset(mpstools::finite::measure::mps_wavefn(*state), sim_name + "/state/psi");
    }
    mpstools::finite::io::write_all_state(*state, *h5pp_file, sim_name);
    t_sto.toc();
    store_algorithm_state_to_file();
    h5pp_file->writeDataset(true, sim_name + "/simOK");
}



void class_algorithm_finite::store_table_entry_site_state(bool force){
    if (not force) {
        if (Math::mod(sim_state.iteration, store_freq()) != 0) { return; }
        if (store_freq() == 0) { return; }
//        if (not state->position_is_the_middle()) { return; }
        if (settings::hdf5::storage_level < StorageLevel::FULL){return;}
    }
    log->trace("Storing chain_entry to file");
    compute_observables();
    t_sto.tic();
    table_dmrg_chain->append_record(
            sim_state.iteration,
            state->get_length(),
            state->get_position(),
            state->measurements.bond_dimension_midchain.value(),
            sim_state.energy_now,
            state->measurements.entanglement_entropy_midchain.value(),
            state->truncation_error[state->get_position()]
    );
    t_sto.toc();
}
//
//void class_algorithm_finite::store_table_entry_progress(bool force){
//    if(not force) {
//        if (Math::mod(sim_state.step, store_freq()) != 0) { return; }
//        if (store_freq()) { return; }
//        if (settings::hdf5::storage_level < StorageLevel::NORMAL){return;}
//
//    }
//    compute_observables();
//    log->trace("Storing table entry to file");
//    t_sto.tic();
////    table_dmrg->append_record(
////            sim_state.step,
////            state->get_length(),
////            state->get_position(),
////            state->measurements.bond_dimension_current.value(),
////            settings::xdmrg::chi_max,
////            state->measurements.energy_per_site.value(),
////            sim_state.energy_min,
////            sim_state.energy_max,
////            sim_state.energy_target,
////            state->measurements.energy_variance_per_site.value(),
////            state->measurements.entanglement_entropy_midchain.value(),
////            state->truncation_error[state->get_position()],
////            t_tot.get_age());
//    t_sto.toc();
//}


void class_algorithm_finite::print_profiling(){
    if (settings::profiling::on) {
        log->trace("Printing profiling information (tot)");
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_ste.print_time_w_percent(t_tot);
        t_prt.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        print_profiling_sim(t_sim);
        mpstools::common::profiling::obs::print_profiling(t_obs);
    }
}

void class_algorithm_finite::print_profiling_sim(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        log->trace("Printing profiling information (sim)");
        std::cout << "\n Simulation breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_opt.print_time_w_percent(t_parent);
        t_eig.print_time_w_percent(t_parent);
        t_ham.print_time_w_percent(t_parent);
        t_svd.print_time_w_percent(t_parent);
        t_env.print_time_w_percent(t_parent);
        t_mps.print_time_w_percent(t_parent);
        t_con.print_time_w_percent(t_parent);
    }
}


void class_algorithm_finite::print_status_update() {
    if (Math::mod(sim_state.iteration, print_freq()) != 0) {return;}
//    if (not state->position_is_the_middle()) {return;}
    if (print_freq() == 0) {return;}
    using namespace std;
    using namespace mpstools::finite::measure;
    compute_observables();
    t_prt.tic();
    std::stringstream report;
    report << setprecision(16) << fixed << left;
    report << left  << sim_name << " ";
    report << left  << "Iter: "                       << setw(6) << sim_state.iteration;
    report << left  << "E: ";

    switch(sim_type) {
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            report << setw(21) << setprecision(16)    << fixed   << energy_per_site(*state);
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }

    if (sim_type == SimulationType::xDMRG){
        report << left  << " ε: "<< setw(8) << setprecision(4) << fixed << sim_state.energy_dens;
    }

    report << left  << "log₁₀ σ²(E): ";
    switch(sim_type) {
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            report << setw(18) << setprecision(10)    << fixed   << std::log10(energy_variance_per_site(*state));
            break;
        default: throw std::runtime_error("Wrong simulation type");

    }


    report << left  << "S: "                          << setw(21) << setprecision(16)    << fixed   << entanglement_entropy_current(*state);
    report << left  << "χmax: "                       << setw(4)  << setprecision(3)     << fixed   << chi_max();
    report << left  << "χ: "                          << setw(4)  << setprecision(3)     << fixed   << bond_dimension_current(*state);
    report << left  << "log₁₀ trunc: "                << setw(10) << setprecision(4)     << fixed   << std::log10(state->truncation_error[state->get_position()]);
    report << left  << "Sites: "                      << setw(6)  << setprecision(1)     << fixed   << state->get_length();
    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            report << left  << "@ site: "                    << setw(5)  << state->get_position();
            break;
        default: throw std::runtime_error("Wrong simulation type");

    }

    report << left  << " Convergence [";
    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            report << left  << " σ²-"  << std::boolalpha << setw(6) << sim_state.variance_mpo_has_converged;
            report << left  << " S-"   << std::boolalpha << setw(6) << sim_state.entanglement_has_converged;
            break;
        default: throw std::runtime_error("Wrong simulation type");

    }
    report << left  << "]";
    report << left  << " Saturation [";
    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            report << left  << " σ²-" << setw(2) << sim_state.variance_mpo_saturated_for << " steps";
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

void class_algorithm_finite::print_status_full(){
    compute_observables();

    using namespace mpstools::common::measure;
    t_prt.tic();
    log->info("--- Final results  --- {} ---", sim_name);
    log->info("Iterations            = {:<16d}"    , sim_state.iteration);
    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            log->info("Energy MPO            = {:<16.16f}" , state->measurements.energy_per_site.value());
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            log->info("log₁₀ σ²(E) MPO       = {:<16.16f}" , log10(state->measurements.energy_variance_per_site.value()));
            break;
        default: throw std::runtime_error("Wrong simulation type");

    }
    log->info("Mid-chain properties");
    log->info("Entanglement Entropy  = {:<16.16f}" , state->measurements.entanglement_entropy_midchain.value());
    log->info("χmax                  = {:<16d}"    , chi_max());
    log->info("χ                     = {:<16d}"    , state->measurements.bond_dimension_midchain.value());
    log->info("log₁₀ truncation:     = {:<16.16f}" , log10(state->truncation_error[state->get_length()/2]));

    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            log->info("state length          = {:<16d}"    , state->measurements.length.value());
            log->info("Sweep                 = {:<16d}"    , state->get_sweeps());
            break;
        default: throw std::runtime_error("Wrong simulation type");

    }

    log->info("Simulation converged  = {:<}"    , sim_state.simulation_has_converged);

    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            log->info("σ² MPO slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_mpo_slope ,sim_state.variance_mpo_has_converged, sim_state.variance_mpo_has_saturated);
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}" , S_slope,sim_state.entanglement_has_converged, sim_state.entanglement_has_saturated);
    log->info("Time                  = {:<16.16f}" , t_tot.get_age());
    log->info("Peak memory           = {:<6.1f} MB" , process_memory_in_mb("VmPeak"));
    t_prt.toc();
}

