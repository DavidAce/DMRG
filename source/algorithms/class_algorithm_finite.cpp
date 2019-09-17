//
// Created by david on 2019-06-24.
//

#include <iomanip>
#include "class_algorithm_finite.h"
#include <state/class_finite_state.h>
#include <io/class_hdf5_log_buffer.h>
#include <math/nmspc_math.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <general/nmspc_random_numbers.h>
#include <h5pp/h5pp.h>
#include <tools/nmspc_tools.h>

class_algorithm_finite::class_algorithm_finite(std::shared_ptr<h5pp::File> h5ppFile_, std::string sim_name, SimulationType sim_type, size_t num_sites)
    : class_algorithm_base(std::move(h5ppFile_), sim_name,sim_type)

{
    log->trace("Constructing class_algorithm_finite");
    state = std::make_unique<class_finite_state>();
    log->trace("Constructing log buffers in finite base");
    log_measurements       = std::make_shared<class_hdf5_log<class_log_finite_dmrg_measurements>> (h5pp_file, sim_name + "/logs", "measurements", sim_name);


    state->set_chi_max(sim_status.chi_max);
    tools::finite::mpo::initialize(*state, num_sites, settings::model::model_type);
    tools::finite::mps::initialize(*state, num_sites);
    tools::finite::mpo::randomize(*state,settings::model::seed_model);
    tools::finite::mps::randomize(*state,settings::model::initial_parity_sector,settings::model::seed_state);
    tools::finite::debug::check_integrity(*state);


    S_mat.resize(state->get_length()+1);
    BS_mat.resize(state->get_length()+1);
    XS_mat.resize(state->get_length()+1);

    tools::finite::print::print_hamiltonians(*state);
    tools::finite::print::print_state(*state);
    tools::finite::io::write_model(*state, *h5pp_file, sim_name);
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
            // We can go ahead and load the state from output
            log->trace("Loading MPS from file");
            try{
                tools::finite::io::load_from_hdf5(*h5pp_file, *state, sim_status, sim_name);
            }
            catch(std::exception &ex){
                log->error("Failed to load from output: {}", ex.what());
                throw std::runtime_error("Failed to resume from file: " + std::string(ex.what()));
            }
            catch(...){log->error("Unknown error when trying to resume from file.");}

            bool convergence_was_reached  = h5pp_file->readDataset<bool>(sim_name + "/sim_status/simulation_has_converged");
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
    log->info("Running {} preprocessing",sim_name);
    t_pre.tic();
    sim_status.chi_max = chi_max();
    state->set_chi_max(sim_status.chi_max);
    t_pre.toc();
    log->info("Finished {} preprocessing", sim_name);
}




void class_algorithm_finite::single_DMRG_step(std::string ritz){
/*!
 * \fn void single_DMRG_step(std::string ritz)
 */
    log->trace("Starting single xDMRG step");
    t_run.tic();
    Eigen::Tensor<Scalar,4> theta = tools::finite::opt::find_ground_state(*state, ritz);
    tools::finite::opt::truncate_theta(theta, *state, sim_status.chi_temp, settings::precision::SVDThreshold);
    move_center_point();
    state->unset_measurements();
    t_run.toc();
    sim_status.wall_time = t_tot.get_age();
    sim_status.simu_time = t_run.get_measured_time();
}

void class_algorithm_finite::run_postprocessing(){

    log->info("Running {} postprocessing",sim_name);
    t_pos.tic();
    tools::finite::debug::check_integrity(*state);
    state->unset_measurements();
    state->do_all_measurements();
    print_status_update();
    tools::finite::io::write_all_measurements(*state, *h5pp_file, sim_name);
    tools::finite::io::write_all_state(*state,*h5pp_file, sim_name);
    tools::finite::io::write_projection_to_closest_parity_sector(*state, *h5pp_file, sim_name,
                                                                 settings::model::target_parity_sector,false);

    //  Write the wavefunction (this is only defined for short enough state ( L < 14 say)
    if(store_wave_function()){
        h5pp_file->writeDataset(tools::finite::measure::mps_wavefn(*state), sim_name + "/state/psi");
    }
    print_status_full();
    print_profiling();
    h5pp_file->writeDataset(true, sim_name + "/simOK");
    t_pos.toc();
    log->info("Finished {} postprocessing",sim_name);
}


void class_algorithm_finite::move_center_point(){
    log->trace("Moving center point ");
    size_t move_steps = state->active_sites.empty() ? 1 : std::max(1ul,state->active_sites.size()-2ul);
    try{
        for(size_t i = 0; i < move_steps;i++){
            tools::finite::mps::move_center_point(*state);
        }
    }catch(std::exception & e){
        tools::finite::print::print_state(*state);
        throw std::runtime_error("Failed to move center point: " + std::string(e.what()));
    }
}

void class_algorithm_finite::reset_to_random_state(const std::string parity_sector, int seed_state) {
    log->trace("Resetting MPS to random product state in parity sector: {} with seed {}", parity_sector,seed_state);
    if (state->get_length() != (size_t)num_sites()) throw std::range_error("System size mismatch");
    // Randomize state
    state->set_chi_max(chi_max());
    tools::finite::mps::randomize(*state,parity_sector,seed_state, settings::model::use_pauli_eigvecs, settings::model::use_seed_state_as_enumeration);
//    tools::finite::mps::project_to_closest_parity_sector(*state, parity_sector);
    clear_saturation_status();
    sim_status.iteration = state->reset_sweeps();
}


void class_algorithm_finite::check_convergence_variance(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    if (not state->position_is_any_edge()){return;}
    log->debug("Checking convergence of variance mpo");
    threshold       = std::isnan(threshold)       ? settings::precision::VarConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::VarSaturationThreshold  : slope_threshold;
    auto report = check_saturation_using_slope(
                    B_mpo_vec,
                    V_mpo_vec,
                    X_mpo_vec,
                    tools::finite::measure::energy_variance_per_site(*state),
                    sim_status.iteration,
                    1,
                    slope_threshold);
    if (report.has_computed){
        V_mpo_slope  = report.slope;
        log->debug("Variance slope details:");
        log->debug(" -- relative slope  = {} %", report.slope);
        log->debug(" -- tolerance       = {} %", slope_threshold);
        log->debug(" -- average         = {} " , report.avgY);
        log->debug(" -- history         = {} " , V_mpo_vec);
        log->debug(" -- has saturated   = {} " , report.has_saturated);
        log->debug(" -- checked from    = {} to {}", report.check_from, X_mpo_vec.size());
    }
    sim_status.variance_mpo_has_saturated = report.has_saturated;
    sim_status.variance_mpo_saturated_for = (int) count(B_mpo_vec.begin(), B_mpo_vec.end(), true);
    sim_status.variance_mpo_has_converged =  state->measurements.energy_variance_per_site.value() < threshold;

}


void class_algorithm_finite::check_convergence_entg_entropy(double slope_threshold) {
    //Based on the the slope of entanglement entanglement_entropy_midchain
    // This one is cheap to compute.
    if (not state->position_is_any_edge()){return;}
    log->debug("Checking convergence of entanglement");

    slope_threshold = std::isnan(slope_threshold) ? settings::precision::EntEntrSaturationThreshold  : slope_threshold;
    auto entropies  = tools::finite::measure::entanglement_entropies(*state);
    std::vector<SaturationReport> reports(entropies.size());
    for (size_t site = 0; site < entropies.size(); site++){
        reports[site] = check_saturation_using_slope(
                BS_mat[site],
                S_mat[site],
                XS_mat[site],
                entropies[site],
                sim_status.iteration,
                1,
                slope_threshold);
    }
    bool all_computed = std::all_of(reports.begin(), reports.end(), [](const SaturationReport r) { return r.has_computed; });
    sim_status.entanglement_has_saturated = false;
    if(all_computed){
        size_t idx = std::distance(reports.begin(),
                     std::max_element(reports.begin(),reports.end(),
                        [](const SaturationReport &r1, const SaturationReport &r2) {return r1.slope < r2.slope;}));
        S_slope = reports[idx].slope;
        sim_status.entanglement_has_saturated = reports[idx].has_saturated;
        sim_status.entanglement_saturated_for = (int) count(BS_mat[idx].begin(), BS_mat[idx].end(), true);
        std::vector<double> all_slopes;
        for (auto &r : reports) all_slopes.push_back(r.slope);
        log->debug("Max slope of entanglement entropy at site {}: {:.8f} %",idx, S_slope);
        log->debug("Entanglement slope details of worst slope:");
        log->debug(" -- site            = {}"  , idx);
        log->debug(" -- relative slope  = {} %", reports[idx].slope);
        log->debug(" -- tolerance       = {} %", slope_threshold);
        log->debug(" -- average         = {} " , reports[idx].avgY);
        log->debug(" -- history         = {} " , S_mat[idx]);
        log->debug(" -- has saturated   = {} " , reports[idx].has_saturated);
        log->debug(" -- checked from    = {} to {}", reports[idx].check_from, S_mat[idx].size());
        log->debug(" -- all slopes      = {}", all_slopes);
    }
    sim_status.entanglement_has_converged = sim_status.entanglement_has_saturated;

}


void class_algorithm_finite::clear_saturation_status(){
    log->trace("Clearing saturation status");
    for(auto &mat : S_mat){mat.clear();}
    for(auto &mat : BS_mat){mat.clear();}
    for(auto &mat : XS_mat){mat.clear();}

    B_mpo_vec.clear();
    V_mpo_vec.clear();
    X_mpo_vec.clear();

    sim_status.entanglement_has_saturated      = false;
    sim_status.variance_mpo_has_saturated      = false;
    sim_status.variance_mpo_saturated_for      = 0;

    sim_status.entanglement_has_converged     = false;
    sim_status.variance_mpo_has_converged     = false;

    sim_status.bond_dimension_has_reached_max = false;
    sim_status.simulation_has_to_stop         = false;
}


void class_algorithm_finite::compute_observables(){
    log->trace("Starting all measurements on current mps");
    state->do_all_measurements();
}


void class_algorithm_finite::write_measurements(bool force){
    if(not force){
        if (math::mod(sim_status.iteration, write_freq()) != 0) {return;}
        if (not state->position_is_any_edge()){return;}
        if (write_freq() == 0){return;}
        if (settings::output::storage_level <= StorageLevel::NONE){return;}
    }
    log->trace("Writing all measurements to file");
    state->unset_measurements();
    compute_observables();
    h5pp_file->writeDataset(false, sim_name + "/simOK");
    tools::finite::io::write_all_measurements(*state, *h5pp_file, sim_name);

    if (settings::output::storage_level >= StorageLevel::NORMAL){
        std::string log_name = sim_name + "/logs/iter_" + std::to_string(sim_status.iteration);
        tools::finite::io::write_all_measurements(*state, *h5pp_file, log_name);
    }
    h5pp_file->writeDataset(true, sim_name + "/simOK");
}

void class_algorithm_finite::write_state(bool force){
    if(not force){
        if (math::mod(sim_status.iteration, write_freq()) != 0) {return;}
        if (not state->position_is_any_edge()){return;}
        if (write_freq() == 0){return;}
        if (settings::output::storage_level <= StorageLevel::NONE){return;}
    }
    log->trace("Writing state to file");
    h5pp_file->writeDataset(false, sim_name + "/simOK");
    tools::finite::io::write_projection_to_closest_parity_sector(*state, *h5pp_file, sim_name,
                                                                 settings::model::target_parity_sector, false);
    //  Write the wavefunction (this is only defined for short enough state ( L < 14 say)
    if(store_wave_function()){
        h5pp_file->writeDataset(tools::finite::measure::mps_wavefn(*state), sim_name + "/state/psi");
    }
    tools::finite::io::write_all_state(*state, *h5pp_file, sim_name);

    if (settings::output::storage_level >= StorageLevel::FULL){
        std::string log_name = sim_name + "/logs/iter_" + std::to_string(sim_status.iteration);
        tools::finite::io::write_all_state(*state, *h5pp_file, log_name);

    }
    //Write the simulation status here as well, since the base has no notion of state edge
    if (settings::output::storage_level >= StorageLevel::NORMAL){
        std::string log_name = sim_name + "/logs/iter_" + std::to_string(sim_status.iteration);
        tools::common::io::write_simulation_status(sim_status, *h5pp_file, log_name);
    }

    h5pp_file->writeDataset(true, sim_name + "/simOK");
}


void class_algorithm_finite::write_status(bool force){
    if (not force){
        if (math::mod(sim_status.iteration, write_freq()) != 0) {return;}
        if (not state->position_is_any_edge()){return;}
        if (write_freq() == 0){return;}
        if (settings::output::storage_level <= StorageLevel::NONE){return;}
    }
    log->trace("Writing simulation status to file");
    h5pp_file->writeDataset(false, sim_name + "/simOK");
    tools::common::io::write_simulation_status(sim_status, *h5pp_file, sim_name);
    h5pp_file->writeDataset(true, sim_name + "/simOK");
}


void class_algorithm_finite::write_logs(bool force){
    if(not force){
        if (not settings::output::save_logs){return;}
        if (math::mod(sim_status.step, write_freq()) != 0) {return;}
        if (settings::output::storage_level < StorageLevel::NORMAL){return;}
    }
    write_log_measurement();
    write_log_sim_status();
    write_log_profiling();

}

void class_algorithm_finite::write_log_sim_status(){
    if (log_sim_status == nullptr){return;}
    log->trace("Appending sim_status log entry");
    log_sim_status->append_record(sim_status);
}



void class_algorithm_finite::write_log_measurement(){
    if (log_measurements == nullptr){return;}
    log->trace("Appending measurement log entry");
    state->do_all_measurements();
    class_log_finite_dmrg_measurements::data_struct measurements_entry;
    measurements_entry.step                            = sim_status.step;
    measurements_entry.iteration                       = sim_status.iteration;
    measurements_entry.position                        = sim_status.position;
    measurements_entry.length                          = state->get_length();
    measurements_entry.bond_dimension_midchain         = state->measurements.bond_dimension_midchain.value();
    measurements_entry.bond_dimension_current          = state->measurements.bond_dimension_current.value();
    measurements_entry.entanglement_entropy_midchain   = state->measurements.entanglement_entropy_midchain.value();
    measurements_entry.entanglement_entropy_current    = state->measurements.entanglement_entropy_current.value();
    measurements_entry.norm                            = state->measurements.norm.value();
    measurements_entry.energy                          = state->measurements.energy.value();
    measurements_entry.energy_per_site                 = state->measurements.energy_per_site.value();
    measurements_entry.energy_variance                 = state->measurements.energy_variance_mpo.value();
    measurements_entry.energy_variance_per_site        = state->measurements.energy_variance_per_site.value();
    measurements_entry.spin_component_sx               = state->measurements.spin_component_sx.value();
    measurements_entry.spin_component_sy               = state->measurements.spin_component_sy.value();
    measurements_entry.spin_component_sz               = state->measurements.spin_component_sz.value();
    measurements_entry.truncation_error                = state->truncation_error[state->get_position()];
    measurements_entry.wall_time                       = t_tot.get_age();
    log_measurements->append_record(measurements_entry);

}

void class_algorithm_finite::write_log_profiling(){
    if (log_profiling == nullptr){return;}
    log->trace("Appending profiling log entry");
    class_log_profiling::data_struct profiling_entry;
    profiling_entry.step            = sim_status.step;
    profiling_entry.iteration       = sim_status.iteration;
    profiling_entry.position        = sim_status.position;
    profiling_entry.t_tot           = t_tot.get_age();
    profiling_entry.t_run           = t_run.get_measured_time();

    profiling_entry.t_eig           = tools::common::profile::t_eig.get_measured_time();
    profiling_entry.t_svd           = tools::common::profile::t_svd.get_measured_time();
    profiling_entry.t_ene           = tools::common::profile::t_ene.get_measured_time();
    profiling_entry.t_var           = tools::common::profile::t_var.get_measured_time();
    profiling_entry.t_ent           = tools::common::profile::t_ent.get_measured_time();
    profiling_entry.t_hdf           = tools::common::profile::t_hdf.get_measured_time();
    profiling_entry.t_prj           = tools::common::profile::t_prj.get_measured_time();
    profiling_entry.t_opt           = tools::common::profile::t_opt.get_measured_time();
    profiling_entry.t_chk           = tools::common::profile::t_chk.get_measured_time();
    profiling_entry.t_ene_mpo       = tools::common::profile::t_ene_mpo.get_measured_time();
    profiling_entry.t_ene_ham       = tools::common::profile::t_ene_ham.get_measured_time();
    profiling_entry.t_ene_mom       = tools::common::profile::t_ene_mom.get_measured_time();
    profiling_entry.t_var_mpo       = tools::common::profile::t_var_mpo.get_measured_time();
    profiling_entry.t_var_ham       = tools::common::profile::t_var_ham.get_measured_time();
    profiling_entry.t_var_mom       = tools::common::profile::t_var_mom.get_measured_time();


    profiling_entry.t_env = 0;
    profiling_entry.t_evo = 0;
    profiling_entry.t_udt = 0;
    profiling_entry.t_ste = 0;
    profiling_entry.t_prt = 0;
    profiling_entry.t_obs = 0;
    profiling_entry.t_mps = 0;
    profiling_entry.t_chi = 0;

    log_profiling->append_record(profiling_entry);
}



void class_algorithm_finite::print_status_update() {
    if (math::mod(sim_status.step, print_freq()) != 0) {return;}
//    if (not state->position_is_the_middle()) {return;}
    if (print_freq() == 0) {return;}
    using namespace std;
    using namespace tools::finite::measure;
    compute_observables();
    t_prt.tic();
    std::stringstream report;
    report << setprecision(16) << fixed << left;
    report << left  << sim_name << " ";
    report << left  << "Iter: "                       << setw(6) << sim_status.iteration;
    report << left  << "E: ";

    switch(sim_type) {
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            report << setw(21) << setprecision(16)    << fixed   << energy_per_site(*state);
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }

    if (sim_type == SimulationType::xDMRG){
        report << left  << " ε: "<< setw(8) << setprecision(4) << fixed << sim_status.energy_dens;
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
            report << left  << " σ²-"  << std::boolalpha << setw(6) << sim_status.variance_mpo_has_converged;
            report << left  << " S-"   << std::boolalpha << setw(6) << sim_status.entanglement_has_converged;
            break;
        default: throw std::runtime_error("Wrong simulation type");

    }
    report << left  << "]";
    report << left  << " Saturation [";
    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            report << left  << "σ²:" << setw(2) << sim_status.variance_mpo_saturated_for;
            report << left  << " S:" << setw(2) << sim_status.entanglement_saturated_for;
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
    t_prt.tic();
    log->info("--- Final results  --- {} ---", sim_name);
    log->info("Iterations            = {:<16d}"    , sim_status.iteration);
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
    log->info("χmax                    = {:<16d}"    , chi_max());
    log->info("χ                       = {}"         , state->measurements.bond_dimensions.value());
    log->info("Entanglement Entropies  = {}"         , state->measurements.entanglement_entropies.value());
    log->info("Truncation Errors       = {}"         , state->truncation_error);

    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            log->info("state length          = {:<16d}"    , state->measurements.length.value());
            log->info("Sweep                 = {:<16d}"    , state->get_sweeps());
            break;
        default: throw std::runtime_error("Wrong simulation type");

    }

    log->info("Simulation converged  = {:<}"    , sim_status.simulation_has_converged);

    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            log->info("σ² MPO slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_mpo_slope ,sim_status.variance_mpo_has_converged, sim_status.variance_mpo_has_saturated);
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}" , S_slope,sim_status.entanglement_has_converged, sim_status.entanglement_has_saturated);
    log->info("Time                  = {:<16.16f}" , t_tot.get_age());
    log->info("Peak memory           = {:<6.1f} MB" , process_memory_in_mb("VmPeak"));
    t_prt.toc();
}

