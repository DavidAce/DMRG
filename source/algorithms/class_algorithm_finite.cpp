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
    state        = std::make_unique<class_finite_state>();
    state_backup = std::make_unique<class_finite_state>();
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

    // Do a backup
    *state_backup = *state;

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
    state_backup->unset_measurements();
    print_status_update();

    double variance_candidate = tools::finite::measure::energy_variance_per_site(*state);
    double variance_champion  = tools::finite::measure::energy_variance_per_site(*state_backup);
    log->trace("Variance candidate = {}", std::log10(variance_candidate));
    log->trace("Variance champion  = {}", std::log10(variance_champion));
    if (variance_champion < variance_candidate){
        log->trace("Replacing the current state with the champion");
        *state = *state_backup;
    }else{
        log->trace("The current state is better than the champion");
        *state_backup = *state;
    }
//    state->do_all_measurements();

    write_measurements(true);
    write_state(true);
    write_status(true);
    write_logs(true);

//    tools::finite::io::write_all_measurements(*state, *h5pp_file, sim_name);
//    tools::finite::io::write_all_state(*state,*h5pp_file, sim_name);
    tools::finite::io::write_projection_to_closest_parity_sector(*state, *h5pp_file, sim_name,
                                                                 settings::model::target_parity_sector);

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
    state->clear_cache();
//    log->debug("Variance check before move               : {:.16f}", std::log10(tools::finite::measure::energy_variance_per_site(*state)));
    try{
        for(size_t i = 0; i < move_steps;i++){
            tools::finite::mps::move_center_point(*state);
//            log->debug("Variance check after move  {:2}          : {:.16f}",i, std::log10(tools::finite::measure::energy_variance_per_site(*state)));
        }
    }catch(std::exception & e){
        tools::finite::print::print_state(*state);
        throw std::runtime_error("Failed to move center point: " + std::string(e.what()));
    }
}

void class_algorithm_finite::update_bond_dimension(){
    sim_status.chi_max = chi_max();
    if(not chi_grow() or sim_status.bond_dimension_has_reached_max or sim_status.chi_temp == chi_max() ){
        sim_status.chi_temp = chi_max();
        sim_status.bond_dimension_has_reached_max = true;
    }
    if(not sim_status.simulation_has_converged
       and sim_status.simulation_has_saturated
       and sim_status.chi_temp < chi_max()){
        log->trace("Updating bond dimension");
        sim_status.chi_temp = std::min(chi_max(), sim_status.chi_temp * 2);
        log->info("New chi = {}", sim_status.chi_temp);
        clear_saturation_status();
    }
    if(sim_status.chi_temp == chi_max()){
        sim_status.bond_dimension_has_reached_max = true;
    }
    state->set_chi_max(sim_status.chi_max);
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

void class_algorithm_finite::backup_best_state(const class_finite_state &state) {
    log->trace("Checking if given state can beat the backup");
    double variance_candidate  = tools::finite::measure::energy_variance_per_site(state);
    double variance_champion   = tools::finite::measure::energy_variance_per_site(*state_backup);
    if (variance_candidate <  variance_champion){
        log->debug("We have a new champion!");
        *state_backup = state;
    }else{
        log->trace("Champion defended his title");
    }
//    if (sim_status.variance_mpo_has_saturated and sim_status.entanglement_has_saturated){
//
//    }else{
//        log->trace("Waiting for saturation before challenging the current champion");
//    }

}




void class_algorithm_finite::check_convergence_variance(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    if (not state->position_is_any_edge()){return;}
    log->debug("Checking convergence of variance mpo");
    threshold       = std::isnan(threshold) ? settings::precision::varianceConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::varianceSlopeThreshold : slope_threshold;
    auto report = check_saturation_using_slope2(
                    V_mpo_vec,
                    X_mpo_vec,
                    tools::finite::measure::energy_variance_per_site(*state_backup), // HEY! THIS IS NEW
                    sim_status.iteration,
                    1,
                    slope_threshold);
    sim_status.variance_mpo_has_converged = tools::finite::measure::energy_variance_per_site(*state_backup) < threshold;
    if (report.has_computed){
        auto last_nonconverged_ptr = std::find_if(V_mpo_vec.begin(),V_mpo_vec.end(), [threshold](auto const& val){ return val > threshold; });
        if (last_nonconverged_ptr  == V_mpo_vec.end()) last_nonconverged_ptr =  V_mpo_vec.begin();
        size_t converged_count = (size_t)  std::count_if(last_nonconverged_ptr, V_mpo_vec.end(),[threshold](auto const& val){ return val <= threshold; });
        sim_status.variance_mpo_has_saturated = report.has_saturated or converged_count >= min_saturation_iters;
        sim_status.variance_mpo_saturated_for = std::max(converged_count, report.saturated_for) ;
        V_mpo_slope  = report.slopes.back();
        log->debug("Variance slope details:");
        log->debug(" -- relative slope    = {} %", report.slopes.back());
        log->debug(" -- tolerance         = {} %", slope_threshold);
        log->debug(" -- last var average  = {} " , report.avgY.back());
        log->debug(" -- var history       = {} " , V_mpo_vec);
        log->debug(" -- slope history     = {} " , report.slopes);
        log->debug(" -- has saturated     = {} " , sim_status.variance_mpo_has_saturated);
        log->debug(" -- has saturated for = {} " , sim_status.variance_mpo_saturated_for);
        log->debug(" -- converged count   = {}", converged_count);
        if (V_mpo_vec.back() < threshold and sim_status.variance_mpo_saturated_for == 0) throw std::logic_error("Variance should have saturated");
        if (V_mpo_vec.back() < threshold and not sim_status.variance_mpo_has_converged ) throw std::logic_error("Variance should have converged");
    }

//    auto report = check_saturation_using_slope(
//            B_mpo_vec,
//            V_mpo_vec,
//            X_mpo_vec,
//            tools::finite::measure::energy_variance_per_site(*state_backup), // HEY! THIS IS NEW
//            sim_status.iteration,
//            1,
//            slope_threshold);
//
//    if (report.has_computed){
//        sim_status.variance_mpo_has_converged =  V_mpo_vec.back() < threshold;
//        auto last_nonconverged_ptr = std::find_if(V_mpo_vec.begin(),V_mpo_vec.end(), [threshold](auto const& val){ return val > threshold; });
//        if (last_nonconverged_ptr  == V_mpo_vec.end()) last_nonconverged_ptr =  V_mpo_vec.begin();
//        size_t converged_count = (size_t)  std::count_if(last_nonconverged_ptr, V_mpo_vec.end(),[threshold](auto const& val){ return val <= threshold; });
//        size_t saturated_count = (size_t)  std::count(B_mpo_vec.begin(), B_mpo_vec.end(), true);
//        sim_status.variance_mpo_has_saturated = report.has_saturated or converged_count >= min_saturation_iters;
//        sim_status.variance_mpo_saturated_for = std::max(converged_count, saturated_count) ;
//        V_mpo_slope  = report.slope;
//        log->debug("Variance slope details:");
//        log->debug(" -- relative slope    = {} %", report.slope);
//        log->debug(" -- tolerance         = {} %", slope_threshold);
//        log->debug(" -- average           = {} " , report.avgY);
//        log->debug(" -- history           = {} " , V_mpo_vec);
//        log->debug(" -- has saturated     = {} " , sim_status.variance_mpo_has_saturated);
//        log->debug(" -- has saturated for = {} B.size() = {} " , sim_status.variance_mpo_saturated_for,B_mpo_vec.size());
//        log->debug(" -- checked from      = {} to {}", report.check_from, X_mpo_vec.size());
//        log->debug(" -- converged count   = {}", converged_count);
//        log->debug(" -- saturated count   = {}", saturated_count);
//        if (V_mpo_vec.back() < threshold and sim_status.variance_mpo_saturated_for == 0) throw std::logic_error("Variance should have saturated");
//        if (V_mpo_vec.back() < threshold and not sim_status.variance_mpo_has_converged ) throw std::logic_error("Variance should have converged");
//    }

}


void class_algorithm_finite::check_convergence_entg_entropy(double slope_threshold) {
    //Based on the the slope of entanglement entanglement_entropy_midchain
    // This one is cheap to compute.
    if (not state->position_is_any_edge()){return;}
    log->debug("Checking convergence of entanglement");

    slope_threshold = std::isnan(slope_threshold) ? settings::precision::entropySlopeThreshold : slope_threshold;
    auto entropies  = tools::finite::measure::entanglement_entropies(*state_backup);
    std::vector<SaturationReport2> reports(entropies.size());

    for (size_t site = 0; site < entropies.size(); site++){
        reports[site] = check_saturation_using_slope2(
                S_mat[site],
                XS_mat[site],
                entropies[site],
                sim_status.iteration,
                1,
                slope_threshold);
    }
    bool all_computed = std::all_of(reports.begin(), reports.end(), [](const SaturationReport2 r) { return r.has_computed; });
    sim_status.entanglement_has_saturated = false;
    if(all_computed){
        // idx_max_slope is the index to the site with maximum slope
        size_t idx_max_slope = std::distance(reports.begin(),
                                             std::max_element(reports.begin(),reports.end(),
                                   [](const SaturationReport2 &r1, const SaturationReport2 &r2)
                                   {return r1.slopes.back() < r2.slopes.back();}));
        // idx_max_slope is the index to the site with maximum slope
        size_t idx_min_satur = std::distance(reports.begin(),
                                             std::min_element(reports.begin(),reports.end(),
                                   [](const SaturationReport2 &r1, const SaturationReport2 &r2)
                                   {return r1.saturated_for < r2.saturated_for;}));
        S_slope = reports[idx_max_slope].slopes.back();
        sim_status.entanglement_has_saturated = reports[idx_max_slope].has_saturated;
        sim_status.entanglement_saturated_for = reports[idx_min_satur].saturated_for;
        std::vector<double> all_avergs;
        std::vector<double> all_slopes;
        for (auto &r : reports) all_avergs.push_back(r.avgY.back());
        for (auto &r : reports) all_slopes.push_back(r.slopes.back());
        log->debug("Max slope of entanglement entropy at site {}: {:.8f} %", idx_max_slope, S_slope);
        log->debug("Entanglement slope details of worst slope:");
        log->debug(" -- site              = {}"  , idx_max_slope);
        log->debug(" -- relative slope    = {} %", reports[idx_max_slope].slopes.back());
        log->debug(" -- tolerance         = {} %", slope_threshold);
        log->debug(" -- ent history       = {} " , S_mat[idx_max_slope]);
        log->debug(" -- slope history     = {} " , reports[idx_max_slope].slopes);
        log->debug(" -- has saturated     = {} " , reports[idx_max_slope].has_saturated);
        log->debug(" -- has saturated for = {} (site {} )" , sim_status.entanglement_saturated_for, idx_min_satur);
        log->debug(" -- all averages      = {} " , all_avergs);
        log->debug(" -- all slopes        = {} " , all_slopes);
//        for(auto&r:reports) log->debug(" avgY : {} " , r.avgY);
//        for(auto&r:reports) log->debug(" slope: {} " , r.slopes);
        //        if(reports[idx_max_slope].slopes.back() == 0 ) throw std::runtime_error("Max slope is zero! Impossible!");
//        if(idx_max_slope == 0 ) throw std::runtime_error("Site 0 has the worst slope! That's impossible!!");
//        if(idx_max_slope ==  entropies.size() - 1) throw std::runtime_error("Last site has the worst slope! That's impossible!!");
        if (reports[idx_max_slope].slopes.back() > slope_threshold and sim_status.entanglement_has_saturated)
            throw std::logic_error("Not supposed to be saturated!!");
    }
    sim_status.entanglement_has_converged = sim_status.entanglement_has_saturated;

//    for (size_t site = 0; site < entropies.size(); site++){
//        reports[site] = check_saturation_using_slope(
//                BS_mat[site],
//                S_mat[site],
//                XS_mat[site],
//                entropies[site],
//                sim_status.iteration,
//                1,
//                slope_threshold);
//    }
//    bool all_computed = std::all_of(reports.begin(), reports.end(), [](const SaturationReport r) { return r.has_computed; });
//    sim_status.entanglement_has_saturated = false;
//    if(all_computed){
//        // idx is the index to the site with maximum slope
//        size_t idx = std::distance(reports.begin(),
//                     std::max_element(reports.begin(),reports.end(),
//                        [](const SaturationReport &r1, const SaturationReport &r2) {return r1.slope < r2.slope;}));
//        S_slope = reports[idx].slope;
//        sim_status.entanglement_has_saturated = reports[idx].has_saturated;
//        int shortest_saturation_length =  (int) count(BS_mat[idx].begin(), BS_mat[idx].end(), true);
//        for(auto & bool_list : BS_mat){
//            shortest_saturation_length = std::min(shortest_saturation_length, (int) count(bool_list.begin(), bool_list.end(), true));
//        }
//        sim_status.entanglement_saturated_for = shortest_saturation_length;
//        std::vector<double> all_slopes;
//        for (auto &r : reports) all_slopes.push_back(r.slope);
//        log->debug("Max slope of entanglement entropy at site {}: {:.8f} %",idx, S_slope);
//        log->debug("Entanglement slope details of worst slope:");
//        log->debug(" -- site              = {}"  , idx);
//        log->debug(" -- relative slope    = {} %", reports[idx].slope);
//        log->debug(" -- tolerance         = {} %", slope_threshold);
//        log->debug(" -- average           = {} " , reports[idx].avgY);
//        log->debug(" -- history           = {} " , S_mat[idx]);
//        log->debug(" -- has saturated     = {} " , reports[idx].has_saturated);
//        log->debug(" -- has saturated for = {} " , sim_status.entanglement_saturated_for);
//        log->debug(" -- checked from      = {} to {}", reports[idx].check_from, S_mat[idx].size());
//        log->debug(" -- all slopes        = {}", all_slopes);
//    }
//    sim_status.entanglement_has_converged = sim_status.entanglement_has_saturated;

}


void class_algorithm_finite::clear_saturation_status(){
    log->trace("Clearing saturation status");
    for(auto &mat : S_mat){mat.clear();}
    for(auto &mat : BS_mat){mat.clear();}
    for(auto &mat : XS_mat){mat.clear();}

    B_mpo_vec.clear();
    V_mpo_vec.clear();
    X_mpo_vec.clear();

    sim_status.entanglement_has_converged     = false;
    sim_status.entanglement_has_saturated     = false;
    sim_status.entanglement_saturated_for     = 0;

    sim_status.variance_mpo_has_converged     = false;
    sim_status.variance_mpo_has_saturated     = false;
    sim_status.variance_mpo_saturated_for     = 0;


    sim_status.bond_dimension_has_reached_max = false;
    sim_status.simulation_has_to_stop         = false;
    sim_status.simulation_has_converged       = false;
    sim_status.simulation_has_saturated       = false;
    sim_status.simulation_has_succeeded       = false;

}


void class_algorithm_finite::compute_observables(){
    log->trace("Starting all measurements on current mps");
    state->do_all_measurements();
}


void class_algorithm_finite::write_measurements(bool force){
    if (settings::output::storage_level == StorageLevel::NONE){return;}
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

    // THIS IS NOT NEEDED ANYMORE SINCE WE STORE INTO TABLES INSTEAD
    //    if (settings::output::storage_level >= StorageLevel::NORMAL){
    //        std::string log_name = sim_name + "/logs/iter_" + std::to_string(sim_status.iteration);
    //        tools::finite::io::write_all_measurements(*state, *h5pp_file, log_name);
    //    }
    h5pp_file->writeDataset(true, sim_name + "/simOK");
}

void class_algorithm_finite::write_state(bool force){
    if (settings::output::storage_level == StorageLevel::NONE){return;}
    if(not force){
        if (math::mod(sim_status.iteration, write_freq()) != 0) {return;}
        if (not state->position_is_any_edge()){return;}
        if (write_freq() == 0){return;}
        if (settings::output::storage_level <= StorageLevel::NONE){return;}
    }
    log->trace("Writing state to file");
    h5pp_file->writeDataset(false, sim_name + "/simOK");
    tools::finite::io::write_projection_to_closest_parity_sector(*state, *h5pp_file, sim_name,
                                                                 settings::model::target_parity_sector);
    //  Write the wavefunction (this is only defined for short enough state ( L < 14 say)
    if(store_wave_function()){
        h5pp_file->writeDataset(tools::finite::measure::mps_wavefn(*state), sim_name + "/state/psi");
    }
    tools::finite::io::write_all_state(*state, *h5pp_file, sim_name);

    if (settings::output::storage_level >= StorageLevel::FULL){
        log->trace("Writing state to logs");
        std::string log_name = sim_name + "/logs/iter_" + std::to_string(sim_status.iteration);
        tools::finite::io::write_all_state(*state, *h5pp_file, log_name);
    }

    h5pp_file->writeDataset(true, sim_name + "/simOK");
}


void class_algorithm_finite::write_status(bool force){
    if (settings::output::storage_level == StorageLevel::NONE){return;}
    if (not force){
        if (math::mod(sim_status.iteration, write_freq()) != 0) {return;}
        if (not state->position_is_any_edge()){return;}
        if (write_freq() == 0){return;}
        if (settings::output::storage_level <= StorageLevel::LIGHT){return;}
    }
    log->trace("Writing simulation status to file");
    h5pp_file->writeDataset(false, sim_name + "/simOK");
    tools::common::io::write_simulation_status(sim_status, *h5pp_file, sim_name);

    // THIS IS NOT NEEDED ANYMORE SINCE WE STORE INTO TABLES INSTEAD
    // Write the simulation status here as well, since the base has no notion of state edge
    //    if (settings::output::storage_level >= StorageLevel::NORMAL){
    //        log->trace("Writing simulation status to logs");
    //        std::string log_name = sim_name + "/logs/iter_" + std::to_string(sim_status.iteration);
    //        tools::common::io::write_simulation_status(sim_status, *h5pp_file, log_name);
    //    }
    h5pp_file->writeDataset(true, sim_name + "/simOK");
}


void class_algorithm_finite::write_logs(bool force){
    if (settings::output::storage_level == StorageLevel::NONE){return;}
    if(not force){
        if (not settings::output::save_logs){return;}
        if (math::mod(sim_status.step, write_freq()) != 0) {return;}
        if (settings::output::storage_level <= StorageLevel::LIGHT){return;}
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
    report << fmt::format("{:<} "                                             ,sim_name);
    report << fmt::format("iter: {:<4} "                                      ,sim_status.iteration);
    report << fmt::format("step: {:<5} "                                      ,sim_status.step);
    report << fmt::format("L: {} l: {:<2} "                                   ,state->get_length(), state->get_position());
    report << fmt::format("E/L: {:<20.16f} "                                  ,tools::finite::measure::energy_per_site(*state));
    if (sim_type == SimulationType::xDMRG){
        report << fmt::format("ε: {:<6.4f} " ,sim_status.energy_dens);
    }
    report << fmt::format("Sₑ(l): {:<10.8f} "                                 ,tools::finite::measure::entanglement_entropy_current(*state));
    report << fmt::format("log₁₀ σ²(E)/L: {:<10.6f} [{:<10.6f}] "             ,std::log10(tools::finite::measure::energy_variance_per_site(*state)), std::log10(tools::finite::measure::energy_variance_per_site(*state_backup)));
    report << fmt::format("χmax: {:<3} χ: {:<3} "                             ,chi_max(), tools::finite::measure::bond_dimension_current(*state));
    report << fmt::format("log₁₀ trunc: {:<6.4f} "                            ,std::log10(state->truncation_error[state->get_position()]));
    report << fmt::format("con: [σ² {:<5} Sₑ {:<5}] "                         ,sim_status.variance_mpo_has_converged,sim_status.entanglement_has_converged);
    report << fmt::format("sat: [σ² {:<2} Sₑ {:<2}] "                         ,sim_status.variance_mpo_saturated_for,sim_status.entanglement_saturated_for);
    report << fmt::format("time: {:<8.2f}s "                                  ,t_tot.get_age());
    report << fmt::format("mem MB: [Rss {:<.1f} Peak {:<.1f} Vm {:<.1f}] "    ,process_memory_in_mb("VmRSS"), process_memory_in_mb("VmHWM") ,process_memory_in_mb("VmPeak"));
    log->info(report.str());
    t_prt.toc();
}



void class_algorithm_finite::print_status_full(){
    t_prt.tic();
    log->info("{:=^60}","");
    log->info("= {: ^56} =","Final results [" + sim_name + "]");
    log->info("{:=^60}","");

    log->info("--- Final results  --- {} ---", sim_name);
    log->info("Sites                              = {}"    , state->get_length());
    log->info("Iterations                         = {}"    , sim_status.iteration);
    log->info("Sweeps                             = {}"    , state->get_sweeps());
    log->info("Simulation time                    = {:<.1f} s = {:<.2f} min" , t_tot.get_age(), t_tot.get_age()/60);
    log->info("Energy per site E/L                = {:<.16f}"   , tools::finite::measure::energy_per_site(*state));
    if (sim_type == SimulationType::xDMRG){
    log->info("Energy density (rescaled 0 to 1) ε = {:<6.4f}" ,sim_status.energy_dens);
    }
    log->info("Variance per site log₁₀ σ²(E)/L    = {:<.16f}"   , std::log10(tools::finite::measure::energy_variance_per_site(*state)));
    log->info("Bond dimension maximum χmax        = {}"         , chi_max());
    log->info("Bond dimensions χ                  = {}"         , tools::finite::measure::bond_dimensions(*state));
    log->info("Entanglement entropies Sₑ          = {}"         , tools::finite::measure::entanglement_entropies(*state));
    log->info("Truncation Errors                  = {}"         , state->truncation_error);
    log->info("Simulation converged               = {:<}"       , sim_status.simulation_has_converged);
    log->info("Simulation saturated               = {:<}"       , sim_status.simulation_has_saturated);
    log->info("Simulation succeeded               = {:<}"       , sim_status.simulation_has_succeeded);
    log->info("Simulation got stuck               = {:<}"       , sim_status.simulation_has_got_stuck);
    log->info("σ² slope                           = {:<8.4f} %   Converged : {:<8}  Saturated: {:<8}" , V_mpo_slope ,sim_status.variance_mpo_has_converged, sim_status.variance_mpo_has_saturated);
    log->info("Sₑ slope                           = {:<8.4f} %   Converged : {:<8}  Saturated: {:<8}" , S_slope     ,sim_status.entanglement_has_converged, sim_status.entanglement_has_saturated);
    log->info("Memory RSS                         = {:<.1f} MB" , process_memory_in_mb("VmRSS"));
    log->info("Memory Peak                        = {:<.1f} MB" , process_memory_in_mb("VmHWM"));
    log->info("Memory Vm                          = {:<.1f} MB" , process_memory_in_mb("VmPeak"));
    t_prt.toc();
}

