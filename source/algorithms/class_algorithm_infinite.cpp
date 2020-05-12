//
// Created by david on 2019-06-24.
//
#include "class_algorithm_infinite.h"
#include <h5pp/h5pp.h>
#include <math/nmspc_math.h>
#include <state/class_state_infinite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/infinite.h>

class_algorithm_infinite::class_algorithm_infinite(std::shared_ptr<h5pp::File> h5ppFile_, SimulationType sim_type)
    : class_algorithm_base(std::move(h5ppFile_), sim_type) {
    state = std::make_unique<class_state_infinite>(sim_type);
    state->set_chi_lim(2); // Can't call chi_init() <-- it's a pure virtual function
    tools::infinite::mpo::initialize(*state, settings::model::model_type);
    tools::infinite::mps::initialize(*state, settings::model::model_type);
    tools::infinite::mpo::randomize(*state);
    tools::infinite::env::initialize(*state);
    tools::infinite::debug::check_integrity(*state);
    state->assert_positions();
}

void class_algorithm_infinite::run() {
    if(not sim_on()) return;
    tools::common::profile::t_tot->tic();
    run_preprocessing();
    run_simulation();
    run_postprocessing();
    tools::common::profile::t_tot->toc();
}

void class_algorithm_infinite::run_old() {
    if(not sim_on()) return;
    tools::common::profile::t_tot->tic();
    run_preprocessing();
    run_simulation();
    run_postprocessing();
    tools::common::profile::t_tot->toc();
}

void class_algorithm_infinite::run_preprocessing() {
    tools::common::profile::t_pre->tic();
    tools::infinite::io::h5table::write_model(*h5pp_file, sim_name + "/model", settings::output::storage_level_model, *state);
    tools::infinite::io::h5dset::write_mpo(*h5pp_file, sim_name + "/model",  settings::output::storage_level_model, *state);
    state->set_chi_max(chi_max());
    sim_status.chi_max = chi_max();
    update_bond_dimension_limit(chi_init());

    tools::common::profile::t_pre->toc();
}

void class_algorithm_infinite::run_postprocessing() {
    tools::common::profile::t_pos->tic();
    write_to_file(StorageReason::FINISHED);
    copy_from_tmp(StorageReason::FINISHED);
    print_status_full();
    tools::common::profile::t_pos->toc();
    tools::common::profile::print_profiling();
}

void class_algorithm_infinite::update_truncation_limit() {
    // Will update SVD threshold iff the state precision is being limited by truncation error
}

void class_algorithm_infinite::update_bond_dimension_limit(std::optional<long> tmp_bond_limit) {
    if(tmp_bond_limit.has_value()) {
        state->set_chi_lim(tmp_bond_limit.value());
        sim_status.chi_lim = tmp_bond_limit.value();
        return;
    }

    try {
        long chi_lim_now = state->get_chi_lim();
        if(chi_lim_now < chi_init()) throw std::logic_error("Chi limit should be larger than chi init");
    } catch(std::exception &error) {
        // If we reached this stage, either
        // 1) chi_lim is not initialized yet
        // 2) chi_lim is initialized, but it is smaller than the init value found in settings
        // Either way, we should set chi_lim to be chi_init, unless chi_init is larger than tmp_bond_limit
        tools::log->info("Setting initial bond dimension limit: {}", chi_init());
        state->set_chi_lim(chi_init());
        sim_status.chi_lim = chi_init();
        return;
    }

    sim_status.chi_lim_has_reached_chi_max = state->get_chi_lim() >= chi_max();
    if(not sim_status.chi_lim_has_reached_chi_max) {
        if(chi_grow()) {
            // Here the settings specify to grow the bond dimension limit progressively during the simulation
            // Only do this if the simulation is stuck.
            if(sim_status.simulation_has_got_stuck) {
                tools::log->debug("Truncation error : {}", state->get_truncation_error());
                tools::log->debug("Bond dimensions  : {}", tools::infinite::measure::bond_dimension(*state));
                if(state->get_truncation_error() > 0.5 * settings::precision::svd_threshold and
                   tools::infinite::measure::bond_dimension(*state) >= state->get_chi_lim()) {
                    // Write results before updating bond dimension chi
                    //                    backup_best_state(*state);
                    write_to_file(StorageReason::CHI_UPDATE);
                    long chi_new_limit = std::min(state->get_chi_max(), state->get_chi_lim() * 2);
                    tools::log->info("Updating bond dimension limit {} -> {}", state->get_chi_lim(), chi_new_limit);
                    state->set_chi_lim(chi_new_limit);
                    clear_saturation_status();
                    sim_status.chi_lim_has_reached_chi_max = state->get_chi_lim() == chi_max();

                    copy_from_tmp();

                } else {
                    tools::log->debug(
                        "chi_grow is ON, and simulation is stuck, but there is no reason to increase bond dimension -> Kept current bond dimension limit {}",
                        state->get_chi_lim());
                }
            } else {
                tools::log->debug("Not stuck -> Kept current bond dimension limit {}", state->get_chi_lim());
            }
        } else {
            // Here the settings specify to just set the limit to maximum chi directly
            tools::log->info("Setting bond dimension limit to maximum = {}", chi_max());
            state->set_chi_lim(chi_max());
        }
    } else {
        tools::log->debug("Chi limit has reached max: {} -> Kept current bond dimension limit {}", chi_max(), state->get_chi_lim());
    }
    sim_status.chi_lim = state->get_chi_lim();
    if(state->get_chi_lim() > state->get_chi_max())
        throw std::runtime_error(fmt::format("chi_lim is larger than chi_max! {} > {}", state->get_chi_lim(), state->get_chi_max()));
}

// void class_algorithm_infinite::update_bond_dimension_limit(std::optional<long> max_bond_dim){
//    if(not max_bond_dim.has_value()) {
//        tools::log->debug("No max bond dim given, setting {}", chi_max());
//        max_bond_dim = chi_max();
//    }
//    try{
//        long chi_lim_now = state->get_chi_lim();
//        if(chi_lim_now < chi_init())
//            throw std::logic_error("Chi limit should be larger than chi init");
//    }catch(std::exception &error){
//        //If we reached this stage, either
//        // 1) chi_lim is not initialized yet
//        // 2) chi_lim is initialized, but it is smaller than the init value found in settings
//        // Either way, we should set chi_lim to be chi_init, unless chi_init is larger than max_bond_dim
//        tools::log->info("Setting initial bond dimension limit: {}", chi_init());
//        state->set_chi_lim(std::min(max_bond_dim.value(),chi_init()));
//        sim_status.chi_max = max_bond_dim.value();
//        sim_status.chi_lim = state->get_chi_lim();
//        return;
//    }
//
//    sim_status.chi_lim_has_reached_chi_max = state->get_chi_lim() == max_bond_dim;
//    if(not sim_status.chi_lim_has_reached_chi_max){
//        if(chi_grow()){
//            // Here the settings specify to grow the bond dimension limit progressively during the simulation
//            // Only do this if the simulation is stuck.
//            if(sim_status.simulation_has_got_stuck){
//                long chi_new_limit = std::min(max_bond_dim.value(), state->get_chi_lim() * 2);
//                tools::log->debug("Updating bond dimension limit {} -> {}", state->get_chi_lim(), chi_new_limit);
//                state->set_chi_lim(chi_new_limit);
//                clear_saturation_status();
//            }else{
//                tools::log->debug("chi_grow is ON but sim is not stuck -> Kept current bond dimension limit {}", state->get_chi_lim());
//            }
//        }else{
//            // Here the settings specify to just set the limit to maximum chi directly
//            tools::log->debug("Setting bond dimension limit to maximum = {}", chi_max());
//            state->set_chi_lim(max_bond_dim.value());
//        }
//    }else{
//        tools::log->debug("Chi limit has reached max: {} -> Kept current bond dimension limit {}", chi_max(),state->get_chi_lim());
//    }
//    sim_status.chi_max = max_bond_dim.value();
//    sim_status.chi_lim = state->get_chi_lim();
//}

void class_algorithm_infinite::reset_to_initial_state() {
    tools::log->trace("Resetting MPS to random product state");
    sim_status.iter = 0;

    // Randomize state
    *state = tools::infinite::mps::set_random_state(*state, settings::strategy::initial_parity_sector);
    clear_saturation_status();
}

void class_algorithm_infinite::reset_to_random_product_state(const std::string &parity) {
    tools::log->trace("Resetting MPS to random product state");
    sim_status.iter = 0;

    // Randomize state
    *state = tools::infinite::mps::set_random_state(*state, parity);
    clear_saturation_status();
}

void class_algorithm_infinite::reset_to_random_current_state([[maybe_unused]] std::optional<double> chi_lim) {
    tools::log->critical("Resetting MPS state based on current is not currently implemented for infinite MPS algorithms");
    throw std::runtime_error("Resetting MPS state based on current is not currently implemented for infinite MPS algorithms");
}

void class_algorithm_infinite::clear_saturation_status() {
    tools::log->trace("Clearing saturation status");

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

    sim_status.entanglement_has_saturated = false;
    sim_status.variance_mpo_has_saturated = false;
    sim_status.variance_ham_has_saturated = false;
    sim_status.variance_mom_has_saturated = false;

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

void class_algorithm_infinite::enlarge_environment() {
    tools::log->trace("Enlarging environment");
    state->enlarge_environment(0);
}

void class_algorithm_infinite::swap() {
    tools::log->trace("Swap AB sites on state");
    state->swap_AB();
}

void class_algorithm_infinite::check_convergence_variance_mpo(double threshold, double slope_threshold) {
    // Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    tools::log->debug("Checking convergence of variance mpo");
    threshold       = std::isnan(threshold) ? settings::precision::variance_convergence_threshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::variance_slope_threshold : slope_threshold;
    //    compute_observables();

    auto report =
        check_saturation_using_slope(V_mpo_vec, X_mpo_vec, tools::infinite::measure::energy_variance_per_site_mpo(*state), sim_status.iter, 1, slope_threshold);
    //    if(report.has_computed) V_mpo_slope  = report.slopes.back(); //TODO: Fix this, changed slope calculation, back is not relevant
    if(report.has_computed) {
        V_mpo_slope = report.slope; // TODO: Fix this, changed slope calculation, back is not relevant
        sim_status.variance_mpo_has_saturated = V_mpo_slope < slope_threshold;
        sim_status.variance_mpo_saturated_for = (size_t) count(B_mpo_vec.begin(), B_mpo_vec.end(), true);
        sim_status.variance_mpo_has_converged = tools::infinite::measure::energy_variance_per_site_mpo(*state) < threshold;
    }
}

void class_algorithm_infinite::check_convergence_variance_ham(double threshold, double slope_threshold) {
    // Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    tools::log->trace("Checking convergence of variance ham");

    threshold       = std::isnan(threshold) ? settings::precision::variance_convergence_threshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::variance_slope_threshold : slope_threshold;
    auto report     = check_saturation_using_slope(
        //            B_ham_vec,
        V_ham_vec, X_ham_vec, tools::infinite::measure::energy_variance_per_site_ham(*state), sim_status.iter, 1, slope_threshold);
    //    if(report.has_computed) V_ham_slope  = report.slopes.back();//TODO: Fix this, changed slope calculation, back is not relevant
    if(report.has_computed) {
        V_ham_slope = report.slope; // TODO: Fix this, changed slope calculation, back is not relevant
        sim_status.variance_ham_has_saturated = V_ham_slope < slope_threshold;
        sim_status.variance_ham_has_converged = tools::infinite::measure::energy_variance_per_site_ham(*state) < threshold;
    }
}

void class_algorithm_infinite::check_convergence_variance_mom(double threshold, double slope_threshold) {
    // Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    tools::log->trace("Checking convergence of variance mom");

    threshold       = std::isnan(threshold) ? settings::precision::variance_convergence_threshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::variance_slope_threshold : slope_threshold;
    auto report     = check_saturation_using_slope(
        //            B_mom_vec,
        V_mom_vec, X_mom_vec, tools::infinite::measure::energy_variance_per_site_mom(*state), sim_status.iter, 1, slope_threshold);
    if(report.has_computed) {
        V_mom_slope = report.slope; // TODO: Fix this, slopes.back() not relevant anymore
        sim_status.variance_mom_has_saturated = V_mom_slope < slope_threshold;
        sim_status.variance_mom_has_converged = tools::infinite::measure::energy_variance_per_site_mom(*state) < threshold;
    }
}

void class_algorithm_infinite::check_convergence_entg_entropy(double slope_threshold) {
    // Based on the the slope of entanglement entanglement_entropy_midchain
    // This one is cheap to compute.
    tools::log->debug("Checking convergence of entanglement");

    slope_threshold = std::isnan(slope_threshold) ? settings::precision::entropy_slope_threshold : slope_threshold;
    auto report     = check_saturation_using_slope(
        //            BS_vec,
        S_vec, XS_vec, tools::infinite::measure::entanglement_entropy(*state), sim_status.iter, 1, slope_threshold);
    if(report.has_computed) {
        S_slope = report.slope; // TODO: Fix this, changed slope calculation, back is not relevant
        sim_status.entanglement_has_saturated = S_slope < slope_threshold;
        sim_status.entanglement_has_converged = sim_status.entanglement_has_saturated;
    }
}

void class_algorithm_infinite::write_to_file(StorageReason storage_reason) {
    StorageLevel      storage_level;
    const std::string table_prefix = sim_name + "/" + state_name;
    std::string       state_prefix = sim_name + "/" + state_name; // May get modified
    std::string       model_prefix = sim_name + "/model";
    switch(storage_reason) {
        case StorageReason::FINISHED: {
            if(sim_status.simulation_has_succeeded)
                storage_level = settings::output::storage_level_good_state;
            else
                storage_level = settings::output::storage_level_fail_state;
            state_prefix.append("/finished");
            break;
        }
        case StorageReason::CHECKPOINT: {
            if(math::mod(sim_status.iter, settings::output::checkpoint_frequency) != 0) return;
            state_prefix.append("/checkpoint");
            storage_level = settings::output::storage_level_checkpoint;
            if(settings::output::checkpoint_keep_newest_only)
                state_prefix.append("/iter_last");
            else
                state_prefix.append("/iter_" + std::to_string(sim_status.iter));
            break;
        }

        case StorageReason::CHI_UPDATE: {
            if(not chi_grow()) return;
            storage_level = settings::output::storage_level_checkpoint;
            state_prefix.append("/checkpoint");
            state_prefix.append("/chi_" + std::to_string(state->get_chi_lim()));
            break;
        }
        case StorageReason::PROJ_STATE: {
            storage_level = settings::output::storage_level_proj_state;
            state_prefix.append("/projection");
            break;
        }
        case StorageReason::INIT_STATE: {
            storage_level = settings::output::storage_level_init_state;
            state_prefix.append("/state_init");
            break;
        }
        case StorageReason::EMIN_STATE: {
            storage_level = settings::output::storage_level_emin_state;
            state_prefix  = sim_name + "/state_emin";
            break;
        }
        case StorageReason::EMAX_STATE: {
            storage_level = settings::output::storage_level_emax_state;
            state_prefix  = sim_name + "/state_emax";
            break;
        }
    }
    if(storage_level == StorageLevel::NONE) return;
    if(state_prefix.empty()) throw std::runtime_error("State prefix is empty");
    tools::infinite::io::h5dset::write_mps(*h5pp_file, state_prefix, storage_level, *state);
    tools::infinite::io::h5dset::write_env(*h5pp_file, state_prefix, storage_level, *state);
    tools::infinite::io::h5dset::write_env2(*h5pp_file, state_prefix, storage_level, *state);
    tools::common::io::h5attr::write_meta(*h5pp_file, sim_name, state_prefix, model_prefix, settings::model::model_type, storage_level, sim_status);


    // The main results have now been written. Next we append data to tables
    // Some storage reasons should not do this however. Like projection.
    // Also we can avoid repeated entries by only allowing fresh step numbers.
    if(storage_reason == StorageReason::PROJ_STATE) return;
    static size_t last_step_written = 0;
    if(sim_status.step == last_step_written and last_step_written > 0) return;


    tools::infinite::io::h5table::write_measurements(*h5pp_file, table_prefix, storage_level, *state, sim_status);
    tools::infinite::io::h5table::write_sim_status(*h5pp_file, table_prefix, storage_level, sim_status);
    tools::infinite::io::h5table::write_profiling(*h5pp_file, table_prefix, storage_level, sim_status);
    tools::infinite::io::h5table::write_mem_usage(*h5pp_file, table_prefix, storage_level, sim_status);
    last_step_written = sim_status.step;
}


void class_algorithm_infinite::copy_from_tmp(StorageReason storage_reason) {
    if(not h5pp_file) return;
    if(not settings::output::use_temp_dir) return;
    switch(storage_reason) {
        case StorageReason::CHECKPOINT:
            if(math::mod(sim_status.iter, settings::output::copy_from_temp_freq) != 0) return; // Check that we write according to the frequency given
        case StorageReason::FINISHED:
        case StorageReason::CHI_UPDATE:
        case StorageReason::PROJ_STATE:
        case StorageReason::INIT_STATE:
        case StorageReason::EMIN_STATE:
        case StorageReason::EMAX_STATE: break;
    }
    tools::common::io::h5tmp::copy_from_tmp(h5pp_file->getFilePath());
}


void class_algorithm_infinite::print_status_update() {
    if(math::mod(sim_status.iter, print_freq()) != 0) {
        return;
    }
    //    if (not state->position_is_the_middle()) {return;}
    if(print_freq() == 0) {
        return;
    }
    //    compute_observables();
    using namespace std;
    std::stringstream report;
    report << setprecision(16) << fixed << left;
    report << left << sim_name << " ";
    report << left << "Iter: " << setw(6) << sim_status.iter;
    report << left << "E: ";

    switch(sim_type) {
        case SimulationType::iDMRG:
            report << setw(21) << setprecision(16) << fixed << tools::infinite::measure::energy_per_site_mpo(*state);
            report << setw(21) << setprecision(16) << fixed << tools::infinite::measure::energy_per_site_ham(*state);
            report << setw(21) << setprecision(16) << fixed << tools::infinite::measure::energy_per_site_mom(*state);
            break;
        case SimulationType::iTEBD:
            report << setw(21) << setprecision(16) << fixed << tools::infinite::measure::energy_per_site_ham(*state);
            report << setw(21) << setprecision(16) << fixed << tools::infinite::measure::energy_per_site_mom(*state);
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }

    report << left << "log₁₀ σ²(E): ";
    switch(sim_type) {
        case SimulationType::iDMRG:
            report << setw(12) << setprecision(4) << fixed << std::log10(tools::infinite::measure::energy_variance_per_site_mpo(*state));
            report << setw(12) << setprecision(4) << fixed << std::log10(tools::infinite::measure::energy_variance_per_site_ham(*state));
            report << setw(12) << setprecision(4) << fixed << std::log10(tools::infinite::measure::energy_variance_per_site_mom(*state));
            break;
        case SimulationType::iTEBD:
            report << setw(12) << setprecision(4) << fixed << std::log10(tools::infinite::measure::energy_variance_per_site_ham(*state));
            report << setw(12) << setprecision(4) << fixed << std::log10(tools::infinite::measure::energy_variance_per_site_mom(*state));
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }

    report << left << "S: " << setw(21) << setprecision(16) << fixed << tools::infinite::measure::entanglement_entropy(*state);
    report << left << "χmax: " << setw(4) << setprecision(3) << fixed << chi_max();
    report << left << "χ: " << setw(4) << setprecision(3) << fixed << tools::infinite::measure::bond_dimension(*state);
    report << left << "log₁₀ trunc: " << setw(10) << setprecision(4) << fixed << std::log10(tools::infinite::measure::truncation_error(*state));
    report << left << "Sites: " << setw(6) << setprecision(1) << fixed << tools::infinite::measure::length(*state);
    switch(sim_type) {
        case SimulationType::iDMRG:
        case SimulationType::iTEBD: break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    report << left << " Convergence [";
    switch(sim_type) {
        case SimulationType::iDMRG:
            report << left << " S-" << std::boolalpha << setw(6) << sim_status.entanglement_has_converged;
            report << left << " σ²-" << std::boolalpha << setw(6) << sim_status.variance_mpo_has_converged;
            break;
        case SimulationType::iTEBD: report << left << " S-" << std::boolalpha << setw(6) << sim_status.entanglement_has_converged; break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    report << left << "]";
    report << left << " Saturation [";
    switch(sim_type) {
        case SimulationType::iDMRG:
            report << left << " σ²- " << setw(2) << sim_status.variance_mpo_saturated_for << " steps";
            report << left << " S-" << std::boolalpha << setw(6) << sim_status.entanglement_has_saturated;
            break;
        case SimulationType::iTEBD: report << left << " S-" << std::boolalpha << setw(6) << sim_status.entanglement_has_saturated; break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    report << left << "]";
    report << left << " Time: " << setw(10) << setprecision(2) << fixed << tools::common::profile::t_tot->get_age();
    //    report << left << " Memory [";
    //    report << left << "Rss: "     << process_memory_in_mb("VmRSS")<< " MB ";
    //    report << left << "RssPeak: "  << process_memory_in_mb("VmHWM")<< " MB ";
    //    report << left << "VmPeak: "  << process_memory_in_mb("VmPeak")<< " MB";
    report << left << "]";
    report << fmt::format("mem MB: [Rss {:<.1f} Peak {:<.1f} Vm {:<.1f}] ", tools::common::profile::mem_rss_in_mb(), tools::common::profile::mem_hwm_in_mb(),
                          tools::common::profile::mem_vm_in_mb());
    tools::log->info(report.str());
}

void class_algorithm_infinite::print_status_full() {
    //    compute_observables();
    state->do_all_measurements();
    using namespace std;
    using namespace tools::infinite::measure;
    tools::log->info("--- Final results  --- {} ---", sim_name);
    tools::log->info("Iterations            = {:<16d}", sim_status.iter);
    switch(sim_type) {
        case SimulationType::iDMRG:
            tools::log->info("Energy MPO            = {:<16.16f}", tools::infinite::measure::energy_per_site_mpo(*state));
            tools::log->info("Energy HAM            = {:<16.16f}", tools::infinite::measure::energy_per_site_ham(*state));
            tools::log->info("Energy MOM            = {:<16.16f}", tools::infinite::measure::energy_per_site_mom(*state));
            break;
        case SimulationType::iTEBD:
            tools::log->info("Energy HAM            = {:<16.16f}", tools::infinite::measure::energy_per_site_ham(*state));
            tools::log->info("Energy MOM            = {:<16.16f}", tools::infinite::measure::energy_per_site_mom(*state));
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    switch(sim_type) {
        case SimulationType::iDMRG:
            tools::log->info("log₁₀ σ²(E) MPO       = {:<16.16f}", std::log10(tools::infinite::measure::energy_variance_per_site_mpo(*state)));
            tools::log->info("log₁₀ σ²(E) HAM       = {:<16.16f}", std::log10(tools::infinite::measure::energy_variance_per_site_ham(*state)));
            tools::log->info("log₁₀ σ²(E) MOM       = {:<16.16f}", std::log10(tools::infinite::measure::energy_variance_per_site_mom(*state)));
            break;
        case SimulationType::iTEBD:
            tools::log->info("log₁₀ σ²(E) HAM       = {:<16.16f}", std::log10(tools::infinite::measure::energy_variance_per_site_ham(*state)));
            tools::log->info("log₁₀ σ²(E) MOM       = {:<16.16f}", std::log10(tools::infinite::measure::energy_variance_per_site_mom(*state)));
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }

    tools::log->info("Entanglement Entropy  = {:<16.16f}", tools::infinite::measure::entanglement_entropy(*state));
    tools::log->info("χmax                  = {:<16d}", chi_max());
    tools::log->info("χ                     = {:<16d}", tools::infinite::measure::bond_dimension(*state));
    tools::log->info("log₁₀ truncation:     = {:<16.16f}", log10(std::log10(tools::infinite::measure::truncation_error(*state))));

    switch(sim_type) {
        case SimulationType::iDMRG: break;
        case SimulationType::iTEBD: tools::log->info("δt                    = {:<16.16f}", sim_status.delta_t); break;

        default: throw std::runtime_error("Wrong simulation type");
    }

    tools::log->info("Simulation saturated  = {:<}", sim_status.simulation_has_saturated);
    tools::log->info("Simulation converged  = {:<}", sim_status.simulation_has_converged);
    tools::log->info("Simulation succeeded  = {:<}", sim_status.simulation_has_succeeded);
    tools::log->info("Simulation got stuck  = {:<}", sim_status.simulation_has_got_stuck);
    switch(sim_type) {
        case SimulationType::iDMRG:
            tools::log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}", S_slope, sim_status.entanglement_has_converged,
                             sim_status.entanglement_has_saturated);
            tools::log->info("σ² MPO slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}", V_mpo_slope, sim_status.variance_mpo_has_converged,
                             sim_status.variance_mpo_has_saturated);
            tools::log->info("σ² HAM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}", V_ham_slope, sim_status.variance_ham_has_converged,
                             sim_status.variance_ham_has_saturated);
            tools::log->info("σ² MOM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}", V_mom_slope, sim_status.variance_mom_has_converged,
                             sim_status.variance_mom_has_saturated);
            break;
        case SimulationType::iTEBD:
            tools::log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}", S_slope, sim_status.entanglement_has_converged,
                             sim_status.entanglement_has_saturated);
            tools::log->info("σ² HAM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}", V_ham_slope, sim_status.variance_ham_has_converged,
                             sim_status.variance_ham_has_saturated);
            tools::log->info("σ² MOM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}", V_mom_slope, sim_status.variance_mom_has_converged,
                             sim_status.variance_mom_has_saturated);
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    tools::log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}", S_slope, sim_status.entanglement_has_converged,
                     sim_status.entanglement_has_saturated);
    tools::log->info("Time                  = {:<16.16f}", tools::common::profile::t_tot->get_age());

    tools::log->info("Peak memory           = {:<6.1f} MB", tools::common::profile::mem_hwm_in_mb());
}
