//
// Created by david on 2019-06-24.
//

#include "class_algorithm_infinite.h"
#include <config/nmspc_settings.h>
#include <math/num.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/infinite/io.h>
#include <tools/infinite/measure.h>
#include <tools/infinite/mps.h>

class_algorithm_infinite::class_algorithm_infinite(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type)
    : class_algorithm_base(std::move(h5ppFile_), algo_type) {
    tools::log->trace("Constructing algorithm infinite");
    tensors.initialize(settings::model::model_type);
    tensors.state->set_algorithm(algo_type);
}

void class_algorithm_infinite::run() {
    if(not cfg_algorithm_is_on()) return;
    auto t_tot = tools::common::profile::t_tot->tic_token();
    run_preprocessing();
    run_simulation();
    run_postprocessing();
}

void class_algorithm_infinite::run_preprocessing() {
    auto t_pre = tools::common::profile::prof[algo_type]["t_pre"]->tic_token();
    status.clear();
    randomize_model(); // First use of random!
    init_bond_dimension_limits();
    write_to_file(StorageReason::MODEL);
}

void class_algorithm_infinite::run_postprocessing() {
    auto t_pos = tools::common::profile::prof[algo_type]["t_pos"]->tic_token();
    write_to_file(StorageReason::FINISHED);
    copy_from_tmp(StorageReason::FINISHED);
    print_status_full();
    tools::common::profile::print_profiling();
}

void class_algorithm_infinite::randomize_model() {
    tools::log->info("Randomizing model");
    tensors.randomize_model();
    clear_convergence_status();
}

void class_algorithm_infinite::update_bond_dimension_limit(std::optional<long> tmp_bond_limit) {
    if(tmp_bond_limit.has_value()) {
        status.chi_lim = tmp_bond_limit.value();
        return;
    }

    try {
        long chi_lim_now = status.chi_lim;
        if(chi_lim_now < cfg_chi_lim_init()) throw std::logic_error("Chi limit should be larger than chi init");
    } catch(std::exception &error) {
        // If we reached this stage, either
        // 1) chi_lim is not initialized yet
        // 2) chi_lim is initialized, but it is smaller than the init value found in settings
        // Either way, we should set chi_lim to be chi_lim_init, unless cfg_chi_lim_init is larger than tmp_bond_limit
        tools::log->info("Setting initial bond dimension limit: {}", cfg_chi_lim_init());
        status.chi_lim      = cfg_chi_lim_init();
        status.chi_lim_init = cfg_chi_lim_init();
        return;
    }

    status.chi_lim_has_reached_chi_max = status.chi_lim >= cfg_chi_lim_max();
    if(not status.chi_lim_has_reached_chi_max) {
        if(cfg_chi_lim_grow()) {
            // Here the settings specify to grow the bond dimension limit progressively during the simulation
            // Only do this if the simulation is stuck.
            if(status.algorithm_has_stuck_for > 0) {
                tools::log->debug("Truncation error : {}", tensors.state->get_truncation_error());
                tools::log->debug("Bond dimensions  : {}", tensors.state->chiC());
                if(tensors.state->get_truncation_error() > 0.5 * settings::precision::svd_threshold and tensors.state->chiC() >= status.chi_lim) {
                    // Write results before updating bond dimension chi
                    //                    backup_best_state(*state);
                    write_to_file(StorageReason::CHI_UPDATE);
                    long chi_new_limit = std::min(status.chi_lim_max, status.chi_lim * 2);
                    tools::log->info("Updating bond dimension limit {} -> {}", status.chi_lim, chi_new_limit);
                    status.chi_lim = chi_new_limit;
                    clear_convergence_status();
                    status.chi_lim_has_reached_chi_max = status.chi_lim == cfg_chi_lim_max();

                    copy_from_tmp();

                } else {
                    tools::log->debug("cfg_chi_lim_grow is ON, and simulation is stuck, but there is no reason to increase bond dimension -> Kept current bond "
                                      "dimension limit {}",
                                      status.chi_lim);
                }
            } else {
                tools::log->debug("Not stuck -> Kept current bond dimension limit {}", status.chi_lim);
            }
        } else {
            // Here the settings specify to just set the limit to maximum chi directly
            tools::log->info("Setting bond dimension limit to maximum = {}", cfg_chi_lim_max());
            status.chi_lim = cfg_chi_lim_max();
        }
    } else {
        tools::log->debug("Chi limit has reached max: {} -> Kept current bond dimension limit {}", cfg_chi_lim_max(), status.chi_lim);
    }
    status.chi_lim = status.chi_lim;
    if(status.chi_lim > status.chi_lim_max)
        throw std::runtime_error(fmt::format("chi_lim is larger than cfg_chi_lim_max! {} > {}", status.chi_lim, status.chi_lim_max));
}

// void class_algorithm_infinite::update_bond_dimension_limit(std::optional<long> max_bond_dim){
//    if(not max_bond_dim.has_value()) {
//        tools::log->debug("No max bond dim given, setting {}", cfg_chi_lim_max());
//        max_bond_dim = cfg_chi_lim_max();
//    }
//    try{
//        long chi_lim_now = status.chi_lim;
//        if(chi_lim_now < cfg_chi_lim_init())
//            throw std::logic_error("Chi limit should be larger than chi init");
//    }catch(std::exception &error){
//        //If we reached this stage, either
//        // 1) chi_lim is not initialized yet
//        // 2) chi_lim is initialized, but it is smaller than the init value found in settings
//        // Either way, we should set chi_lim to be chi_lim_init, unless cfg_chi_lim_init is larger than max_bond_dim
//        tools::log->info("Setting initial bond dimension limit: {}", cfg_chi_lim_init());
//        tensors.state->set_chi_lim(std::min(max_bond_dim.value(),cfg_chi_lim_init()));
//        status.cfg_chi_lim_max = max_bond_dim.value();
//        status.chi_lim = status.chi_lim;
//        return;
//    }
//
//    status.chi_lim_has_reached_chi_max = status.chi_lim == max_bond_dim;
//    if(not status.chi_lim_has_reached_chi_max){
//        if(cfg_chi_lim_grow()){
//            // Here the settings specify to grow the bond dimension limit progressively during the simulation
//            // Only do this if the simulation is stuck.
//            if(status.algorithm_has_got_stuck){
//                long chi_new_limit = std::min(max_bond_dim.value(), status.chi_lim * 2);
//                tools::log->debug("Updating bond dimension limit {} -> {}", status.chi_lim, chi_new_limit);
//                tensors.state->set_chi_lim(chi_new_limit);
//                clear_convergence_status();
//            }else{
//                tools::log->debug("cfg_chi_lim_grow is ON but sim is not stuck -> Kept current bond dimension limit {}", status.chi_lim);
//            }
//        }else{
//            // Here the settings specify to just set the limit to maximum chi directly
//            tools::log->debug("Setting bond dimension limit to maximum = {}", cfg_chi_lim_max());
//            tensors.state->set_chi_lim(max_bond_dim.value());
//        }
//    }else{
//        tools::log->debug("Chi limit has reached max: {} -> Kept current bond dimension limit {}", cfg_chi_lim_max(),state->get_chi_lim());
//    }
//    status.cfg_chi_lim_max = max_bond_dim.value();
//    status.chi_lim = status.chi_lim;
//}

void class_algorithm_infinite::randomize_state(ResetReason reason, std::optional<std::string> sector, std::optional<long> bitfield,
                                               std::optional<bool> use_eigenspinors) {
    tools::log->trace("Resetting to random product state");
    if(reason == ResetReason::SATURATED) {
        if(status.num_resets >= settings::strategy::max_resets)
            return tools::log->warn("Skipped reset: num resets {} >= max resets {}", status.num_resets, settings::strategy::max_resets);
        else
            status.num_resets++;
    }

    if(not sector) sector = settings::strategy::target_sector;
    if(not bitfield) bitfield = settings::input::bitfield;
    if(not use_eigenspinors) use_eigenspinors = settings::strategy::use_eigenspinors;

    status.iter = 0;
    tensors.reset_to_random_product_state(sector.value(), bitfield.value(), use_eigenspinors.value());
    // Randomize state
    clear_convergence_status();
}


void class_algorithm_infinite::clear_convergence_status() {
    tools::log->trace("Clearing saturation status");

    var_mpo_iter.clear();
    var_ham_iter.clear();
    var_mom_iter.clear();
    entropy_iter.clear();
    status.algorithm_has_finished      = false;
    status.algorithm_has_succeeded     = false;
    status.algorithm_has_to_stop       = false;
    status.algorithm_has_stuck_for     = 0;
    status.algorithm_saturated_for     = 0;
    status.algorithm_converged_for     = 0;
    status.entanglement_converged_for  = 0;
    status.entanglement_saturated_for  = 0;
    status.variance_mpo_converged_for  = 0;
    status.variance_mpo_saturated_for  = 0;
    status.variance_ham_converged_for  = 0;
    status.variance_ham_saturated_for  = 0;
    status.variance_mom_converged_for  = 0;
    status.variance_mom_saturated_for  = 0;
    status.chi_lim_has_reached_chi_max = false;
}

// void class_algorithm_infinite::enlarge_environment() {
//    tools::log->trace("Enlarging environment");
//    tensors.insert_site_pair();
//}
//
// void class_algorithm_infinite::swap() {
//    tools::log->trace("Swap AB sites on state");
//    tensors.state->enlarge();
//}

void class_algorithm_infinite::check_convergence_variance_mpo(std::optional<double> threshold, std::optional<double> sensitivity) {
    tools::log->debug("Checking convergence of variance mpo");
    if (not threshold) threshold = settings::precision::variance_convergence_threshold;
    if (not sensitivity) sensitivity = settings::precision::variance_saturation_sensitivity;
    var_mpo_iter.emplace_back(tools::infinite::measure::energy_variance_per_site_mpo(tensors));
    status.variance_mpo_converged_for = count_convergence(var_mpo_iter,threshold.value());
    auto report = check_saturation(var_mpo_iter, sensitivity.value());
    if(report.has_computed) {
        status.variance_mpo_saturated_for = report.saturated_count;
    }
}

void class_algorithm_infinite::check_convergence_variance_ham(std::optional<double> threshold, std::optional<double> sensitivity) {
    tools::log->trace("Checking convergence of variance ham");
    if (not threshold) threshold = settings::precision::variance_convergence_threshold;
    if (not sensitivity) sensitivity = settings::precision::variance_saturation_sensitivity;
    var_ham_iter.emplace_back(tools::infinite::measure::energy_variance_per_site_ham(tensors));
    status.variance_ham_converged_for = count_convergence(var_ham_iter,threshold.value());
    auto report     = check_saturation(var_ham_iter, sensitivity.value());
    if(report.has_computed) {
        status.variance_ham_saturated_for = report.saturated_count;
    }
}

void class_algorithm_infinite::check_convergence_variance_mom(std::optional<double> threshold, std::optional<double> sensitivity) {
    tools::log->trace("Checking convergence of variance mom");
    if (not threshold) threshold = settings::precision::variance_convergence_threshold;
    if (not sensitivity) sensitivity = settings::precision::variance_saturation_sensitivity;
    var_mom_iter.emplace_back(tools::infinite::measure::energy_variance_per_site_mom(tensors));
    status.variance_mom_converged_for = count_convergence(var_mom_iter,threshold.value());
    auto report     = check_saturation(var_mom_iter, sensitivity.value());
    if(report.has_computed) {
        status.variance_mom_saturated_for = report.saturated_count;
    }
}

void class_algorithm_infinite::check_convergence_entg_entropy(std::optional<double> sensitivity) {
    tools::log->debug("Checking convergence of entanglement");
    if (not sensitivity) sensitivity = settings::precision::entropy_saturation_sensitivity;
    entropy_iter.emplace_back(tools::infinite::measure::entanglement_entropy(*tensors.state));
    auto report     = check_saturation(entropy_iter, sensitivity.value());
    if(report.has_computed) {
        status.entanglement_saturated_for = report.saturated_count;
        status.entanglement_converged_for = report.saturated_count;
    }
}

void class_algorithm_infinite::write_to_file(StorageReason storage_reason, std::optional<CopyPolicy> copy_policy) {
    StorageLevel             storage_level;
    std::string              state_prefix = fmt::format("{}/{}", algo_name, tensors.state->get_name()); // May get modified
    std::string              model_prefix = fmt::format("{}/{}", algo_name, "model");
    std::vector<std::string> table_prefxs = {fmt::format("{}/{}", state_prefix, "tables")}; // Common tables

    switch(storage_reason) {
        case StorageReason::FINISHED: {
            if(status.algorithm_has_succeeded)
                storage_level = settings::output::storage_level_good_state;
            else
                storage_level = settings::output::storage_level_fail_state;
            state_prefix += "/finished";
            table_prefxs.emplace_back(state_prefix); // Appends to its own table as well as the common ones
            break;
        }
        case StorageReason::SAVEPOINT: {
            if(num::mod(status.iter, settings::output::savepoint_frequency) != 0) return;
            state_prefix += "/savepoint";
            storage_level = settings::output::storage_level_savepoint;
            if(settings::output::savepoint_keep_newest_only)
                state_prefix += "/iter_last";
            else
                state_prefix += fmt::format("/iter_{}", status.iter);
            table_prefxs.emplace_back(state_prefix); // Appends to its own table as well as the common ones
            break;
        }
        case StorageReason::CHECKPOINT: {
            if(num::mod(status.iter, settings::output::checkpoint_frequency) != 0) return;
            state_prefix += "/checkpoint";
            storage_level = settings::output::storage_level_checkpoint;
            if(settings::output::checkpoint_keep_newest_only)
                state_prefix += "/iter_last";
            else
                state_prefix += fmt::format("/iter_{}", status.iter);
            table_prefxs.emplace_back(state_prefix); // Appends to its own table as well as the common ones
            break;
        }
        case StorageReason::CHI_UPDATE: {
            if(not cfg_chi_lim_grow()) return;
            storage_level = settings::output::storage_level_checkpoint;
            state_prefix += fmt::format("/checkpoint/chi_{}", status.chi_lim);
            table_prefxs = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::PROJ_STATE: {
            storage_level = settings::output::storage_level_proj_state;
            state_prefix += "/projection";
            table_prefxs = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::INIT_STATE: {
            storage_level = settings::output::storage_level_init_state;
            state_prefix += "/state_init";
            table_prefxs = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::EMIN_STATE: {
            storage_level = settings::output::storage_level_emin_state;
            state_prefix  = algo_name + "/state_emin";
            table_prefxs  = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::EMAX_STATE: {
            storage_level = settings::output::storage_level_emax_state;
            state_prefix  = algo_name + "/state_emax";
            table_prefxs  = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::MODEL: {
            storage_level = settings::output::storage_level_model;
            tools::infinite::io::h5table::save_model(*h5pp_file, model_prefix + "/hamiltonian", storage_level, *tensors.model);
            tools::infinite::io::h5dset::save_model(*h5pp_file, model_prefix + "/mpo", storage_level, *tensors.model);
            copy_from_tmp(storage_reason, CopyPolicy::TRY);
            return;
        }
    }
    if(storage_level == StorageLevel::NONE) return;
    if(state_prefix.empty()) throw std::runtime_error("State prefix is empty");
    tools::log->info("Writing to file: Reason [{}] | Level [{}] | hdf5 prefix [{}]", enum2str(storage_reason), enum2str(storage_level), state_prefix);
    // Start saving tensors and metadata
    tools::infinite::io::h5dset::save_state(*h5pp_file, state_prefix, storage_level, *tensors.state, status);
    tools::infinite::io::h5dset::save_edges(*h5pp_file, state_prefix, storage_level, *tensors.edges);
    tools::common::io::h5attr::save_meta(*h5pp_file, storage_level, storage_reason, settings::model::model_type, settings::model::model_size, algo_type,
                                         tensors.state->get_name(), state_prefix, model_prefix, status);
    // Some storage reasons should not go further. Like projection.
    if(storage_reason == StorageReason::PROJ_STATE) return;

    // The main results have now been written. Next we append data to tables
    for(const auto &table_prefix : table_prefxs) {
        tools::infinite::io::h5table::save_measurements(*h5pp_file, table_prefix + "/measurements", storage_level, tensors, status);
        tools::infinite::io::h5table::save_sim_status(*h5pp_file, table_prefix + "/status", storage_level, status);
        tools::infinite::io::h5table::save_profiling(*h5pp_file, table_prefix + "/profiling", storage_level, status);
        tools::infinite::io::h5table::save_mem_usage(*h5pp_file, table_prefix + "/mem_usage", storage_level, status);
    }
    // Copy from temporary location to destination depending on given policy
    copy_from_tmp(storage_reason, copy_policy);
}

void class_algorithm_infinite::print_status_update() {
    if(num::mod(status.iter, cfg_print_freq()) != 0) { return; }
    if(cfg_print_freq() == 0) { return; }
    //    compute_observables();
    std::string report;
    report += fmt::format("{:<} ", tensors.state->get_name());
    report += fmt::format("iter: {:<4} ", status.iter);
    report += fmt::format("step: {:<5} ", status.step);

    switch(algo_type) {
        case AlgorithmType::iDMRG:
            report += fmt::format("E/L: mpo {:<20.16f} ham {:<20.16f} mom {:<20.16f}", tools::infinite::measure::energy_per_site_mpo(tensors),
                                  tools::infinite::measure::energy_per_site_ham(tensors), tools::infinite::measure::energy_per_site_mom(tensors));
            break;
        case AlgorithmType::iTEBD:
            report += fmt::format("E/L: ham {:<20.16f} mom {:<20.16f}", tools::infinite::measure::energy_per_site_ham(tensors),
                                  tools::infinite::measure::energy_per_site_mom(tensors));
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }

    switch(algo_type) {
        case AlgorithmType::iDMRG:
            report +=
                fmt::format("log₁₀ σ²(E): mpo {:<10.6f} ham {:<10.6f} mom {:<10.6f}", tools::infinite::measure::energy_variance_per_site_mpo(tensors),
                            tools::infinite::measure::energy_variance_per_site_ham(tensors), tools::infinite::measure::energy_variance_per_site_mom(tensors));
            break;
        case AlgorithmType::iTEBD:
            report += fmt::format("log₁₀ σ²(E): ham {:<10.6f} mom {:<10.6f}", tools::infinite::measure::energy_variance_per_site_ham(tensors),
                                  tools::infinite::measure::energy_variance_per_site_mom(tensors));
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    report += fmt::format("Sₑ(l): {:<10.8f} ", tools::infinite::measure::entanglement_entropy(*tensors.state));
    report += fmt::format("χmax: {:<3} χlim: {:<3} χ: {:<3} ", cfg_chi_lim_max(), status.chi_lim, tools::infinite::measure::bond_dimension(*tensors.state));
    report += fmt::format("log₁₀ trunc: {:<10.4f} ", std::log10(tools::infinite::measure::truncation_error(*tensors.state)));
    report += fmt::format("Sites: {:6}", tensors.get_length());

    report += fmt::format("stk: {:<1} ", status.algorithm_has_stuck_for);
    report += fmt::format("sat: [σ² {:<1} Sₑ {:<1}] ", status.variance_mpo_saturated_for, status.entanglement_saturated_for);
    switch(algo_type) {
        case AlgorithmType::iDMRG:
            report += fmt::format("sat: [σ² {:<1} Sₑ {:<1}] ", status.variance_mpo_saturated_for, status.entanglement_saturated_for);
            break;
        case AlgorithmType::iTEBD: report += fmt::format("sat: [Sₑ {:<1}] ", status.entanglement_saturated_for); break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    report += fmt::format("con: {:<4} ", status.algorithm_converged_for);
    report += fmt::format("time:{:>8.2f}s ", tools::common::profile::t_tot->get_measured_time());
    report += fmt::format("mem: [rss {:<.1f} peak {:<.1f} vm {:<.1f}] MB ", tools::common::profile::mem_rss_in_mb(), tools::common::profile::mem_hwm_in_mb(),
                          tools::common::profile::mem_vm_in_mb());
    tools::log->info(report);
}

void class_algorithm_infinite::print_status_full() {
    //    compute_observables();
    tensors.state->do_all_measurements();
    using namespace std;
    using namespace tools::infinite::measure;
    tools::log->info("--- Final results  --- {} ---", algo_name);
    tools::log->info("Iterations            = {:<16d}", status.iter);
    switch(algo_type) {
        case AlgorithmType::iDMRG:
            tools::log->info("Energy MPO            = {:<16.16f}", tools::infinite::measure::energy_per_site_mpo(tensors));
            tools::log->info("Energy HAM            = {:<16.16f}", tools::infinite::measure::energy_per_site_ham(tensors));
            tools::log->info("Energy MOM            = {:<16.16f}", tools::infinite::measure::energy_per_site_mom(tensors));
            break;
        case AlgorithmType::iTEBD:
            tools::log->info("Energy HAM            = {:<16.16f}", tools::infinite::measure::energy_per_site_ham(tensors));
            tools::log->info("Energy MOM            = {:<16.16f}", tools::infinite::measure::energy_per_site_mom(tensors));
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }
    switch(algo_type) {
        case AlgorithmType::iDMRG:
            tools::log->info("log₁₀ σ²(E) MPO       = {:<16.16f}", std::log10(tools::infinite::measure::energy_variance_per_site_mpo(tensors)));
            tools::log->info("log₁₀ σ²(E) HAM       = {:<16.16f}", std::log10(tools::infinite::measure::energy_variance_per_site_ham(tensors)));
            tools::log->info("log₁₀ σ²(E) MOM       = {:<16.16f}", std::log10(tools::infinite::measure::energy_variance_per_site_mom(tensors)));
            break;
        case AlgorithmType::iTEBD:
            tools::log->info("log₁₀ σ²(E) HAM       = {:<16.16f}", std::log10(tools::infinite::measure::energy_variance_per_site_ham(tensors)));
            tools::log->info("log₁₀ σ²(E) MOM       = {:<16.16f}", std::log10(tools::infinite::measure::energy_variance_per_site_mom(tensors)));
            break;
        default: throw std::runtime_error("Wrong simulation type");
    }

    tools::log->info("Entanglement Entropy  = {:<16.16f}", tools::infinite::measure::entanglement_entropy(*tensors.state));
    tools::log->info("χmax                  = {:<16d}", cfg_chi_lim_max());
    tools::log->info("χ                     = {:<16d}", tools::infinite::measure::bond_dimension(*tensors.state));
    tools::log->info("log₁₀ truncation:     = {:<16.16f}", log10(std::log10(tools::infinite::measure::truncation_error(*tensors.state))));

    switch(algo_type) {
        case AlgorithmType::iDMRG: break;
        case AlgorithmType::iTEBD: tools::log->info("δt                    = {:<16.16f}", status.delta_t); break;

        default: throw std::runtime_error("Wrong simulation type");
    }

    tools::log->info("Algorithm succeeded      = {:<}", status.algorithm_has_succeeded);
    tools::log->info("Algorithm saturated for  = {:<}", status.algorithm_saturated_for);
    tools::log->info("Algorithm converged for  = {:<}", status.algorithm_converged_for);
    tools::log->info("Algorithm has stuck for  = {:<}", status.algorithm_has_stuck_for);
    tools::log->info("σ² MPO                   = Converged : {:<4}  Saturated: {:<4}", status.variance_mpo_converged_for, status.variance_mpo_saturated_for);
    tools::log->info("σ² HAM                   = Converged : {:<4}  Saturated: {:<4}", status.variance_ham_converged_for, status.variance_ham_saturated_for);
    tools::log->info("σ² MOM                   = Converged : {:<4}  Saturated: {:<4}", status.variance_mom_converged_for, status.variance_mom_saturated_for);
    tools::log->info("Sₑ                       = Converged : {:<4}  Saturated: {:<4}", status.entanglement_converged_for, status.entanglement_saturated_for);
    tools::log->info("Time                     = {:<16.16f}", tools::common::profile::t_tot->get_age());

    tools::log->info("Peak memory           = {:<6.1f} MB", tools::common::profile::mem_hwm_in_mb());
}
