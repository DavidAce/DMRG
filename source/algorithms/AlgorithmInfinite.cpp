#include "AlgorithmInfinite.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "debug/info.h"
#include "math/num.h"
#include "tensors/model/ModelInfinite.h"
#include "tensors/state/StateInfinite.h"
#include "tid/tid.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include "tools/infinite/h5.h"
#include "tools/infinite/measure.h"

AlgorithmInfinite::AlgorithmInfinite(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type) : AlgorithmBase(std::move(h5ppFile_), algo_type) {
    tools::log->trace("Constructing algorithm infinite");
    tensors.initialize(settings::model::model_type);
    tensors.state->set_algorithm(status.algo_type);
}

void AlgorithmInfinite::run() {
    if(not settings::algorithm_is_on(status.algo_type)) return;
    auto t_tot = tid::get_unscoped("t_tot").tic_token();
    run_preprocessing();
    run_algorithm();
    run_postprocessing();
}

void AlgorithmInfinite::run_preprocessing() {
    status.clear();
    randomize_model(); // First use of random!
    init_bond_dimension_limits();
    write_to_file(StorageReason::MODEL);
}

void AlgorithmInfinite::run_postprocessing() {
    auto t_pos = tid::tic_scope("post");
    write_to_file(StorageReason::FINISHED);
    copy_from_tmp(StorageReason::FINISHED);
    print_status_full();
}

void AlgorithmInfinite::randomize_model() {
    tools::log->info("Randomizing model");
    auto t_rnd = tid::tic_scope("rnd_model");
    tensors.randomize_model();
    clear_convergence_status();
}

void AlgorithmInfinite::update_variance_max_digits(std::optional<double> energy) {
    if(not energy) energy = tools::infinite::measure::energy_per_site_mpo(tensors) * static_cast<double>(status.iter);
    if(settings::precision::use_mpo_energy_shift)
        status.energy_variance_max_digits =
            static_cast<size_t>(std::floor(std::numeric_limits<double>::digits10 - std::max(0.0, std::log10(std::abs(energy.value())))));
    else
        status.energy_variance_max_digits =
            static_cast<size_t>(std::floor(std::numeric_limits<double>::digits10 - std::max(0.0, std::log10(std::pow(energy.value(), 2)))));

    status.energy_variance_prec_limit = std::pow(10, -status.energy_variance_max_digits);
}

void AlgorithmInfinite::update_bond_dimension_limit() {
    status.bond_max                   = settings::get_bond_max(status.algo_type);
    status.bond_limit_has_reached_max = status.bond_lim >= status.bond_max;
    if(settings::strategy::bond_increase_when == UpdateWhen::NEVER) {
        status.bond_lim = status.bond_max;
        return;
    }
    if(status.bond_limit_has_reached_max) return;
    auto tic = tid::tic_scope("bond_grow");

    // If we got here we want to increase the bond dimension limit progressively during the simulation
    bool is_saturated      = status.algorithm_saturated_for > 1; // Allow one round while saturated so that extra efforts get a chance.
    bool is_has_stuck      = status.algorithm_has_stuck_for > 1;
    bool is_truncated      = tensors.state->is_limited_by_bond(status.bond_lim) or tensors.state->is_truncated(status.trnc_lim);
    bool grow_if_truncated = settings::strategy::bond_increase_when == UpdateWhen::TRUNCATED;
    bool grow_if_saturated = settings::strategy::bond_increase_when == UpdateWhen::SATURATED;
    bool grow_if_has_stuck = settings::strategy::bond_increase_when == UpdateWhen::STUCK;

    if(grow_if_truncated and not is_truncated) {
        tools::log->info("State is not limited by its bond dimension. Kept current bond limit {}", status.bond_lim);
        return;
    }
    if(grow_if_saturated and not is_saturated) {
        tools::log->info("Algorithm is not saturated. Kept current bond limit {}", status.bond_lim);
        return;
    }
    if(grow_if_has_stuck and not is_has_stuck) {
        tools::log->info("Algorithm is not stuck. Kept current bond limit {}", status.bond_lim);
        return;
    }

    // If we got to this point we will update the bond dimension by a factor
    auto grow_rate = settings::strategy::bond_increase_rate;
    if(grow_rate <= 1.0) throw except::runtime_error("Error: bond_increase_rate == {:.3f} | must be larger than one", grow_rate);

    // Write current results before updating bond dimension
    write_to_file(StorageReason::BOND_INCREASE);

    // If we got to this point we will update the bond dimension by a factor
    auto factor = settings::strategy::bond_increase_rate;
    if(factor <= 1.0) throw except::logic_error("Error: bond_increase_rate == {:.3f} | must be larger than one", factor);

    auto bond_new = static_cast<double>(status.bond_lim);
    if(grow_rate <= 2.0 and grow_rate > 1.0) {
        bond_new = std::ceil(bond_new * grow_rate);
        bond_new = num::round_up_to_multiple_of<double>(bond_new, 4);
    } else if(grow_rate > 2.0) {
        bond_new = bond_new + grow_rate;
    } else
        throw except::logic_error("Expected grow_rate > 1.0. Got {}", grow_rate);
    bond_new = std::min(bond_new, static_cast<double>(status.bond_max));

    tools::log->info("Updating bond dimension limit {} -> {} | truncated {} | saturated {}", status.bond_lim, bond_new, is_truncated, is_saturated);
    status.bond_lim                   = static_cast<long>(bond_new);
    status.bond_limit_has_reached_max = status.bond_lim == status.bond_max;
    status.algorithm_has_stuck_for    = 0;
    status.algorithm_saturated_for    = 0;
    // Last sanity check before leaving here
    if(status.bond_lim > status.bond_max) throw except::runtime_error("bond_lim is larger than get_bond_max! {} > {}", status.bond_lim, status.bond_max);
}

void AlgorithmInfinite::update_truncation_error_limit() {
    if(status.trnc_lim == 0.0) throw std::runtime_error("trnc_lim is zero!");
    status.trnc_min                   = settings::precision::svd_truncation_lim;
    status.trnc_limit_has_reached_min = status.trnc_lim <= status.trnc_min;
    if(settings::strategy::trnc_decrease_when == UpdateWhen::NEVER or settings::strategy::trnc_decrease_rate == 0.0) {
        status.trnc_lim                   = status.trnc_min;
        status.trnc_limit_has_reached_min = true;
        return;
    }
    if(status.trnc_limit_has_reached_min) return;
    auto tic = tid::tic_scope("trnc_down");

    // If we got here we want to decrease the truncation error limit progressively during the simulation
    bool is_saturated      = status.algorithm_saturated_for > 1; // Allow one round while saturated so that extra efforts get a chance.
    bool is_has_stuck      = status.algorithm_has_stuck_for > 1; // Allow one round while saturated so that extra efforts get a chance.
    bool is_truncated      = tensors.state->is_limited_by_bond(status.bond_lim) or tensors.state->is_truncated(status.trnc_lim);
    bool drop_if_truncated = settings::strategy::trnc_decrease_when == UpdateWhen::TRUNCATED;
    bool drop_if_saturated = settings::strategy::trnc_decrease_when == UpdateWhen::SATURATED;
    bool drop_if_has_stuck = settings::strategy::trnc_decrease_when == UpdateWhen::STUCK;

    if(drop_if_truncated and not is_truncated) {
        tools::log->info("State is not truncated. Kept current truncation error limit {:8.2e}", status.trnc_lim);
        return;
    }
    if(drop_if_saturated and not is_saturated) {
        tools::log->info("Algorithm is not saturated. Kept current truncation error limit {:8.2e}", status.trnc_lim);
        return;
    }
    if(drop_if_has_stuck and not is_has_stuck) {
        tools::log->info("Algorithm is not stuck. Kept current truncation error limit {:8.2e}", status.trnc_lim);
        return;
    }

    // Write current results before updating the truncation error limit
    write_to_file(StorageReason::TRNC_DECREASE);

    // If we got to this point we will update the truncation error limit by a factor
    auto rate = settings::strategy::trnc_decrease_rate;
    if(rate > 1.0 or rate < 0) throw except::runtime_error("Error: trnc_decrease_rate == {:8.2e} | must be in [0, 1]");

    auto trnc_new = std::max(status.trnc_min, status.trnc_lim * rate);

    tools::log->info("Updating truncation error limit {:8.2e} -> {:8.2e} | truncated {} | saturated {}", status.trnc_lim, trnc_new, is_truncated, is_saturated);
    status.trnc_lim                   = trnc_new;
    status.trnc_limit_has_reached_min = status.trnc_lim == status.trnc_min;
    status.algorithm_has_stuck_for    = 0;
    status.algorithm_saturated_for    = 0;

    // Last sanity check before leaving here
    if(status.trnc_lim < status.trnc_min) throw except::logic_error("trnc_lim is smaller than trnc_min ! {:8.2e} > {:8.2e}", status.trnc_lim, status.trnc_min);
}

void AlgorithmInfinite::randomize_state(ResetReason reason, std::optional<std::string> sector, std::optional<bool> use_eigenspinors,
                                        std::optional<size_t> bitfield) {
    tools::log->trace("Resetting to random product state");
    if(reason == ResetReason::SATURATED) {
        if(status.num_resets >= settings::strategy::max_resets)
            return tools::log->warn("Skipped reset: num resets {} >= max resets {}", status.num_resets, settings::strategy::max_resets);
        else
            status.num_resets++;
    }

    if(not sector) sector = settings::strategy::initial_sector;
    if(not bitfield) bitfield = settings::input::bitfield;
    if(not use_eigenspinors) use_eigenspinors = settings::strategy::use_eigenspinors;

    status.iter = 0;
    tensors.reset_to_random_product_state(sector.value(), bitfield.value(), use_eigenspinors.value());
    // Randomize state
    clear_convergence_status();
}

void AlgorithmInfinite::clear_convergence_status() {
    tools::log->trace("Clearing saturation status");
    var_mpo_iter.clear();
    var_ham_iter.clear();
    var_mom_iter.clear();
    entropy_iter.clear();
    status.algo_stop                  = AlgorithmStop::NONE;
    status.algorithm_has_finished     = false;
    status.algorithm_has_succeeded    = false;
    status.algorithm_has_to_stop      = false;
    status.algorithm_has_stuck_for    = 0;
    status.algorithm_saturated_for    = 0;
    status.algorithm_converged_for    = 0;
    status.entanglement_converged_for = 0;
    status.entanglement_saturated_for = 0;
    status.variance_mpo_converged_for = 0;
    status.variance_mpo_saturated_for = 0;
    status.variance_ham_converged_for = 0;
    status.variance_ham_saturated_for = 0;
    status.variance_mom_converged_for = 0;
    status.variance_mom_saturated_for = 0;
    status.bond_limit_has_reached_max = false;
    status.trnc_limit_has_reached_min = false;
}

void AlgorithmInfinite::check_convergence_variance_mpo(std::optional<double> threshold, std::optional<double> sensitivity) {
    tools::log->debug("Checking convergence of variance mpo");
    if(not threshold) threshold = settings::precision::variance_convergence_threshold;
    if(not sensitivity) sensitivity = settings::precision::variance_saturation_sensitivity;
    var_mpo_iter.emplace_back(tools::infinite::measure::energy_variance_per_site_mpo(tensors));
    status.variance_mpo_converged_for = count_convergence(var_mpo_iter, threshold.value());
    auto report                       = check_saturation(var_mpo_iter, sensitivity.value(), SaturationScale::log);
    if(report.has_computed) { status.variance_mpo_saturated_for = report.saturated_count; }
}

void AlgorithmInfinite::check_convergence_variance_ham(std::optional<double> threshold, std::optional<double> sensitivity) {
    tools::log->trace("Checking convergence of variance ham");
    if(not threshold) threshold = settings::precision::variance_convergence_threshold;
    if(not sensitivity) sensitivity = settings::precision::variance_saturation_sensitivity;
    var_ham_iter.emplace_back(tools::infinite::measure::energy_variance_per_site_ham(tensors));
    status.variance_ham_converged_for = count_convergence(var_ham_iter, threshold.value());
    auto report                       = check_saturation(var_ham_iter, sensitivity.value(), SaturationScale::log);
    if(report.has_computed) { status.variance_ham_saturated_for = report.saturated_count; }
}

void AlgorithmInfinite::check_convergence_variance_mom(std::optional<double> threshold, std::optional<double> sensitivity) {
    tools::log->trace("Checking convergence of variance mom");
    if(not threshold) threshold = settings::precision::variance_convergence_threshold;
    if(not sensitivity) sensitivity = settings::precision::variance_saturation_sensitivity;
    var_mom_iter.emplace_back(tools::infinite::measure::energy_variance_per_site_mom(tensors));
    status.variance_mom_converged_for = count_convergence(var_mom_iter, threshold.value());
    auto report                       = check_saturation(var_mom_iter, sensitivity.value(), SaturationScale::log);
    if(report.has_computed) { status.variance_mom_saturated_for = report.saturated_count; }
}

void AlgorithmInfinite::check_convergence_entg_entropy(std::optional<double> sensitivity) {
    tools::log->debug("Checking convergence of entanglement");
    if(not sensitivity) sensitivity = settings::precision::entropy_saturation_sensitivity;
    entropy_iter.emplace_back(tools::infinite::measure::entanglement_entropy(*tensors.state));
    auto report = check_saturation(entropy_iter, sensitivity.value(), SaturationScale::lin);
    if(report.has_computed) {
        status.entanglement_saturated_for = report.saturated_count;
        status.entanglement_converged_for = report.saturated_count;
    }
}

void AlgorithmInfinite::write_to_file(StorageReason storage_reason, std::optional<CopyPolicy> copy_policy) {
    StorageLevel             storage_level;
    std::string              state_prefix = fmt::format("{}/{}", status.algo_type_sv(), tensors.state->get_name()); // May get modified
    std::string              model_prefix = fmt::format("{}/{}", status.algo_type_sv(), "model");
    std::vector<std::string> table_prefxs = {fmt::format("{}/{}", state_prefix, "tables")}; // Common tables

    switch(storage_reason) {
        case StorageReason::FINISHED: {
            storage_level = settings::storage::storage_level_finished;
            state_prefix += "/finished";
            table_prefxs.emplace_back(state_prefix); // Appends to its own table as well as the common ones
            break;
        }
        case StorageReason::SAVEPOINT: {
            if(settings::storage::savepoint_frequency == 0 or num::mod(status.iter, settings::storage::savepoint_frequency) != 0) return;
            state_prefix += "/savepoint";
            storage_level = settings::storage::storage_level_savepoint;
            if(settings::storage::savepoint_keep_newest_only)
                state_prefix += "/iter_last";
            else
                state_prefix += fmt::format("/iter_{}", status.iter);
            table_prefxs.emplace_back(state_prefix); // Appends to its own table as well as the common ones
            break;
        }
        case StorageReason::CHECKPOINT: {
            if(settings::storage::checkpoint_frequency == 0 or num::mod(status.iter, settings::storage::checkpoint_frequency) != 0) return;
            state_prefix += "/checkpoint";
            storage_level = settings::storage::storage_level_checkpoint;
            if(settings::storage::checkpoint_keep_newest_only)
                state_prefix += "/iter_last";
            else
                state_prefix += fmt::format("/iter_{}", status.iter);
            table_prefxs.emplace_back(state_prefix); // Appends to its own table as well as the common ones
            break;
        }
        case StorageReason::BOND_INCREASE: {
            if(settings::strategy::bond_increase_when == UpdateWhen::NEVER) return;
            storage_level = settings::storage::storage_level_bond_state;
            state_prefix += "/bond";
            table_prefxs = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::TRNC_DECREASE: {
            if(settings::strategy::trnc_decrease_when == UpdateWhen::NEVER) return;
            storage_level = settings::storage::storage_level_trnc_state;
            state_prefix += "/trnc";
            table_prefxs = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::FES: {
            if(settings::strategy::fes_rate == 0) return;
            storage_level = settings::storage::storage_level_fes_state;
            state_prefix += "/fes";
            table_prefxs = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::PROJ_STATE: {
            storage_level = settings::storage::storage_level_proj_state;
            state_prefix += "/projection";
            table_prefxs = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::INIT_STATE: {
            storage_level = settings::storage::storage_level_init_state;
            state_prefix += "/state_init";
            table_prefxs = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::EMIN_STATE: {
            storage_level = settings::storage::storage_level_emin_state;
            state_prefix  = fmt::format("{}/state_emin", status.algo_type_sv());
            table_prefxs  = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::EMAX_STATE: {
            storage_level = settings::storage::storage_level_emax_state;
            state_prefix  = fmt::format("{}/state_emax", status.algo_type_sv());
            table_prefxs  = {state_prefix}; // Should not pollute tables other than its own
            break;
        }
        case StorageReason::MODEL: {
            storage_level = settings::storage::storage_level_model;
            tools::infinite::h5::save::model(*h5file, model_prefix + "/hamiltonian", storage_level, *tensors.model);
            tools::infinite::h5::save::mpo(*h5file, model_prefix + "/mpo", storage_level, *tensors.model);
            copy_from_tmp(storage_reason, CopyPolicy::TRY);
            return;
        }
    }
    if(storage_level == StorageLevel::NONE) return;
    if(state_prefix.empty()) throw except::runtime_error("State prefix is empty");
    tools::log->debug("Writing to file: Reason [{}] | Level [{}] | hdf5 prefix [{}]", enum2sv(storage_reason), enum2sv(storage_level), state_prefix);
    // Start saving tensors and metadata
    tools::infinite::h5::save::state(*h5file, state_prefix, storage_level, *tensors.state, status);
    tools::infinite::h5::save::edges(*h5file, state_prefix, storage_level, *tensors.edges);
    tools::common::h5::save::meta(*h5file, storage_level, storage_reason, settings::model::model_type, settings::model::model_size, tensors.state->get_name(),
                                  state_prefix, model_prefix, table_prefxs, status);
    // Some storage reasons should not go further. Like projection.
    if(storage_reason == StorageReason::PROJ_STATE) return;

    // The main results have now been written. Next we append data to tables
    for(const auto &table_prefix : table_prefxs) {
        tools::infinite::h5::save::measurements(*h5file, table_prefix, storage_level, tensors, status);
        tools::common::h5::save::status(*h5file, table_prefix, storage_level, status);
        tools::common::h5::save::timer(*h5file, table_prefix, storage_level, status);
        tools::common::h5::save::mem(*h5file, table_prefix, storage_level, status);
    }
    // Copy from temporary location to destination depending on given policy
    copy_from_tmp(storage_reason, copy_policy);
}

void AlgorithmInfinite::print_status() {
    if(num::mod(status.iter, settings::print_freq(status.algo_type)) != 0) { return; }
    if(settings::print_freq(status.algo_type) == 0) { return; }
    //    compute_observables();
    std::string report;
    report += fmt::format("{:<} ", tensors.state->get_name());
    report += fmt::format("iter: {:<4} ", status.iter);
    report += fmt::format("step: {:<5} ", status.step);

    switch(status.algo_type) {
        case AlgorithmType::iDMRG:
            report += fmt::format("E/L: mpo {:<20.16f} ham {:<20.16f} mom {:<20.16f}", tools::infinite::measure::energy_per_site_mpo(tensors),
                                  tools::infinite::measure::energy_per_site_ham(tensors), tools::infinite::measure::energy_per_site_mom(tensors));
            break;
        case AlgorithmType::iTEBD:
            report += fmt::format("E/L: ham {:<20.16f} mom {:<20.16f}", tools::infinite::measure::energy_per_site_ham(tensors),
                                  tools::infinite::measure::energy_per_site_mom(tensors));
            break;
        default: throw except::runtime_error("Wrong simulation type");
    }

    switch(status.algo_type) {
        case AlgorithmType::iDMRG:
            report +=
                fmt::format("lg σ²H: mpo {:<10.6f} ham {:<10.6f} mom {:<10.6f}", tools::infinite::measure::energy_variance_per_site_mpo(tensors),
                            tools::infinite::measure::energy_variance_per_site_ham(tensors), tools::infinite::measure::energy_variance_per_site_mom(tensors));
            break;
        case AlgorithmType::iTEBD:
            report += fmt::format("lg σ²H: ham {:<10.6f} mom {:<10.6f}", tools::infinite::measure::energy_variance_per_site_ham(tensors),
                                  tools::infinite::measure::energy_variance_per_site_mom(tensors));
            break;
        default: throw except::runtime_error("Wrong simulation type");
    }
    report += fmt::format("ε: {:<8.2e} ", tools::infinite::measure::truncation_error(*tensors.state));
    report += fmt::format("Sₑ(l): {:<10.8f} ", tools::infinite::measure::entanglement_entropy(*tensors.state));
    report += fmt::format("χmax: {:<3} χlim: {:<3} χ: {:<3} ", settings::get_bond_max(status.algo_type), status.bond_lim,
                          tools::infinite::measure::bond_dimension(*tensors.state));
    report += fmt::format("Sites: {:6}", tensors.get_length());

    report += fmt::format("stk: {:<1} ", status.algorithm_has_stuck_for);
    report += fmt::format("sat: [σ² {:<1} Sₑ {:<1}] ", status.variance_mpo_saturated_for, status.entanglement_saturated_for);
    switch(status.algo_type) {
        case AlgorithmType::iDMRG:
            report += fmt::format("sat: [σ² {:<1} Sₑ {:<1}] ", status.variance_mpo_saturated_for, status.entanglement_saturated_for);
            break;
        case AlgorithmType::iTEBD: report += fmt::format("sat: [Sₑ {:<1}] ", status.entanglement_saturated_for); break;
        default: throw except::runtime_error("Wrong simulation type");
    }
    report += fmt::format("con: {:<4} ", status.algorithm_converged_for);
    report += fmt::format("time:{:>8.2f}s ", tid::get_unscoped("t_tot").get_time());
    report += fmt::format("mem: [rss {:<.1f} peak {:<.1f} vm {:<.1f}] MB ", debug::mem_rss_in_mb(), debug::mem_hwm_in_mb(), debug::mem_vm_in_mb());
    tools::log->info(report);
}

void AlgorithmInfinite::print_status_full() {
    //    compute_observables();
    tensors.state->do_all_measurements();
    using namespace std;
    using namespace tools::infinite::measure;
    tools::log->info("{:=^60}", "");
    tools::log->info("= {: ^56} =", fmt::format("Completed [{}][{}]", status.algo_type_sv(), tensors.state->get_name()));
    tools::log->info("Iterations            = {:<16d}", status.iter);
    switch(status.algo_type) {
        case AlgorithmType::iDMRG:
            tools::log->info("Energy MPO            = {:<16.16f}", tools::infinite::measure::energy_per_site_mpo(tensors));
            tools::log->info("Energy HAM            = {:<16.16f}", tools::infinite::measure::energy_per_site_ham(tensors));
            tools::log->info("Energy MOM            = {:<16.16f}", tools::infinite::measure::energy_per_site_mom(tensors));
            break;
        case AlgorithmType::iTEBD:
            tools::log->info("Energy HAM            = {:<16.16f}", tools::infinite::measure::energy_per_site_ham(tensors));
            tools::log->info("Energy MOM            = {:<16.16f}", tools::infinite::measure::energy_per_site_mom(tensors));
            break;
        default: throw except::runtime_error("Wrong simulation type");
    }
    switch(status.algo_type) {
        case AlgorithmType::iDMRG:
            tools::log->info("lg σ²H MPO         = {:<8.2e}", tools::infinite::measure::energy_variance_per_site_mpo(tensors));
            tools::log->info("lg σ²H HAM         = {:<8.2e}", tools::infinite::measure::energy_variance_per_site_ham(tensors));
            tools::log->info("lg σ²H MOM         = {:<8.2e}", tools::infinite::measure::energy_variance_per_site_mom(tensors));
            break;
        case AlgorithmType::iTEBD:
            tools::log->info("lg σ²H HAM         = {:<8.2e}", tools::infinite::measure::energy_variance_per_site_ham(tensors));
            tools::log->info("lg σ²H MOM         = {:<8.2e}", tools::infinite::measure::energy_variance_per_site_mom(tensors));
            break;
        default: throw except::runtime_error("Wrong simulation type");
    }
    tools::log->info("Truncation error      = {:<8.2e}", tools::infinite::measure::truncation_error(*tensors.state));
    tools::log->info("Entanglement Entropy  = {:<16.16f}", tools::infinite::measure::entanglement_entropy(*tensors.state));
    tools::log->info("χmax                  = {:<16d}", settings::get_bond_max(status.algo_type));
    tools::log->info("χ                     = {:<16d}", tools::infinite::measure::bond_dimension(*tensors.state));

    switch(status.algo_type) {
        case AlgorithmType::iDMRG: break;
        case AlgorithmType::iTEBD: tools::log->info("δt                    = {:<16.16f}", status.delta_t); break;

        default: throw except::runtime_error("Wrong simulation type");
    }

    tools::log->info("Algorithm succeeded      = {:<}", status.algorithm_has_succeeded);
    tools::log->info("Algorithm saturated for  = {:<}", status.algorithm_saturated_for);
    tools::log->info("Algorithm converged for  = {:<}", status.algorithm_converged_for);
    tools::log->info("Algorithm has stuck for  = {:<}", status.algorithm_has_stuck_for);
    tools::log->info("σ² MPO                   = Converged : {:<4}  Saturated: {:<4}", status.variance_mpo_converged_for, status.variance_mpo_saturated_for);
    tools::log->info("σ² HAM                   = Converged : {:<4}  Saturated: {:<4}", status.variance_ham_converged_for, status.variance_ham_saturated_for);
    tools::log->info("σ² MOM                   = Converged : {:<4}  Saturated: {:<4}", status.variance_mom_converged_for, status.variance_mom_saturated_for);
    tools::log->info("Sₑ                       = Converged : {:<4}  Saturated: {:<4}", status.entanglement_converged_for, status.entanglement_saturated_for);
    tools::log->info("Time                     = {:<16.16f}", tid::get_unscoped("t_tot").get_time());
    tools::log->info("Mem RSS                  = {:<.1f} MB", debug::mem_rss_in_mb());
    tools::log->info("Mem Peak                 = {:<.1f} MB", debug::mem_hwm_in_mb());
    tools::log->info("Mem VM                   = {:<.1f} MB", debug::mem_vm_in_mb());
    tools::log->info("{:=^60}", "");
}
