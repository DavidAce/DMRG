//
// Created by david on 2019-11-07.
//
#include <config/nmspc_settings.h>
#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/io.h>
#include <tools/finite/measure.h>

namespace tools::finite::io::h5table {
    void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile, const std::string &link) {
        if(save_log.empty()) {
            try {
                if(h5ppFile.linkExists(link)) {
                    auto step      = h5ppFile.readAttribute<uint64_t>("step", link);
                    auto iter      = h5ppFile.readAttribute<uint64_t>("iteration", link);
                    save_log[link] = std::make_pair(iter, step);
                }
            } catch(const std::exception &ex) { tools::log->warn("Could not bootstrap save_log: {}", ex.what()); }
        }
    }
}

/*! Write down the Hamiltonian model type and site info as attributes */
void tools::finite::io::h5table::save_model(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                            const class_model_finite &model) {
    if(storage_level == StorageLevel::NONE) return;
    if(h5ppFile.linkExists(table_path)) return tools::log->debug("The hamiltonian has already been written to [{}]", table_path);

    tools::log->trace("Storing table: [{}]", table_path);
    auto t_hdf = tools::common::profile::get_default_prof()["t_hdf"]->tic_token();
    for(auto site = 0ul; site < model.get_length(); site++) model.get_mpo(site).save_hamiltonian(h5ppFile, table_path);
    h5ppFile.writeAttribute(enum2str(settings::model::model_type), "model_type", table_path);
    h5ppFile.writeAttribute(settings::model::model_size, "model_size", table_path);
}

void tools::finite::io::h5table::save_measurements(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                   const class_tensors_finite &tensors, const class_algorithm_status &status, AlgorithmType algo_type){
    save_measurements(h5ppFile, table_path, storage_level, *tensors.state, tensors, status, algo_type);
}

void tools::finite::io::h5table::save_measurements(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                   const class_state_finite & state, const class_tensors_finite &tensors, const class_algorithm_status &status, AlgorithmType algo_type) {
    if(storage_level == StorageLevel::NONE) return;
    // Check if the current entry has already been appended
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    bootstrap_save_log(save_log, h5ppFile, table_path);
    auto save_point = std::make_pair(status.iter, status.step);
    if(save_log[table_path] == save_point) return;

    log->trace("Appending to table: {}", table_path);
    h5pp_table_measurements_finite::register_table_type();
    if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_measurements_finite::h5_type, table_path, "measurements");

    h5pp_table_measurements_finite::table measurement_entry{};
    measurement_entry.step                          = static_cast<uint64_t>(status.step);
    measurement_entry.iter                          = static_cast<uint64_t>(status.iter);
    measurement_entry.position                      = static_cast<long>(status.position);
    measurement_entry.length                        = static_cast<uint64_t>(tools::finite::measure::length(state));
    measurement_entry.bond_dimension_midchain       = static_cast<long>(tools::finite::measure::bond_dimension_midchain(state));
    measurement_entry.bond_dimension_current        = static_cast<long>(tools::finite::measure::bond_dimension_current(state));
    measurement_entry.bond_dimension_limit          = status.chi_lim;
    measurement_entry.bond_dimension_maximum        = status.chi_lim_max;
    measurement_entry.entanglement_entropy_midchain = tools::finite::measure::entanglement_entropy_midchain(state);
    measurement_entry.entanglement_entropy_current  = tools::finite::measure::entanglement_entropy_current(state);
    measurement_entry.number_entropy_midchain       = tools::finite::measure::number_entropy_midchain(state);
    measurement_entry.number_entropy_current        = tools::finite::measure::number_entropy_current(state);
    measurement_entry.norm                          = tools::finite::measure::norm(state);
    measurement_entry.energy                        = tools::finite::measure::energy(state,tensors);
    measurement_entry.energy_per_site               = tools::finite::measure::energy_per_site(state,tensors);
    if(algo_type != AlgorithmType::fLBIT) {
        measurement_entry.energy_variance                 = tools::finite::measure::energy_variance(state,tensors);
        measurement_entry.energy_variance_per_site        = tools::finite::measure::energy_variance_per_site(state,tensors);
        measurement_entry.energy_variance_lowest          = status.energy_variance_lowest;
        measurement_entry.energy_variance_per_site_lowest = status.energy_variance_lowest / static_cast<double>(tensors.get_length());
    }
    measurement_entry.spin_components  = tools::finite::measure::spin_components(state);
    measurement_entry.truncation_error = state.get_truncation_error_midchain();
    measurement_entry.total_time       = status.wall_time;
    measurement_entry.algorithm_time   = status.algo_time;
    measurement_entry.physical_time    = status.phys_time;
    auto t_hdf = tools::common::profile::get_default_prof()["t_hdf"]->tic_token();
    h5ppFile.appendTableRecords(measurement_entry, table_path);
    h5ppFile.writeAttribute(status.iter, "iteration", table_path);
    h5ppFile.writeAttribute(status.step, "step", table_path);
    save_log[table_path] = save_point;
}

void tools::finite::io::h5table::save_sim_status(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level,
                                                 const class_algorithm_status &status) {
    tools::common::io::h5table::save_sim_status(h5ppFile, table_prefix, storage_level, status);
}

void tools::finite::io::h5table::save_profiling(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level,
                                                const class_algorithm_status &status) {
    tools::common::io::h5table::save_profiling(h5ppFile, table_prefix, storage_level, status);
}

void tools::finite::io::h5table::save_mem_usage(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level,
                                                const class_algorithm_status &status) {
    tools::common::io::h5table::save_mem_usage(h5ppFile, table_prefix, storage_level, status);
}