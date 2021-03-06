//
// Created by david on 2019-11-07.
//

#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <tensors/class_tensors_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/infinite/io.h>
#include <tools/infinite/measure.h>

/*! Write down the Hamiltonian model type and site info as attributes */
void tools::infinite::io::h5table::save_model(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                              const class_model_infinite &model) {
    if(storage_level < StorageLevel::LIGHT) return;
    tools::log->trace("Writing Hamiltonian model");
    if(h5ppFile.linkExists(table_path)) return tools::log->debug("The hamiltonian has already been written to [{}]", table_path);

    tools::log->trace("Storing table: [{}]", table_path);
    auto t_hdf = tools::common::profile::get_default_prof()["t_hdf"]->tic_token();
    model.get_mpo_siteA().save_hamiltonian(h5ppFile, table_path);
    model.get_mpo_siteB().save_hamiltonian(h5ppFile, table_path);
    h5ppFile.writeAttribute(enum2str(settings::model::model_type), "model_type", table_path);
    h5ppFile.writeAttribute(settings::model::model_size, "model_size", table_path);
}

void tools::infinite::io::h5table::save_measurements(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                     const class_tensors_infinite &tensors, const class_algorithm_status &status) {
    if(storage_level == StorageLevel::NONE) return;

    // Check if the current entry has already been appended
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    auto                                                                  save_point = std::make_pair(status.iter, status.step);
    if(save_log[table_path] == save_point) return;

    log->trace("Appending to table: {}", table_path);
    h5pp_table_measurements_infinite::register_table_type();
    if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_measurements_infinite::h5_type, table_path, "measurements");

    h5pp_table_measurements_infinite::table measurement_entry{};
    const auto &                            state = *tensors.state;

    measurement_entry.step                   = static_cast<uint64_t>(status.step);
    measurement_entry.iter                   = static_cast<uint64_t>(status.iter);
    measurement_entry.position               = static_cast<uint64_t>(status.position);
    measurement_entry.length                 = static_cast<uint64_t>(tools::infinite::measure::length(tensors));
    measurement_entry.bond_dimension         = tools::infinite::measure::bond_dimension(state);
    measurement_entry.bond_dimension_limit   = status.chi_lim;
    measurement_entry.bond_dimension_maximum = status.chi_lim_max;
    measurement_entry.entanglement_entropy   = tools::infinite::measure::entanglement_entropy(state);
    measurement_entry.norm                   = tools::infinite::measure::norm(*tensors.state);
    if(std::abs(status.delta_t) == 0) {
        // MPO calculations do not make sense for iTEBD, and we know delta_t != 0 on iTEBD.
        measurement_entry.energy_mpo                   = tools::infinite::measure::energy_mpo(tensors);
        measurement_entry.energy_per_site_mpo          = tools::infinite::measure::energy_per_site_mpo(tensors);
        measurement_entry.energy_variance_mpo          = tools::infinite::measure::energy_variance_mpo(tensors);
        measurement_entry.energy_variance_per_site_mpo = tools::infinite::measure::energy_variance_per_site_mpo(tensors);
    }
    measurement_entry.energy_per_site_ham          = tools::infinite::measure::energy_per_site_ham(tensors);
    measurement_entry.energy_per_site_mom          = tools::infinite::measure::energy_per_site_mom(tensors);
    measurement_entry.energy_variance_per_site_ham = tools::infinite::measure::energy_variance_per_site_ham(tensors);
    measurement_entry.energy_variance_per_site_mom = tools::infinite::measure::energy_variance_per_site_mom(tensors);
    measurement_entry.truncation_error             = tools::infinite::measure::truncation_error(state);
    measurement_entry.wall_time                    = status.wall_time;
    measurement_entry.phys_time                    = status.phys_time;
    measurement_entry.time_step                    = status.delta_t;
    auto t_hdf = tools::common::profile::get_default_prof()["t_hdf"]->tic_token();
    h5ppFile.appendTableRecords(measurement_entry, table_path);
    save_log[table_path] = save_point;
}

void tools::infinite::io::h5table::save_sim_status(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                   const class_algorithm_status &status) {
    tools::common::io::h5table::save_sim_status(h5ppFile, table_path, storage_level, status);
}

void tools::infinite::io::h5table::save_profiling(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                  const class_algorithm_status &status) {
    tools::common::io::h5table::save_profiling(h5ppFile, table_path, storage_level, status);
}

void tools::infinite::io::h5table::save_mem_usage(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                  const class_algorithm_status &status) {
    tools::common::io::h5table::save_mem_usage(h5ppFile, table_path, storage_level, status);
}