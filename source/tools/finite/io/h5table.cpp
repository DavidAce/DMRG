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

/*! Write down the Hamiltonian model type and site info as attributes */
void tools::finite::io::h5table::write_model(h5pp::File &h5ppFile, const std::string &model_prefix, const StorageLevel &storage_level,
                                             const class_model_finite &model) {
    if(storage_level == StorageLevel::NONE) return;
    std::string table_path = model_prefix + "/Hamiltonian";
    if(h5ppFile.linkExists(table_path)) return tools::log->debug("The hamiltonian has already been written to [{}]", table_path);

    tools::log->trace("Storing table: [{}]", table_path);
    tools::common::profile::t_hdf->tic();
    for(auto site = 0ul; site < model.get_length(); site++) model.get_mpo(site).write_hamiltonian(h5ppFile, table_path);
    h5ppFile.writeAttribute(enum2str(settings::model::model_type), "model_type", table_path);
    h5ppFile.writeAttribute(settings::model::model_size, "model_size", table_path);

    tools::common::profile::t_hdf->toc();
}

void tools::finite::io::h5table::write_measurements(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level,
                                                    const class_tensors_finite &tensors, const class_algorithm_status &status) {
    if(storage_level == StorageLevel::NONE) return;
    std::string table_path = table_prefix + "/measurements";
    log->trace("Appending to table: {}", table_path);
    h5pp_table_measurements_finite::register_table_type();
    if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_measurements_finite::h5_type, table_path, "measurements");

    h5pp_table_measurements_finite::table measurement_entry{};
    measurement_entry.step                            = static_cast<uint64_t>(status.step);
    measurement_entry.iter                            = static_cast<uint64_t>(status.iter);
    measurement_entry.position                        = static_cast<uint64_t>(status.position);
    measurement_entry.length                          = static_cast<uint64_t>(tools::finite::measure::length(tensors));
    measurement_entry.bond_dimension_midchain         = static_cast<long>(tools::finite::measure::bond_dimension_midchain(*tensors.state));
    measurement_entry.bond_dimension_current          = static_cast<long>(tools::finite::measure::bond_dimension_current(*tensors.state));
    measurement_entry.bond_dimension_limit            = tensors.state->get_chi_lim();
    measurement_entry.bond_dimension_maximum          = tensors.state->get_chi_max();
    measurement_entry.entanglement_entropy_midchain   = tools::finite::measure::entanglement_entropy_midchain(*tensors.state);
    measurement_entry.entanglement_entropy_current    = tools::finite::measure::entanglement_entropy_current(*tensors.state);
    measurement_entry.norm                            = tools::finite::measure::norm(*tensors.state);
    measurement_entry.energy                          = tools::finite::measure::energy(tensors);
    measurement_entry.energy_per_site                 = tools::finite::measure::energy_per_site(tensors);
    measurement_entry.energy_variance                 = tools::finite::measure::energy_variance(tensors);
    measurement_entry.energy_variance_per_site        = tools::finite::measure::energy_variance_per_site(tensors);
    measurement_entry.energy_variance_lowest          = tensors.state->lowest_recorded_variance;
    measurement_entry.energy_variance_per_site_lowest = tensors.state->lowest_recorded_variance / static_cast<double>(tensors.get_length());
    measurement_entry.spin_components                 = tools::finite::measure::spin_components(*tensors.state);
    measurement_entry.truncation_error                = tensors.state->get_truncation_error();
    measurement_entry.wall_time                       = status.wall_time;
    tools::common::profile::t_hdf->tic();
    h5ppFile.appendTableEntries(measurement_entry, table_path);
    tools::common::profile::t_hdf->toc();
}

void tools::finite::io::h5table::write_sim_status(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level,
                                                  const class_algorithm_status &status) {
    tools::common::io::h5table::write_sim_status(h5ppFile, table_prefix, storage_level, status);
}

void tools::finite::io::h5table::write_profiling(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level,
                                                 const class_algorithm_status &status) {
    tools::common::io::h5table::write_profiling(h5ppFile, table_prefix, storage_level, status);
}

void tools::finite::io::h5table::write_mem_usage(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level,
                                                 const class_algorithm_status &status) {
    tools::common::io::h5table::write_mem_usage(h5ppFile, table_prefix, storage_level, status);
}