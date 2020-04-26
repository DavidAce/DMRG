//
// Created by david on 2019-11-07.
//
#include <general/nmspc_quantum_mechanics.h>
#include <io/table_types.h>
#include <state/class_state_finite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/io.h>
#include <tools/finite/measure.h>

#include <h5pp/h5pp.h>
#include <simulation/nmspc_settings.h>

/*! Write down the Hamiltonian model type and site info as attributes */
void tools::finite::io::h5table::write_model(h5pp::File &h5ppFile, const std::string &prefix, const StorageLevel &storage_level,
                                             const class_state_finite &state) {
    if(storage_level < StorageLevel::LIGHT) return;
    std::string table_path = prefix.substr(0,prefix.find('/',1)) +  + "/model/Hamiltonian";
    if(h5ppFile.linkExists(table_path))
        return tools::log->debug("The hamiltonian has already been written to [{}]", table_path);

    tools::log->trace("Writing Hamiltonian model");
    tools::common::profile::t_hdf->tic();
    for(auto site = 0ul; site < state.get_length(); site++) state.get_MPO(site).write_parameters(h5ppFile, table_path);
    h5ppFile.writeAttribute(settings::model::model_type, "model_type", table_path);
    h5ppFile.writeAttribute(state.get_length(), "sites", table_path);
    h5ppFile.writeAttribute(state.get_position(), "position", table_path);
    tools::common::profile::t_hdf->toc();
}

void tools::finite::io::h5table::write_measurements(h5pp::File &h5ppFile, const std::string &prefix, const StorageLevel &storage_level,
                                                    const class_simulation_status &sim_status, const class_state_finite &state) {
    if(storage_level < StorageLevel::LIGHT) return;
    std::string table_path = prefix + "/measurements";
    log->trace("Appending measurement entry to table: {}", table_path);
    h5pp_table_measurements_finite::register_table_type();
    if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_measurements_finite::h5_type, table_path, "measurements");

    h5pp_table_measurements_finite::table measurement_entry;
    measurement_entry.step                            = (uint64_t) sim_status.step;
    measurement_entry.iter                            = (uint64_t) sim_status.iter;
    measurement_entry.position                        = (uint64_t) sim_status.position;
    measurement_entry.length                          = tools::finite::measure::length(state);
    measurement_entry.bond_dimension_midchain         = (long) tools::finite::measure::bond_dimension_midchain(state);
    measurement_entry.bond_dimension_current          = (long) tools::finite::measure::bond_dimension_current(state);
    measurement_entry.bond_dimension_limit            = state.get_chi_lim();
    measurement_entry.bond_dimension_maximum          = state.get_chi_max();
    measurement_entry.entanglement_entropy_midchain   = tools::finite::measure::entanglement_entropy_midchain(state);
    measurement_entry.entanglement_entropy_current    = tools::finite::measure::entanglement_entropy_current(state);
    measurement_entry.norm                            = tools::finite::measure::norm(state);
    measurement_entry.energy                          = tools::finite::measure::energy(state);
    measurement_entry.energy_per_site                 = tools::finite::measure::energy_per_site(state);
    measurement_entry.energy_variance                 = tools::finite::measure::energy_variance(state);
    measurement_entry.energy_variance_per_site        = tools::finite::measure::energy_variance_per_site(state);
    measurement_entry.energy_variance_lowest          = state.lowest_recorded_variance;
    measurement_entry.energy_variance_per_site_lowest = state.lowest_recorded_variance / (double)state.get_length();
    measurement_entry.spin_components                 = tools::finite::measure::spin_components(state);
    measurement_entry.truncation_error                = state.get_truncation_error();
    measurement_entry.wall_time                       = sim_status.wall_time;
    tools::common::profile::t_hdf->tic();
    h5ppFile.appendTableEntries(measurement_entry, table_path);
    tools::common::profile::t_hdf->toc();
}

void tools::finite::io::h5table::write_sim_status(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                  const class_simulation_status &sim_status) {
    tools::common::io::h5table::write_sim_status(h5ppFile, table_path, storage_level, sim_status);
}

void tools::finite::io::h5table::write_profiling(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                 const class_simulation_status &sim_status) {
    tools::common::io::h5table::write_profiling(h5ppFile, table_path, storage_level, sim_status);
}

void tools::finite::io::h5table::write_mem_usage(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                 const class_simulation_status &sim_status) {
    tools::common::io::h5table::write_mem_usage(h5ppFile, table_path, storage_level, sim_status);
}