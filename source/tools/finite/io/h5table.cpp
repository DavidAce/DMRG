//
// Created by david on 2019-11-07.
//
#include <tools/finite/io.h>
#include <tools/finite/measure.h>
#include <tools/common/log.h>
#include <tools/common/io.h>

#include <io/table_types.h>
#include <state/class_state_finite.h>
#include <general/nmspc_quantum_mechanics.h>

#include <h5pp/h5pp.h>

void tools::finite::io::h5table::write_measurements(h5pp::File & h5ppFile, const std::string &path, const StorageLevel & storage_level, const class_simulation_status & sim_status, const class_state_finite & state) {
    log->trace("Appending measurement entry to table: {}...", path);
    h5pp_table_measurements_finite::register_table_type();
    if(not h5ppFile.linkExists(path))
        h5ppFile.createTable(h5pp_table_measurements_finite::h5_type, path, "measurements");

    h5pp_table_measurements_finite::table measurement_entry;
    measurement_entry.step                            = sim_status.step;
    measurement_entry.iteration                       = sim_status.iteration;
    measurement_entry.position                        = sim_status.position;
    measurement_entry.length                          = tools::finite::measure::length(state);
    measurement_entry.bond_dimension_midchain         = tools::finite::measure::bond_dimension_midchain(state);
    measurement_entry.bond_dimension_current          = tools::finite::measure::bond_dimension_current(state);
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
    measurement_entry.energy_variance_per_site_lowest = state.lowest_recorded_variance/state.get_length();
    measurement_entry.spin_component_sx               = tools::finite::measure::spin_component(state,qm::spinOneHalf::sx);
    measurement_entry.spin_component_sy               = tools::finite::measure::spin_component(state,qm::spinOneHalf::sy);
    measurement_entry.spin_component_sz               = tools::finite::measure::spin_component(state,qm::spinOneHalf::sz);
    measurement_entry.truncation_error                = state.get_truncation_error();
    measurement_entry.wall_time                       = sim_status.wall_time;

    h5ppFile.appendTableEntries(measurement_entry, path);
    log->trace("Appending measurement entry to table: {}... OK", path);
}



void tools::finite::io::h5table::write_sim_status(h5pp::File & h5ppFile, const std::string & table_path, const StorageLevel & storage_level, const class_simulation_status & sim_status) {
    tools::common::io::h5table::write_sim_status(h5ppFile, table_path, storage_level, sim_status);
}

void tools::finite::io::h5table::write_profiling(h5pp::File & h5ppFile, const std::string & table_path, const StorageLevel & storage_level, const class_simulation_status & sim_status) {
    tools::common::io::h5table::write_profiling(h5ppFile, table_path, storage_level, sim_status);
}
