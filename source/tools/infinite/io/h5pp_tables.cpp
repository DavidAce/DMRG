//
// Created by david on 2019-11-07.
//

#include <tools/infinite/io.h>
#include <tools/infinite/measure.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <state/class_state_infinite.h>
#include <io/table_types.h>
#include <h5pp/h5pp.h>

void tools::infinite::io::h5table::write_measurements(const class_state_infinite &state, const class_simulation_status &sim_status, h5pp::File & h5ppFile, const std::string & table_path) {
    log->trace("Appending measurement entry to table: {}...",table_path);
    h5pp_table_measurements_infinite::register_table_type();
    if(not h5ppFile.linkExists(table_path))
        h5ppFile.createTable(h5pp_table_measurements_infinite::h5_type,table_path, "measurements");

    h5pp_table_measurements_infinite::table measurement_entry;
    measurement_entry.step                            = sim_status.step;
    measurement_entry.iteration                       = sim_status.iteration;
    measurement_entry.position                        = sim_status.position;
    measurement_entry.length                          = tools::infinite::measure::length(state);
    measurement_entry.bond_dimension                  = tools::infinite::measure::bond_dimension(state);
    measurement_entry.bond_dimension_limit            = state.get_chi_lim();
    measurement_entry.bond_dimension_maximum          = state.get_chi_max();
    measurement_entry.entanglement_entropy            = tools::infinite::measure::entanglement_entropy(state);
    measurement_entry.norm                            = tools::infinite::measure::norm(state);
    measurement_entry.energy_mpo                      = tools::infinite::measure::energy_mpo(state);
    measurement_entry.energy_per_site_mpo             = tools::infinite::measure::energy_per_site_mpo(state);
    measurement_entry.energy_per_site_ham             = tools::infinite::measure::energy_per_site_ham(state);
    measurement_entry.energy_per_site_mom             = tools::infinite::measure::energy_per_site_mom(state);
    measurement_entry.energy_variance_mpo             = tools::infinite::measure::energy_variance_mpo(state);
    measurement_entry.energy_variance_per_site_mpo    = tools::infinite::measure::energy_variance_per_site_mpo(state);
    measurement_entry.energy_variance_per_site_ham    = tools::infinite::measure::energy_variance_per_site_ham(state);
    measurement_entry.energy_variance_per_site_mom    = tools::infinite::measure::energy_variance_per_site_mom(state);
    measurement_entry.truncation_error                = tools::infinite::measure::truncation_error(state);
    measurement_entry.wall_time                       = sim_status.wall_time;
    measurement_entry.phys_time                       = sim_status.phys_time;
    measurement_entry.time_step                       = sim_status.delta_t;
    h5ppFile.appendTableEntries(measurement_entry,table_path);
    log->trace("Appending measurement entry to table: {}... OK",table_path);
}




void tools::infinite::io::h5table::write_sim_status(const class_simulation_status &sim_status, h5pp::File & h5ppFile, const std::string & table_path) {
    tools::common::io::h5table::write_sim_status(sim_status,h5ppFile,table_path);
}

void tools::infinite::io::h5table::write_profiling(const class_simulation_status &sim_status, h5pp::File & h5ppFile, const std::string & table_path) {
    tools::common::io::h5table::write_profiling(sim_status,h5ppFile,table_path);
}
