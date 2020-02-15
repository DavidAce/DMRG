//
// Created by david on 2019-11-07.
//
#include <tools/finite/io.h>
#include <tools/finite/measure.h>
#include <tools/common/log.h>
#include <tools/common/io.h>

#include <io/table_types.h>
#include <io/class_h5table_buffer.h>
#include <state/class_state_finite.h>
#include <general/nmspc_quantum_mechanics.h>
void tools::finite::io::h5table::write_measurements(const class_state_finite &state, const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_measurements_finite> &h5tbuf) {
    log->trace("Appending measurement entry to table: {}...",h5tbuf.get_table_name());
    class_h5table_measurements_finite::data_struct measurements_entry;
    measurements_entry.step                            = sim_status.step;
    measurements_entry.iteration                       = sim_status.iteration;
    measurements_entry.position                        = sim_status.position;
    measurements_entry.length                          = tools::finite::measure::length(state);
    measurements_entry.bond_dimension_midchain         = tools::finite::measure::bond_dimension_midchain(state);
    measurements_entry.bond_dimension_current          = tools::finite::measure::bond_dimension_current(state);
    measurements_entry.bond_dimension_limit            = state.get_chi_lim();
    measurements_entry.bond_dimension_maximum          = state.get_chi_max();
    measurements_entry.entanglement_entropy_midchain   = tools::finite::measure::entanglement_entropy_midchain(state);
    measurements_entry.entanglement_entropy_current    = tools::finite::measure::entanglement_entropy_current(state);
    measurements_entry.norm                            = tools::finite::measure::norm(state);
    measurements_entry.energy                          = tools::finite::measure::energy(state);
    measurements_entry.energy_per_site                 = tools::finite::measure::energy_per_site(state);
    measurements_entry.energy_variance                 = tools::finite::measure::energy_variance(state);
    measurements_entry.energy_variance_per_site        = tools::finite::measure::energy_variance_per_site(state);
    measurements_entry.energy_variance_lowest          = state.lowest_recorded_variance;
    measurements_entry.energy_variance_per_site_lowest = state.lowest_recorded_variance/state.get_length();
    measurements_entry.spin_component_sx               = tools::finite::measure::spin_component(state,qm::spinOneHalf::sx);
    measurements_entry.spin_component_sy               = tools::finite::measure::spin_component(state,qm::spinOneHalf::sy);
    measurements_entry.spin_component_sz               = tools::finite::measure::spin_component(state,qm::spinOneHalf::sz);
    measurements_entry.truncation_error                = state.get_truncation_error();
    measurements_entry.wall_time                       = sim_status.wall_time;
    h5tbuf.append_record(measurements_entry);
    log->trace("Appending measurement entry to table: {}... OK",h5tbuf.get_table_name());
}

void tools::finite::io::h5table::write_sim_status(const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_simulation_status> &h5tbuf) {
    tools::common::io::h5table::write_sim_status(sim_status,h5tbuf);
}

void tools::finite::io::h5table::write_profiling(const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_profiling> &h5tbuf) {
    tools::common::io::h5table::write_profiling(sim_status,h5tbuf);
}
