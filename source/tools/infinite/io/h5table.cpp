//
// Created by david on 2019-11-07.
//

#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <model/class_mpo_base.h>
#include <state/class_state_infinite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/infinite/io.h>
#include <tools/infinite/measure.h>

/*! Write down the Hamiltonian model type and site info as attributes */
void tools::infinite::io::h5table::write_model(h5pp::File &h5ppFile, const std::string &model_prefix, const StorageLevel &storage_level,const class_state_infinite &state){
    if(storage_level < StorageLevel::LIGHT) return;
    tools::log->trace("Writing Hamiltonian model");
    std::string table_path = model_prefix + "/Hamiltonian";
    if(h5ppFile.linkExists(table_path)) return tools::log->debug("The hamiltonian has already been written to [{}]", table_path);

    tools::log->trace("Storing table: [{}]", table_path);
    tools::common::profile::t_hdf->tic();
    state.HA->write_hamiltonian(h5ppFile, table_path);
    state.HB->write_hamiltonian(h5ppFile, table_path);
    h5ppFile.writeAttribute(enum2str(settings::model::model_type), "model_type", table_path);
    h5ppFile.writeAttribute(settings::model::model_size, "model_size", table_path);
    tools::common::profile::t_hdf->toc();
}


void tools::infinite::io::h5table::write_measurements(h5pp::File & h5ppFile, const std::string &table_prefix, const StorageLevel & storage_level, const class_state_infinite & state, const class_simulation_status & sim_status) {
    if(storage_level == StorageLevel::NONE) return;
    std::string table_path = table_prefix + "/measurements";
    log->trace("Appending to table: {}", table_path);
    h5pp_table_measurements_finite::register_table_type();
    if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_measurements_infinite::h5_type, table_path, "measurements");

    h5pp_table_measurements_infinite::table measurement_entry{};
    measurement_entry.step                            = (uint64_t) sim_status.step;
    measurement_entry.iter                            = (uint64_t) sim_status.iter;
    measurement_entry.position                        = (uint64_t) sim_status.position;
    measurement_entry.length                          = (uint64_t) tools::infinite::measure::length(state);
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
    tools::common::profile::t_hdf->tic();
    h5ppFile.appendTableEntries(measurement_entry, table_path);
    tools::common::profile::t_hdf->toc();

}




void tools::infinite::io::h5table::write_sim_status(h5pp::File & h5ppFile, const std::string &table_prefix, const StorageLevel & storage_level, const class_simulation_status &sim_status) {
    tools::common::io::h5table::write_sim_status(h5ppFile, table_prefix,storage_level, sim_status);
}

void tools::infinite::io::h5table::write_profiling(h5pp::File & h5ppFile, const std::string &table_prefix, const StorageLevel & storage_level, const class_simulation_status &sim_status) {
    tools::common::io::h5table::write_profiling(h5ppFile, table_prefix,storage_level,sim_status);
}

void tools::infinite::io::h5table::write_mem_usage(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level,
                                                 const class_simulation_status &sim_status) {
    tools::common::io::h5table::write_mem_usage(h5ppFile, table_prefix, storage_level, sim_status);
}