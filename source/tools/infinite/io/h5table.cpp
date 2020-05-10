//
// Created by david on 2019-11-07.
//

#include <tools/infinite/io.h>
#include <tools/infinite/measure.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <model/class_model_base.h>
#include <state/class_state_infinite.h>
#include <io/table_types.h>
#include <h5pp/h5pp.h>

/*! Write down the Hamiltonian model type and site info as attributes */
void tools::infinite::io::h5table::write_model(h5pp::File & h5ppFile, const std::string &prefix, const StorageLevel & storage_level, const class_state_infinite & state) {
    if(storage_level < StorageLevel::LIGHT) return;
    tools::log->trace("Writing Hamiltonian model");
    tools::common::profile::t_hdf->tic();
    std::string table_name = prefix + "/model/Hamiltonian";
    state.HA->write_hamiltonian(h5ppFile, table_name);
    state.HB->write_hamiltonian(h5ppFile, table_name);
    h5ppFile.writeAttribute(settings::model::model_type, "model_type", prefix + "/model/Hamiltonian");
    h5ppFile.writeAttribute(settings::model::model_size, "model_size", prefix + "/model/Hamiltonian");
    tools::common::profile::t_hdf->toc();
}


void tools::infinite::io::h5table::write_measurements(h5pp::File & h5ppFile, const std::string &prefix, const StorageLevel & storage_level, const class_simulation_status & sim_status, const class_state_infinite & state) {
    if(storage_level == StorageLevel::NONE) return;
    log->trace("Appending measurement entry to table: {}", prefix);
    h5pp_table_measurements_infinite::register_table_type();
    if(not h5ppFile.linkExists(prefix))
        h5ppFile.createTable(h5pp_table_measurements_infinite::h5_type, prefix, "measurements");

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
    h5ppFile.appendTableEntries(measurement_entry, prefix);
}




void tools::infinite::io::h5table::write_sim_status(h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_simulation_status &sim_status) {
    tools::common::io::h5table::write_sim_status(h5ppFile, path,storage_level, sim_status);
}

void tools::infinite::io::h5table::write_profiling(h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_simulation_status &sim_status) {
    tools::common::io::h5table::write_profiling(h5ppFile,path,storage_level,sim_status);
}
