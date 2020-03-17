//
// Created by david on 2019-11-07.
//
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <io/table_types.h>
#include <h5pp/h5pp.h>
#include <simulation/enums.h>

void tools::common::io::h5table::write_sim_status(h5pp::File & h5ppFile, const std::string &path, const StorageLevel &storage_level, const class_simulation_status &sim_status) {
    if(storage_level == StorageLevel::NONE) return;
    log->trace("Appending simulation status table: {}...", path);
    h5pp_table_sim_status::register_table_type();
    if(not h5ppFile.linkExists(path))
        h5ppFile.createTable(h5pp_table_sim_status::h5_type, path, "simulation status");
    h5ppFile.appendTableEntries(sim_status, path);
    log->trace("Appending simulation status table: {}... OK", path);
}


void tools::common::io::h5table::write_profiling(h5pp::File & h5ppFile, const std::string &path, const StorageLevel &storage_level, const class_simulation_status &sim_status) {
    log->trace("Appending profiling data to table: {}...", path);
    h5pp_table_profiling::register_table_type();
    if(not h5ppFile.linkExists(path))
        h5ppFile.createTable(h5pp_table_profiling::h5_type, path, "profiling");

    h5pp_table_profiling::table profiling_entry;
    profiling_entry.step            = sim_status.step;
    profiling_entry.iteration       = sim_status.iteration;
    profiling_entry.position        = sim_status.position;
    profiling_entry.t_tot           = tools::common::profile::t_tot->get_age();
    profiling_entry.t_pre           = tools::common::profile::t_pre->get_measured_time();
    profiling_entry.t_pos           = tools::common::profile::t_pos->get_measured_time();
    profiling_entry.t_sim           = tools::common::profile::t_sim->get_measured_time();
    profiling_entry.t_con           = tools::common::profile::t_con->get_measured_time();
    profiling_entry.t_eig           = tools::common::profile::t_eig->get_measured_time();
    profiling_entry.t_svd           = tools::common::profile::t_svd->get_measured_time();
    profiling_entry.t_ham           = tools::common::profile::t_ham->get_measured_time();
    profiling_entry.t_ham_sq        = tools::common::profile::t_ham_sq->get_measured_time();
    profiling_entry.t_mpo           = tools::common::profile::t_mpo->get_measured_time();
    profiling_entry.t_opt           = tools::common::profile::t_opt->get_measured_time();
    profiling_entry.t_evo           = tools::common::profile::t_evo->get_measured_time();
    profiling_entry.t_env           = tools::common::profile::t_env->get_measured_time();
    profiling_entry.t_ent           = tools::common::profile::t_ent->get_measured_time();
    profiling_entry.t_ene           = tools::common::profile::t_ene->get_measured_time();
    profiling_entry.t_var           = tools::common::profile::t_var->get_measured_time();
    profiling_entry.t_prj           = tools::common::profile::t_prj->get_measured_time();
    profiling_entry.t_chk           = tools::common::profile::t_chk->get_measured_time();
    profiling_entry.t_hdf           = tools::common::profile::t_hdf->get_measured_time();
    profiling_entry.t_ene_ham       = tools::common::profile::t_ene_ham->get_measured_time();
    profiling_entry.t_ene_mom       = tools::common::profile::t_ene_mom->get_measured_time();
    profiling_entry.t_var_ham       = tools::common::profile::t_var_ham->get_measured_time();
    profiling_entry.t_var_mom       = tools::common::profile::t_var_mom->get_measured_time();
    profiling_entry.t_vH2v          = tools::common::profile::t_vH2v->get_measured_time();
    profiling_entry.t_vHv           = tools::common::profile::t_vHv->get_measured_time();
    profiling_entry.t_vH2           = tools::common::profile::t_vH2->get_measured_time();
    profiling_entry.t_vH            = tools::common::profile::t_vH->get_measured_time();
    profiling_entry.t_op            = tools::common::profile::t_op->get_measured_time();
    h5ppFile.appendTableEntries(profiling_entry, path);
    tools::log->trace("Appending profiling data to table: {}... OK", path);

}
