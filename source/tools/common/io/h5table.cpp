//
// Created by david on 2019-11-07.
//
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <io/table_types.h>
#include <h5pp/h5pp.h>
#include <simulation/enums.h>

void tools::common::io::h5table::write_sim_status(h5pp::File & h5ppFile, const std::string &prefix, const StorageLevel &storage_level, const class_simulation_status &sim_status) {
    if(storage_level < StorageLevel::LIGHT) return;
    std::string table_path = prefix + "/sim_status";
    log->trace("Appending to table: {}", table_path);
    h5pp_table_sim_status::register_table_type();
    if(not h5ppFile.linkExists(table_path))
        h5ppFile.createTable(h5pp_table_sim_status::h5_type, table_path, "simulation status");
    h5ppFile.appendTableEntries(sim_status, table_path);
}


void tools::common::io::h5table::write_profiling(h5pp::File & h5ppFile, const std::string &prefix, const StorageLevel &storage_level, const class_simulation_status &sim_status) {
    if(storage_level < StorageLevel::LIGHT) return;
    std::string table_path = prefix + "/profiling";
    log->trace("Appending to table: {}", table_path);

    h5pp_table_profiling::register_table_type();
    if(not h5ppFile.linkExists(table_path))
        h5ppFile.createTable(h5pp_table_profiling::h5_type, table_path, "profiling");

    h5pp_table_profiling::table profiling_entry;
    profiling_entry.iter            = sim_status.iter;
    profiling_entry.step            = sim_status.step;
    profiling_entry.position        = sim_status.position;
    profiling_entry.t_tot           = tools::common::profile::t_tot->get_age();
    profiling_entry.t_pre           = tools::common::profile::t_pre->get_measured_time();
    profiling_entry.t_pos           = tools::common::profile::t_pos->get_measured_time();
    profiling_entry.t_sim           = tools::common::profile::t_sim->get_measured_time();
    profiling_entry.t_con           = tools::common::profile::t_con->get_measured_time();
    profiling_entry.t_eig           = tools::common::profile::t_eig->get_measured_time();
    profiling_entry.t_svd           = tools::common::profile::t_svd->get_measured_time();
    profiling_entry.t_ham           = tools::common::profile::t_ham->get_measured_time();
    profiling_entry.t_hsq           = tools::common::profile::t_hsq->get_measured_time();
    profiling_entry.t_mps           = tools::common::profile::t_mps->get_measured_time();
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
    h5ppFile.appendTableEntries(profiling_entry, table_path);
}

void tools::common::io::h5table::write_mem_usage(h5pp::File & h5ppFile, const std::string &prefix, const StorageLevel &storage_level, const class_simulation_status &sim_status) {
    if(storage_level < StorageLevel::LIGHT) return;
    std::string table_path = prefix + "/mem_usage";
    log->trace("Appending to table: {}", table_path);

    h5pp_table_memory_usage::register_table_type();
    if(not h5ppFile.linkExists(table_path))
        h5ppFile.createTable(h5pp_table_memory_usage::h5_type, table_path, "memory usage");

    h5pp_table_memory_usage::table mem_usage_entry{};
    mem_usage_entry.iter           = sim_status.iter;
    mem_usage_entry.step           = sim_status.step;
    mem_usage_entry.rss            = tools::common::profile::mem_rss_in_mb();
    mem_usage_entry.hwm            = tools::common::profile::mem_hwm_in_mb();
    mem_usage_entry.vm             = tools::common::profile::mem_vm_in_mb();
    h5ppFile.appendTableEntries(mem_usage_entry, table_path);
}



