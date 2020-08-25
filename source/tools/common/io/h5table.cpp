//
// Created by david on 2019-11-07.
//
#include <config/enums.h>
#include <config/nmspc_settings.h>
#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

void tools::common::io::h5table::save_sim_status(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                 const class_algorithm_status &status) {
    if(storage_level < StorageLevel::LIGHT) return;
    tools::log->trace("Appending to table: {}", table_path);
    h5pp_table_algorithm_status::register_table_type();
    if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_algorithm_status::h5_type, table_path, "Algorithm Status");
    h5ppFile.appendTableRecords(status, table_path);
}

void tools::common::io::h5table::save_profiling(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                const class_algorithm_status &status) {
    if(storage_level < StorageLevel::LIGHT) return;
    tools::log->trace("Appending to table: {}", table_path);

    h5pp_table_profiling::register_table_type();
    if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_profiling::h5_type, table_path, "Profiling");

    h5pp_table_profiling::table profiling_entry;
    /* clang-format off */
    profiling_entry.iter            = status.iter;
    profiling_entry.step            = status.step;
    profiling_entry.position        = status.position;
    profiling_entry.t_tot           = tools::common::profile::t_tot->get_age();
    profiling_entry.t_pre           = tools::common::profile::t_pre->get_measured_time();
    profiling_entry.t_pos           = tools::common::profile::t_pos->get_measured_time();
    profiling_entry.t_sim           = tools::common::profile::t_sim->get_measured_time();
    profiling_entry.t_con           = tools::common::profile::t_con->get_measured_time();
    profiling_entry.t_eig           = tools::common::profile::t_eig->get_measured_time();
    profiling_entry.t_svd           = tools::common::profile::t_svd->get_measured_time();
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
    profiling_entry.t_mps           = tools::common::profile::t_mps->get_measured_time();
    profiling_entry.t_mpo           = tools::common::profile::t_mpo->get_measured_time();
    profiling_entry.t_opt           = tools::common::profile::t_opt->get_measured_time();
    profiling_entry.t_opt_dir       = tools::common::profile::t_opt_dir->get_measured_time();
    profiling_entry.t_opt_dir_bfgs  = tools::common::profile::t_opt_dir_bfgs->get_measured_time();
    profiling_entry.t_opt_dir_vH2   = tools::common::profile::t_opt_dir_vH2->get_measured_time();
    profiling_entry.t_opt_dir_vH2v  = tools::common::profile::t_opt_dir_vH2v->get_measured_time();
    profiling_entry.t_opt_dir_vH    = tools::common::profile::t_opt_dir_vH->get_measured_time();
    profiling_entry.t_opt_dir_vHv   = tools::common::profile::t_opt_dir_vHv->get_measured_time();
    profiling_entry.t_opt_sub       = tools::common::profile::t_opt_sub->get_measured_time();
    profiling_entry.t_opt_sub_ham   = tools::common::profile::t_opt_sub_ham->get_measured_time();
    profiling_entry.t_opt_sub_hsq   = tools::common::profile::t_opt_sub_hsq->get_measured_time();
    profiling_entry.t_opt_sub_lu    = tools::common::profile::t_opt_sub_lu->get_measured_time();
    profiling_entry.t_opt_sub_eig   = tools::common::profile::t_opt_sub_eig->get_measured_time();
    profiling_entry.t_opt_sub_bfgs  = tools::common::profile::t_opt_sub_bfgs->get_measured_time();
    profiling_entry.t_opt_sub_vH2   = tools::common::profile::t_opt_sub_vH2->get_measured_time();
    profiling_entry.t_opt_sub_vH2v  = tools::common::profile::t_opt_sub_vH2v->get_measured_time();
    profiling_entry.t_opt_sub_vH    = tools::common::profile::t_opt_sub_vH->get_measured_time();
    profiling_entry.t_opt_sub_vHv   = tools::common::profile::t_opt_sub_vHv->get_measured_time();
    /* clang-format on */
    h5ppFile.appendTableRecords(profiling_entry, table_path);
}

void tools::common::io::h5table::save_mem_usage(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                const class_algorithm_status &status) {
    if(storage_level < StorageLevel::LIGHT) return;
    log->trace("Appending to table: {}", table_path);

    h5pp_table_memory_usage::register_table_type();
    if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_memory_usage::h5_type, table_path, "memory usage");

    h5pp_table_memory_usage::table mem_usage_entry{};
    mem_usage_entry.iter = status.iter;
    mem_usage_entry.step = status.step;
    mem_usage_entry.rss  = tools::common::profile::mem_rss_in_mb();
    mem_usage_entry.hwm  = tools::common::profile::mem_hwm_in_mb();
    mem_usage_entry.vm   = tools::common::profile::mem_vm_in_mb();
    h5ppFile.appendTableRecords(mem_usage_entry, table_path);
}

void tools::common::io::h5table::load_sim_status(const h5pp::File &h5ppFile, const std::string &state_prefix, class_algorithm_status &status) {
    tools::common::profile::t_hdf->tic();
    std::string table_path = fmt::format("{}/status",state_prefix);
    if(h5ppFile.linkExists(table_path)) {
        tools::log->info("Loading status from table: [{}]", table_path);
        h5ppFile.readTableRecords(status, table_path); // Reads the last entry by default
    }else{
        throw std::runtime_error(fmt::format("Could not find table [status] in file [{}] at prefix [{}] at path [{}]", h5ppFile.getFilePath(), state_prefix, table_path));
    }
    tools::common::profile::t_hdf->toc();
}

void tools::common::io::h5table::load_profiling(const h5pp::File &h5ppFile, const std::string &state_prefix) {
    if(not settings::profiling::on) return;
    h5pp_table_profiling::table prof_entry;
    tools::common::profile::t_hdf->tic();
    std::string table_path = state_prefix + std::string("/profiling");
    if(h5ppFile.linkExists(table_path)) {
        tools::log->info("Loading profiling from table: [{}]", table_path);
        h5ppFile.readTableRecords(prof_entry, table_path); // Reads the last entry by default
    }else
        throw std::runtime_error(fmt::format("Could not find table [profiling] in file [{}] at prefix [{}] at path [{}]", h5ppFile.getFilePath(),state_prefix, table_path));
    tools::common::profile::t_hdf->toc();

    /* clang-format off */
    *tools::common::profile::t_tot          = prof_entry.t_tot;
    *tools::common::profile::t_pre          = prof_entry.t_pre;
    *tools::common::profile::t_pos          = prof_entry.t_pos;
    *tools::common::profile::t_sim          = prof_entry.t_sim;
    *tools::common::profile::t_con          = prof_entry.t_con;
    *tools::common::profile::t_eig          = prof_entry.t_eig;
    *tools::common::profile::t_svd          = prof_entry.t_svd;
    *tools::common::profile::t_evo          = prof_entry.t_evo;
    *tools::common::profile::t_env          = prof_entry.t_env;
    *tools::common::profile::t_ent          = prof_entry.t_ent;
    *tools::common::profile::t_ene          = prof_entry.t_ene;
    *tools::common::profile::t_var          = prof_entry.t_var;
    *tools::common::profile::t_prj          = prof_entry.t_prj;
    *tools::common::profile::t_chk          = prof_entry.t_chk;
    *tools::common::profile::t_hdf          = prof_entry.t_hdf;
    *tools::common::profile::t_ene_ham      = prof_entry.t_ene_ham;
    *tools::common::profile::t_ene_mom      = prof_entry.t_ene_mom;
    *tools::common::profile::t_var_ham      = prof_entry.t_var_ham;
    *tools::common::profile::t_var_mom      = prof_entry.t_var_mom;
    *tools::common::profile::t_mps          = prof_entry.t_mps;
    *tools::common::profile::t_mpo          = prof_entry.t_mpo;
    *tools::common::profile::t_opt          = prof_entry.t_opt;
    *tools::common::profile::t_opt_dir      = prof_entry.t_opt_dir;
    *tools::common::profile::t_opt_dir_bfgs = prof_entry.t_opt_dir_bfgs;
    *tools::common::profile::t_opt_dir_vH2  = prof_entry.t_opt_dir_vH2;
    *tools::common::profile::t_opt_dir_vH2v = prof_entry.t_opt_dir_vH2v;
    *tools::common::profile::t_opt_dir_vH   = prof_entry.t_opt_dir_vH;
    *tools::common::profile::t_opt_dir_vHv  = prof_entry.t_opt_dir_vHv;
    *tools::common::profile::t_opt_sub      = prof_entry.t_opt_sub;
    *tools::common::profile::t_opt_sub_ham  = prof_entry.t_opt_sub_ham;
    *tools::common::profile::t_opt_sub_hsq  = prof_entry.t_opt_sub_hsq;
    *tools::common::profile::t_opt_sub_lu   = prof_entry.t_opt_sub_lu;
    *tools::common::profile::t_opt_sub_eig  = prof_entry.t_opt_sub_eig;
    *tools::common::profile::t_opt_sub_bfgs = prof_entry.t_opt_sub_bfgs;
    *tools::common::profile::t_opt_sub_vH2  = prof_entry.t_opt_sub_vH2;
    *tools::common::profile::t_opt_sub_vH2v = prof_entry.t_opt_sub_vH2v;
    *tools::common::profile::t_opt_sub_vH   = prof_entry.t_opt_sub_vH;
    *tools::common::profile::t_opt_sub_vHv  = prof_entry.t_opt_sub_vHv;
    /* clang-format off */
}
