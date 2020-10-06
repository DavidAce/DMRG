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
    if(storage_level == StorageLevel::NONE) return;
    // Check if the current entry has already been appended
    // Status is special, flags can be updated without changing iter or step
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    auto                                                                  save_point = std::make_pair(status.iter, status.step);
    tools::log->trace("Appending to table: {}", table_path);
    h5pp_table_algorithm_status::register_table_type();
    if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_algorithm_status::h5_type, table_path, "Algorithm Status");
    tools::common::profile::get_default_prof()["t_hdf"]->tic();
    if(save_log[table_path] == save_point) {
        auto tableInfo = h5ppFile.getTableInfo(table_path);
        h5ppFile.writeTableRecords(status, table_path, tableInfo.numRecords.value() - 1);
    } else
        h5ppFile.appendTableRecords(status, table_path);

    tools::common::profile::get_default_prof()["t_hdf"]->toc();
    save_log[table_path] = save_point;
}

void tools::common::io::h5table::save_profiling(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                const class_algorithm_status &status, std::optional<AlgorithmType> algo_type) {
    if(storage_level == StorageLevel::NONE) return;
    // Check if the current entry has already been appended
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    auto                                                                  save_point = std::make_pair(status.iter, status.step);
    if(save_log[table_path] == save_point) return;

    tools::log->trace("Appending to table: {}", table_path);
    if(not algo_type) algo_type = tools::common::profile::get_current_algo_type();
    switch(algo_type.value()) {
        case(AlgorithmType::xDMRG): {
            h5pp_table_xdmrg_profiling::register_table_type();
            if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_xdmrg_profiling::h5_type, table_path, "xDMRG Profiling");
            h5pp_table_xdmrg_profiling::table profiling_entry;
            /* clang-format off */
            profiling_entry.iter            = status.iter;
            profiling_entry.step            = status.step;
            profiling_entry.position        = status.position;
            profiling_entry.t_tot           = tools::common::profile::t_tot->get_age();
            profiling_entry.t_pre           = tools::common::profile::prof[algo_type.value()]["t_pre"]->get_measured_time();
            profiling_entry.t_rnd           = tools::common::profile::prof[algo_type.value()]["t_rnd"]->get_measured_time();
            profiling_entry.t_pos           = tools::common::profile::prof[algo_type.value()]["t_pos"]->get_measured_time();
            profiling_entry.t_sim           = tools::common::profile::prof[algo_type.value()]["t_sim"]->get_measured_time();
            profiling_entry.t_con           = tools::common::profile::prof[algo_type.value()]["t_con"]->get_measured_time();
            profiling_entry.t_eig           = tools::common::profile::prof[algo_type.value()]["t_eig"]->get_measured_time();
            profiling_entry.t_svd           = tools::common::profile::prof[algo_type.value()]["t_svd"]->get_measured_time();
            profiling_entry.t_env           = tools::common::profile::prof[algo_type.value()]["t_env"]->get_measured_time();
            profiling_entry.t_ent           = tools::common::profile::prof[algo_type.value()]["t_ent"]->get_measured_time();
            profiling_entry.t_ene           = tools::common::profile::prof[algo_type.value()]["t_ene"]->get_measured_time();
            profiling_entry.t_var           = tools::common::profile::prof[algo_type.value()]["t_var"]->get_measured_time();
            profiling_entry.t_prj           = tools::common::profile::prof[algo_type.value()]["t_prj"]->get_measured_time();
            profiling_entry.t_chk           = tools::common::profile::prof[algo_type.value()]["t_chk"]->get_measured_time();
            profiling_entry.t_hdf           = tools::common::profile::prof[algo_type.value()]["t_hdf"]->get_measured_time();
            profiling_entry.t_mps           = tools::common::profile::prof[algo_type.value()]["t_mps"]->get_measured_time();
            profiling_entry.t_mpo           = tools::common::profile::prof[algo_type.value()]["t_mpo"]->get_measured_time();
            profiling_entry.t_opt           = tools::common::profile::prof[algo_type.value()]["t_opt"]->get_measured_time();
            profiling_entry.t_opt_dir       = tools::common::profile::prof[algo_type.value()]["t_opt_dir"]->get_measured_time();
            profiling_entry.t_opt_dir_bfgs  = tools::common::profile::prof[algo_type.value()]["t_opt_dir_bfgs"]->get_measured_time();
            profiling_entry.t_opt_dir_vH2   = tools::common::profile::prof[algo_type.value()]["t_opt_dir_vH2"]->get_measured_time();
            profiling_entry.t_opt_dir_vH2v  = tools::common::profile::prof[algo_type.value()]["t_opt_dir_vH2v"]->get_measured_time();
            profiling_entry.t_opt_dir_vH    = tools::common::profile::prof[algo_type.value()]["t_opt_dir_vH"]->get_measured_time();
            profiling_entry.t_opt_dir_vHv   = tools::common::profile::prof[algo_type.value()]["t_opt_dir_vHv"]->get_measured_time();
            profiling_entry.t_opt_sub       = tools::common::profile::prof[algo_type.value()]["t_opt_sub"]->get_measured_time();
            profiling_entry.t_opt_sub_ham   = tools::common::profile::prof[algo_type.value()]["t_opt_sub_ham"]->get_measured_time();
            profiling_entry.t_opt_sub_hsq   = tools::common::profile::prof[algo_type.value()]["t_opt_sub_hsq"]->get_measured_time();
            profiling_entry.t_opt_sub_lu    = tools::common::profile::prof[algo_type.value()]["t_opt_sub_lu"]->get_measured_time();
            profiling_entry.t_opt_sub_eig   = tools::common::profile::prof[algo_type.value()]["t_opt_sub_eig"]->get_measured_time();
            profiling_entry.t_opt_sub_bfgs  = tools::common::profile::prof[algo_type.value()]["t_opt_sub_bfgs"]->get_measured_time();
            profiling_entry.t_opt_sub_vH2   = tools::common::profile::prof[algo_type.value()]["t_opt_sub_vH2"]->get_measured_time();
            profiling_entry.t_opt_sub_vH2v  = tools::common::profile::prof[algo_type.value()]["t_opt_sub_vH2v"]->get_measured_time();
            profiling_entry.t_opt_sub_vH    = tools::common::profile::prof[algo_type.value()]["t_opt_sub_vH"]->get_measured_time();
            profiling_entry.t_opt_sub_vHv   = tools::common::profile::prof[algo_type.value()]["t_opt_sub_vHv"]->get_measured_time();
            /* clang-format on */
            tools::common::profile::get_default_prof()["t_hdf"]->tic();
            h5ppFile.appendTableRecords(profiling_entry, table_path);
            tools::common::profile::get_default_prof()["t_hdf"]->toc();
            break;
        }
        case(AlgorithmType::fDMRG): {
            h5pp_table_fdmrg_profiling::register_table_type();
            if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_fdmrg_profiling::h5_type, table_path, "fDMRG Profiling");
            h5pp_table_fdmrg_profiling::table profiling_entry;
            /* clang-format off */
            profiling_entry.iter            = status.iter;
            profiling_entry.step            = status.step;
            profiling_entry.position        = status.position;
            profiling_entry.t_tot           = tools::common::profile::t_tot->get_age();
            profiling_entry.t_pre           = tools::common::profile::prof[algo_type.value()]["t_pre"]->get_measured_time();
            profiling_entry.t_rnd           = tools::common::profile::prof[algo_type.value()]["t_rnd"]->get_measured_time();
            profiling_entry.t_pos           = tools::common::profile::prof[algo_type.value()]["t_pos"]->get_measured_time();
            profiling_entry.t_sim           = tools::common::profile::prof[algo_type.value()]["t_sim"]->get_measured_time();
            profiling_entry.t_con           = tools::common::profile::prof[algo_type.value()]["t_con"]->get_measured_time();
            profiling_entry.t_eig           = tools::common::profile::prof[algo_type.value()]["t_eig"]->get_measured_time();
            profiling_entry.t_svd           = tools::common::profile::prof[algo_type.value()]["t_svd"]->get_measured_time();
            profiling_entry.t_env           = tools::common::profile::prof[algo_type.value()]["t_env"]->get_measured_time();
            profiling_entry.t_ent           = tools::common::profile::prof[algo_type.value()]["t_ent"]->get_measured_time();
            profiling_entry.t_ene           = tools::common::profile::prof[algo_type.value()]["t_ene"]->get_measured_time();
            profiling_entry.t_var           = tools::common::profile::prof[algo_type.value()]["t_var"]->get_measured_time();
            profiling_entry.t_chk           = tools::common::profile::prof[algo_type.value()]["t_chk"]->get_measured_time();
            profiling_entry.t_hdf           = tools::common::profile::prof[algo_type.value()]["t_hdf"]->get_measured_time();
            profiling_entry.t_mps           = tools::common::profile::prof[algo_type.value()]["t_mps"]->get_measured_time();
            profiling_entry.t_mpo           = tools::common::profile::prof[algo_type.value()]["t_mpo"]->get_measured_time();
            /* clang-format on */
            tools::common::profile::get_default_prof()["t_hdf"]->tic();
            h5ppFile.appendTableRecords(profiling_entry, table_path);
            tools::common::profile::get_default_prof()["t_hdf"]->toc();
            break;
        }
        case(AlgorithmType::iDMRG): {
            h5pp_table_idmrg_profiling::register_table_type();
            if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_idmrg_profiling::h5_type, table_path, "Profiling iDMRG");
            h5pp_table_idmrg_profiling::table profiling_entry;
            /* clang-format off */
            profiling_entry.iter            = status.iter;
            profiling_entry.step            = status.step;
            profiling_entry.position        = status.position;
            profiling_entry.t_tot           = tools::common::profile::t_tot->get_age();
            profiling_entry.t_pre           = tools::common::profile::prof[algo_type.value()]["t_pre"]->get_measured_time();
            profiling_entry.t_rnd           = tools::common::profile::prof[algo_type.value()]["t_rnd"]->get_measured_time();
            profiling_entry.t_pos           = tools::common::profile::prof[algo_type.value()]["t_pos"]->get_measured_time();
            profiling_entry.t_sim           = tools::common::profile::prof[algo_type.value()]["t_sim"]->get_measured_time();
            profiling_entry.t_con           = tools::common::profile::prof[algo_type.value()]["t_con"]->get_measured_time();
            profiling_entry.t_eig           = tools::common::profile::prof[algo_type.value()]["t_eig"]->get_measured_time();
            profiling_entry.t_svd           = tools::common::profile::prof[algo_type.value()]["t_svd"]->get_measured_time();
            profiling_entry.t_env           = tools::common::profile::prof[algo_type.value()]["t_env"]->get_measured_time();
            profiling_entry.t_ent           = tools::common::profile::prof[algo_type.value()]["t_ent"]->get_measured_time();
            profiling_entry.t_ene           = tools::common::profile::prof[algo_type.value()]["t_ene"]->get_measured_time();
            profiling_entry.t_var           = tools::common::profile::prof[algo_type.value()]["t_var"]->get_measured_time();
            profiling_entry.t_chk           = tools::common::profile::prof[algo_type.value()]["t_chk"]->get_measured_time();
            profiling_entry.t_hdf           = tools::common::profile::prof[algo_type.value()]["t_hdf"]->get_measured_time();
            profiling_entry.t_ene_ham       = tools::common::profile::prof[algo_type.value()]["t_ene_ham"]->get_measured_time();
            profiling_entry.t_ene_mom       = tools::common::profile::prof[algo_type.value()]["t_ene_mom"]->get_measured_time();
            profiling_entry.t_var_ham       = tools::common::profile::prof[algo_type.value()]["t_var_ham"]->get_measured_time();
            profiling_entry.t_var_mom       = tools::common::profile::prof[algo_type.value()]["t_var_mom"]->get_measured_time();
            /* clang-format on */
            tools::common::profile::get_default_prof()["t_hdf"]->tic();
            h5ppFile.appendTableRecords(profiling_entry, table_path);
            tools::common::profile::get_default_prof()["t_hdf"]->toc();
            break;
        }
        case(AlgorithmType::iTEBD): {
            h5pp_table_itebd_profiling::register_table_type();
            if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_itebd_profiling::h5_type, table_path, "iTEBD Profiling");
            h5pp_table_itebd_profiling::table profiling_entry;
            /* clang-format off */
            profiling_entry.iter            = status.iter;
            profiling_entry.step            = status.step;
            profiling_entry.position        = status.position;
            profiling_entry.t_tot           = tools::common::profile::t_tot->get_age();
            profiling_entry.t_pre           = tools::common::profile::prof[algo_type.value()]["t_pre"]->get_measured_time();
            profiling_entry.t_pos           = tools::common::profile::prof[algo_type.value()]["t_pos"]->get_measured_time();
            profiling_entry.t_sim           = tools::common::profile::prof[algo_type.value()]["t_sim"]->get_measured_time();
            profiling_entry.t_con           = tools::common::profile::prof[algo_type.value()]["t_con"]->get_measured_time();
            profiling_entry.t_svd           = tools::common::profile::prof[algo_type.value()]["t_svd"]->get_measured_time();
            profiling_entry.t_evo           = tools::common::profile::prof[algo_type.value()]["t_evo"]->get_measured_time();
            profiling_entry.t_ent           = tools::common::profile::prof[algo_type.value()]["t_ent"]->get_measured_time();
            profiling_entry.t_chk           = tools::common::profile::prof[algo_type.value()]["t_chk"]->get_measured_time();
            profiling_entry.t_hdf           = tools::common::profile::prof[algo_type.value()]["t_hdf"]->get_measured_time();
            profiling_entry.t_ene_ham       = tools::common::profile::prof[algo_type.value()]["t_ene_ham"]->get_measured_time();
            profiling_entry.t_ene_mom       = tools::common::profile::prof[algo_type.value()]["t_ene_mom"]->get_measured_time();
            profiling_entry.t_var_ham       = tools::common::profile::prof[algo_type.value()]["t_var_ham"]->get_measured_time();
            profiling_entry.t_var_mom       = tools::common::profile::prof[algo_type.value()]["t_var_mom"]->get_measured_time();
            /* clang-format on */
            tools::common::profile::get_default_prof()["t_hdf"]->tic();
            h5ppFile.appendTableRecords(profiling_entry, table_path);
            tools::common::profile::get_default_prof()["t_hdf"]->toc();
            break;
        }
    }
    save_log[table_path] = save_point;
}

void tools::common::io::h5table::save_mem_usage(h5pp::File &h5ppFile, const std::string &table_path, const StorageLevel &storage_level,
                                                const class_algorithm_status &status) {
    if(storage_level == StorageLevel::NONE) return;
    // Check if the current entry has already been appended
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    auto                                                                  save_point = std::make_pair(status.iter, status.step);
    if(save_log[table_path] == save_point) return;
    log->trace("Appending to table: {}", table_path);

    h5pp_table_memory_usage::register_table_type();
    if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_memory_usage::h5_type, table_path, "memory usage");

    h5pp_table_memory_usage::table mem_usage_entry{};
    mem_usage_entry.iter = status.iter;
    mem_usage_entry.step = status.step;
    mem_usage_entry.rss  = tools::common::profile::mem_rss_in_mb();
    mem_usage_entry.hwm  = tools::common::profile::mem_hwm_in_mb();
    mem_usage_entry.vm   = tools::common::profile::mem_vm_in_mb();
    tools::common::profile::get_default_prof()["t_hdf"]->tic();
    h5ppFile.appendTableRecords(mem_usage_entry, table_path);
    tools::common::profile::get_default_prof()["t_hdf"]->toc();
    save_log[table_path] = save_point;
}

void tools::common::io::h5table::load_sim_status(const h5pp::File &h5ppFile, const std::string &state_prefix, class_algorithm_status &status) {
    tools::common::profile::get_default_prof()["t_hdf"]->tic();
    std::string table_path = fmt::format("{}/status", state_prefix);
    if(h5ppFile.linkExists(table_path)) {
        tools::log->info("Loading status from table: [{}]", table_path);
        h5ppFile.readTableRecords(status, table_path); // Reads the last entry by default
    } else {
        throw std::runtime_error(
            fmt::format("Could not find table [status] in file [{}] at prefix [{}] at path [{}]", h5ppFile.getFilePath(), state_prefix, table_path));
    }
    tools::common::profile::get_default_prof()["t_hdf"]->toc();
}

void tools::common::io::h5table::load_profiling(const h5pp::File &h5ppFile, const std::string &state_prefix, std::optional<AlgorithmType> algo_type) {
    if(not settings::profiling::on) return;
    if(not algo_type) algo_type = tools::common::profile::get_current_algo_type();
    std::string table_path = state_prefix + std::string("/profiling");
    if(h5ppFile.linkExists(table_path)) tools::log->info("Loading profiling from table: [{}]", table_path);
    else
        throw std::runtime_error(
            fmt::format("Could not find table [profiling] in file [{}] at prefix [{}] at path [{}]", h5ppFile.getFilePath(), state_prefix, table_path));

    auto table_info = h5ppFile.getTableInfo(table_path);
    for(const auto &t_name : table_info.fieldNames.value()) {
        if(tools::common::profile::prof[algo_type.value()].find(t_name) == tools::common::profile::prof[algo_type.value()].end()) continue;
        *tools::common::profile::prof[algo_type.value()][t_name] = h5ppFile.readTableField<double>(table_path, t_name, h5pp::TableSelection::LAST);
    }

    *tools::common::profile::t_tot = h5ppFile.readTableField<double>(table_path, "t_tot", h5pp::TableSelection::LAST);
}
