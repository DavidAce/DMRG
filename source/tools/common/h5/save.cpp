
#include <algorithms/AlgorithmStatus.h>
#include <config/enums.h>
#include <config/settings.h>
#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <string>
#include <tid/tid.h>
#include <tools/common/h5.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

namespace tools::common::h5 {

    void save::bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile,
                                  const std::string &link) {
        // from h5table
        if(save_log.empty()) {
            try {
                if(h5ppFile.linkExists(link)) {
                    auto step      = h5ppFile.readAttribute<uint64_t>("step", link);
                    auto iter      = h5ppFile.readAttribute<uint64_t>("iteration", link);
                    save_log[link] = std::make_pair(iter, step);
                }
            } catch(const std::exception &ex) { tools::log->warn("Could not bootstrap save_log: {}", ex.what()); }
        }
    }

    void save::bootstrap_meta_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile,
                                  const std::string &state_prefix) {
        if(save_log.empty()) {
            try {
                uint64_t step = 0;
                uint64_t iter = 0;
                if(h5ppFile.linkExists("common/step")) step = h5ppFile.readAttribute<uint64_t>(state_prefix, "common/step");
                if(h5ppFile.linkExists("common/iteration")) iter = h5ppFile.readAttribute<uint64_t>(state_prefix, "common/iteration");
                save_log[state_prefix] = std::make_pair(iter, step);

            } catch(const std::exception &ex) { tools::log->warn("Could not bootstrap save_log for {}: {}", state_prefix, ex.what()); }
        }
    }

    template<typename AttrType>
    void save::attr(h5pp::File &file, const AttrType &attrData, const std::string &attrName, const std::string &linkPath, const std::string &linkText) {
        tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", attrName, attrData, linkPath);
        if(not file.linkExists(linkPath)) file.writeDataset(linkText, linkPath);
        file.writeAttribute(attrData, attrName, linkPath);
    }

    void save::status(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level, const AlgorithmStatus &status) {
        if(storage_level == StorageLevel::NONE) return;
        // Check if the current entry has already been appended
        // Status is special, flags can be updated without changing iter or step
        std::string                                                           table_path = fmt::format("{}/status", table_prefix);
        auto                                                                  t_hdf      = tid::tic_scope("status");
        static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
        bootstrap_save_log(save_log, h5ppFile, table_path);
        auto save_point = std::make_pair(status.iter, status.step);
        tools::log->trace("Appending to table: {}", table_path);
        h5pp_table_algorithm_status::register_table_type();
        if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_algorithm_status::h5_type, table_path, "Algorithm Status");
        if(save_log[table_path] == save_point) {
            // The table has been saved at this iteration, so we overwrite the last entry.
            auto tableInfo = h5ppFile.getTableInfo(table_path);
            h5ppFile.writeTableRecords(status, table_path, tableInfo.numRecords.value() - 1);
        } else
            h5ppFile.appendTableRecords(status, table_path);
        h5ppFile.writeAttribute(status.iter, "iteration", table_path);
        h5ppFile.writeAttribute(status.step, "step", table_path);

        save_log[table_path] = save_point;
    }

    void save::mem(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level, const AlgorithmStatus &status) {
        if(storage_level == StorageLevel::NONE) return;
        // Check if the current entry has already been appended
        std::string                                                           table_path = fmt::format("{}/mem_usage", table_prefix);
        static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
        bootstrap_save_log(save_log, h5ppFile, table_path);
        auto save_point = std::make_pair(status.iter, status.step);
        if(save_log[table_path] == save_point) return;
        log->trace("Appending to table: {}", table_path);
        auto t_mem = tid::tic_scope("mem");
        h5pp_table_memory_usage::register_table_type();
        if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_memory_usage::h5_type, table_path, "memory usage");

        h5pp_table_memory_usage::table mem_usage_entry{};
        mem_usage_entry.iter = status.iter;
        mem_usage_entry.step = status.step;
        mem_usage_entry.rss  = tools::common::profile::mem_rss_in_mb();
        mem_usage_entry.hwm  = tools::common::profile::mem_hwm_in_mb();
        mem_usage_entry.vm   = tools::common::profile::mem_vm_in_mb();
        h5ppFile.appendTableRecords(mem_usage_entry, table_path);
        h5ppFile.writeAttribute(status.iter, "iteration", table_path);
        h5ppFile.writeAttribute(status.step, "step", table_path);
        save_log[table_path] = save_point;
    }

    void save::meta(h5pp::File &h5ppFile, const StorageLevel &storage_level, const StorageReason &storage_reason, const ModelType &model_type,
                    size_t model_size, const AlgorithmType &algo_type, const std::string &state_name, const std::string &state_prefix,
                    const std::string &model_prefix, const AlgorithmStatus &status) {
        // Write metadata into /common so that we can find state paths later.
        // Under /common we have the following dummy datasets defining attributes which map state_prefix to useful information:
        //
        // -- finished       | true or 1 if state_prefix is finished
        // -- storage_level  | state_prefix -> storage_level
        // -- storage_reason | state_prefix -> storage_reason
        // -- state_root     | state_prefix -> state_root
        // -- hamiltonian    | state_prefix -> path to hamiltonian table
        // -- mps_prefix     | state_prefix -> mps_prefix
        // -- mpo_prefix     | state_prefix -> mpo_prefix
        // -- model_type     | state_prefix -> model_type
        // -- model_size     | state_prefix -> model_size
        // -- algo_type      | state_prefix -> algo_type
        // -- state_name     | state_prefix -> state_name
        // -- iteration      | state_prefix -> iteration
        // -- step           | state_prefix -> step
        // -- position       | state_prefix -> position of the mps

        // Checks if the current entries have already been written
        auto                                                                  t_meta = tid::tic_scope("meta");
        static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
        bootstrap_meta_log(save_log, h5ppFile, state_prefix);

        auto save_point = std::make_pair(status.iter, status.step);
        if(save_log[state_prefix] == save_point) return;

        std::string storage_level_str(enum2str(storage_level));
        std::string storage_reason_str(enum2str(storage_reason));
        std::string model_name(enum2str(model_type));
        std::string algo_name(enum2str(algo_type));
        std::string state_root  = algo_name + '/' + state_name;
        std::string hamiltonian = model_prefix + "/hamiltonian";
        std::string mpo_prefix  = model_prefix + "/mpo";
        std::string mps_prefix  = state_prefix + "/mps";

        save::attr(h5ppFile, status.algorithm_has_finished, state_prefix, "common/finished", "Maps state_prefix -> finished");
        save::attr(h5ppFile, storage_level_str, state_prefix, "common/storage_level", "Maps state_prefix -> storage_level");
        save::attr(h5ppFile, storage_reason_str, state_prefix, "common/storage_reason", "Maps state_prefix -> storage_reason");
        save::attr(h5ppFile, state_root, state_prefix, "common/state_root", "Maps state_prefix -> state_root");
        save::attr(h5ppFile, model_prefix, state_prefix, "common/model_prefix", "Maps state_prefix -> model_prefix");
        save::attr(h5ppFile, hamiltonian, state_prefix, "common/hamiltonian", "Maps state_prefix -> hamiltonian");
        save::attr(h5ppFile, mpo_prefix, state_prefix, "common/mpo_prefix", "Maps state_prefix -> mpo_prefix");
        save::attr(h5ppFile, mps_prefix, state_prefix, "common/mps_prefix", "Maps state_prefix -> mps_prefix");
        save::attr(h5ppFile, model_name, state_prefix, "common/model_type", "Maps state_prefix -> model_type");
        save::attr(h5ppFile, model_size, state_prefix, "common/model_size", "Maps state_prefix -> model_size");
        save::attr(h5ppFile, algo_name, state_prefix, "common/algo_type", "Maps state_prefix -> algo_type");
        save::attr(h5ppFile, state_name, state_prefix, "common/state_name", "Maps state_prefix -> state_name");
        save::attr(h5ppFile, status.iter, state_prefix, "common/iteration", "Maps state_prefix -> iteration");
        save::attr(h5ppFile, status.step, state_prefix, "common/step", "Maps state_prefix -> step");
        save::attr(h5ppFile, status.position, state_prefix, "common/position", "Maps state_prefix -> position");
        save_log[state_prefix] = save_point;
    }

    void save::timer(h5pp::File &h5ppFile, const std::string &state_prefix, const StorageLevel &storage_level, const AlgorithmStatus &status) {
        if(not settings::profiling::on) return;
        if(storage_level == StorageLevel::NONE) return;
        // Check if the current entry has already been appended
        static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
        bootstrap_save_log(save_log, h5ppFile, state_prefix);
        auto save_point = std::make_pair(status.iter, status.step);
        if(save_log[state_prefix] == save_point) return;
        auto prof_prefix = fmt::format("{}/profiling", state_prefix);

        tools::log->trace("Appending profiling to: {}", prof_prefix);
        auto t_prof = tid::tic_token("prof");

        h5pp_ur::register_table_type();
        h5ppFile.setKeepFileOpened();
        for(const auto &t : tid::search(status.algo_type_sv())) {
            auto dsetname = fmt::format("{}/{}", prof_prefix, t.key);
            auto dsetdata = h5pp_ur::item{t->get_time(), t->get_time_avg(), t->get_tic_count()};
            h5ppFile.writeDataset(dsetdata, dsetname, h5pp_ur::h5_type);
        }
        h5ppFile.setKeepFileClosed();
        save_log[state_prefix] = save_point;
    }

}
