
#include "algorithms/AlgorithmStatus.h"
#include "config/enums.h"
#include "config/settings.h"
#include "debug/info.h"
#include "general/iter.h"
#include "io/hdf5_types.h"
#include "tid/tid.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include <h5pp/h5pp.h>
#include <string>

namespace tools::common::h5 {
    void save::bootstrap_meta_log(std::unordered_map<std::string, AlgorithmStatus> &save_log, const h5pp::File &h5file, std::string_view state_prefix) {
        if(save_log.empty()) {
            try {
                AlgorithmStatus status;
                if(h5file.linkExists("common/status")) status = h5file.readAttribute<AlgorithmStatus>("common/status", state_prefix);
                save_log.insert(std::make_pair(state_prefix, status));
            } catch(const std::exception &ex) { tools::log->warn("Could not bootstrap save_log for {}: {}", state_prefix, ex.what()); }
        }
    }

    template<typename AttrType>
    void save::attr(h5pp::File &h5file, const AttrType &attrData, std::string_view linkPath, std::string_view attrName, std::string_view linkText,
                    std::optional<h5pp::hid::h5t> h5type) {
        if(not h5file.linkExists(attrName)) {
            tools::log->trace("link [{}] does not exist. Returning ...", attrName);
            return; // This means that there is nothing written to the state_prefix in attrName.
        }
        if(not h5file.linkExists(linkPath)) h5file.writeDataset(linkText, linkPath);
        tools::log->trace("Link {:<32} | attribute -- {: <40}", linkPath, attrName);
        h5file.writeAttribute(attrData, linkPath, attrName, std::nullopt, h5type);
    }

    void save::status(h5pp::File &h5file, const StorageInfo &sinfo, const AlgorithmStatus &status) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        // Check if the current entry has already been appended
        // Status is special, flags can be updated without changing iter or step
        std::string table_path    = fmt::format("{}/status", sinfo.get_state_prefix());
        auto        t_hdf         = tid::tic_scope("status", tid::level::detailed);
        auto        h5_save_point = save::get_last_save_point(h5file, table_path);
        auto        save_point    = std::make_pair(status.iter, status.step);
        h5pp_table_algorithm_status::register_table_type();
        if(not h5_save_point) h5file.createTable(h5pp_table_algorithm_status::h5_type, table_path, "Algorithm Status");
        auto info   = h5file.getTableInfo(table_path);
        auto offset = h5_save_point and h5_save_point.value() == save_point ? info.numRecords.value() - 1 : info.numRecords.value();
        h5pp::hdf5::writeTableRecords(status, info, offset, 1);
        h5file.writeAttribute(status.iter, table_path, "iter");
        h5file.writeAttribute(status.step, table_path, "step");
        h5file.writeAttribute(status.bond_lim, table_path, "bond_lim");
        h5file.writeAttribute(status.bond_max, table_path, "bond_max");
    }

    void save::mem(h5pp::File &h5file, const StorageInfo &sinfo) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        // Check if the current entry has already been appended
        std::string                    table_path = fmt::format("{}/mem_usage", sinfo.get_state_prefix());
        h5pp_table_memory_usage::table mem_usage_entry{};
        mem_usage_entry.iter     = sinfo.iter;
        mem_usage_entry.step     = sinfo.step;
        mem_usage_entry.bond_lim = sinfo.bond_lim;
        mem_usage_entry.rss      = debug::mem_rss_in_mb();
        mem_usage_entry.hwm      = debug::mem_hwm_in_mb();
        mem_usage_entry.vm       = debug::mem_vm_in_mb();

        auto h5_save_point = save::get_last_save_point(h5file, table_path);
        auto save_point    = std::make_pair(sinfo.iter, sinfo.step);
        h5pp_table_memory_usage::register_table_type();
        if(not h5_save_point) h5file.createTable(h5pp_table_memory_usage::h5_type, table_path, "Memory usage");
        auto info   = h5file.getTableInfo(table_path);
        auto offset = h5_save_point and h5_save_point.value() == save_point ? info.numRecords.value() - 1 : info.numRecords.value();
        h5pp::hdf5::writeTableRecords(mem_usage_entry, info, offset, 1);
        h5file.writeAttribute(sinfo.iter, table_path, "iter");
        h5file.writeAttribute(sinfo.step, table_path, "step");
        h5file.writeAttribute(sinfo.bond_lim, table_path, "bond_lim");
        h5file.writeAttribute(sinfo.bond_max, table_path, "bond_max");
    }

    std::optional<std::string> find::find_duplicate_save(const h5pp::File &h5file, std::string_view state_prefix, const AlgorithmStatus &status) {
        if(not h5file.linkExists("common/status")) return std::nullopt;
        if(not h5file.linkExists("common/state_root")) return std::nullopt;
        h5pp::Options opt_status;
        opt_status.h5Type   = h5pp_table_algorithm_status::h5_type;
        opt_status.linkPath = "common/status";
        AlgorithmStatus h5_status;
        for(const auto &h5_state_prefix : h5file.getAttributeNames("common/status")) {
            if(state_prefix == h5_state_prefix) continue; // Skip self
            auto h5_state_root_exists = h5pp::hdf5::checkIfAttrExists(h5file.openFileHandle(), "common/state_root", h5_state_prefix);
            if(not h5_state_root_exists) return std::nullopt;

            auto h5_state_root = h5file.readAttribute<std::string>("common/state_root", h5_state_prefix);
            if(state_prefix.find(h5_state_root) != std::string_view::npos) {
                // We have a match! A state on file with the same root as the given state prefix
                // A state_prefix is something like xDMRG/state_0/checkpoint/iter_last
                // A state_root   is something like xDMRG/state_0
                // If roots match they are trying to describe the same state.

                // Now read the status of the found state prefix
                opt_status.attrName = h5_state_prefix;
                h5file.readAttribute(h5_status, opt_status);

                if(status.algo_type != h5_status.algo_type) continue;
                if(status.iter != h5_status.iter) continue;
                if(status.step != h5_status.step) continue;
                if(status.position != h5_status.position) continue;
                if(status.bond_lim != h5_status.bond_lim) continue;
                // We have a match: Saving state_prefix would create a duplicate of h5_state_prefix!
                return h5_state_prefix;
            }
        }
        return std::nullopt;
    }

    std::optional<std::pair<uint64_t, uint64_t>> save::get_last_save_point(const h5pp::File &h5file, std::string_view link_path) {
        auto link_info = h5file.getLinkInfo(link_path);
        if(link_info.linkExists.value()) {
            auto iter_exists = h5pp::hdf5::checkIfAttrExists(link_info.getLocId(), link_path, "iter");
            auto step_exists = h5pp::hdf5::checkIfAttrExists(link_info.getLocId(), link_path, "step");
            if(iter_exists and step_exists) {
                auto iter = h5file.readAttribute<uint64_t>(link_path, "iter");
                auto step = h5file.readAttribute<uint64_t>(link_path, "step");
                return std::make_pair(iter, step);
            }
        }
        return std::nullopt;
    }

    //    void save::meta(h5pp::File &h5file, const StorageLevel &storage_level, const StorageEvent &storage_reason, const ModelType &model_type, size_t
    //    model_size,
    //                    std::string_view state_name, std::string_view state_prefix, std::string_view model_prefix, const AlgorithmStatus &status) {
    //        // Write metadata into /common so that we can find state paths later.
    //        // Under /common we have the following dummy datasets defining attributes which map state_prefix to useful information:
    //        //
    //        // -- finished       | true or 1 if state_prefix is finished
    //        // -- storage_level  | state_prefix -> storage_level
    //        // -- storage_reason | state_prefix -> storage_reason
    //        // -- state_root     | state_prefix -> state_root
    //        // -- hamiltonian    | state_prefix -> path to hamiltonian table
    //        // -- table_prefix   | state_prefix -> path to tables. There can be multiple
    //        // -- timer_prefix   | state_prefix -> path to timer data
    //        // -- mps_prefix     | state_prefix -> mps_prefix
    //        // -- mpo_prefix     | state_prefix -> mpo_prefix
    //        // -- model_type     | state_prefix -> model_type
    //        // -- model_size     | state_prefix -> model_size
    //        // -- algo_type      | state_prefix -> algo_type
    //        // -- state_name     | state_prefix -> state_name
    //        // -- iteration      | state_prefix -> iteration
    //        // -- step           | state_prefix -> step
    //        // -- position       | state_prefix -> position of the mps
    //
    //        if(storage_level == StorageLevel::NONE) return;
    //        if(storage_reason == StorageEvent::MODEL) return;
    //        //        if(not h5file.linkExists(state_prefix)) return; // No point in saving metadata for non-existing state prefixes
    //        tools::log->trace("Writing metadata to /common");
    //        auto t_meta = tid::tic_scope("meta", tid::level::detailed);
    //        // Checks if the current entries have already been written
    //        static std::unordered_map<std::string, AlgorithmStatus> save_log;
    //        bootstrap_meta_log(save_log, h5file, state_prefix);
    //        if(save_log.count(std::string(state_prefix)) and save_log.at(std::string(state_prefix)) == status) {
    //            tools::log->debug("Returning because status hasn't changed (?)");
    //            return;
    //        }
    //        if(not h5file.linkExists(state_prefix)) {
    //            tools::log->debug("Returning because state_prefix does not exist: {}", state_prefix);
    //            return;
    //        }
    //
    //        tools::log->trace("Writing attribute metadata for state_prefix: [{}]", state_prefix);
    //        auto storage_level_sv  = enum2sv(storage_level);
    //        auto storage_reason_sv = enum2sv(storage_reason);
    //        auto model_name_sv     = enum2sv(model_type);
    //        auto state_root        = fmt::format("{}/{}", status.algo_type_sv(), state_name);
    //        auto hamiltonian       = fmt::format("{}/hamiltonian", model_prefix);
    //        auto mpo_prefix        = fmt::format("{}/mpo", model_prefix);
    //        auto mps_prefix        = fmt::format("{}/mps", state_prefix);
    //        auto bonds_prefix      = fmt::format("{}/bonds", state_prefix);
    //
    //        h5pp_table_algorithm_status::register_table_type();
    //        save::attr(h5file, status, "common/status", state_prefix, "Maps state_prefix -> status", h5pp_table_algorithm_status::h5_type);
    //        save::attr(h5file, status.algorithm_has_finished, "common/finished", state_prefix, "Maps state_prefix -> finished");
    //        save::attr(h5file, storage_level_sv, "common/storage_level", state_prefix, "Maps state_prefix -> storage_level");
    //        save::attr(h5file, storage_reason_sv, "common/storage_reason", state_prefix, "Maps state_prefix -> storage_reason");
    //        save::attr(h5file, state_root, "common/state_root", state_prefix, "Maps state_prefix -> state_root");
    //        save::attr(h5file, model_prefix, "common/model_prefix", state_prefix, "Maps state_prefix -> model_prefix");
    //        save::attr(h5file, hamiltonian, "common/hamiltonian", state_prefix, "Maps state_prefix -> hamiltonian table");
    //        save::attr(h5file, mps_prefix, "common/mps_prefix", state_prefix, "Maps state_prefix -> mps_prefix");
    //        save::attr(h5file, bonds_prefix, "common/bonds_prefix", state_prefix, "Maps state_prefix -> bonds_prefix");
    //        save::attr(h5file, model_name_sv, "common/model_type", state_prefix, "Maps state_prefix -> model_type");
    //        save::attr(h5file, model_size, "common/model_size", state_prefix, "Maps state_prefix -> model_size");
    //        save::attr(h5file, status.algo_type_sv(), "common/algo_type", state_prefix, "Maps state_prefix -> algo_type");
    //        save::attr(h5file, state_name, "common/state_name", state_prefix, "Maps state_prefix -> state_name");
    //        save::attr(h5file, status.iter, "common/iter", state_prefix, "Maps state_prefix -> iter");
    //        save::attr(h5file, status.step, "common/step", state_prefix, "Maps state_prefix -> step");
    //        save::attr(h5file, status.position, "common/position", state_prefix, "Maps state_prefix -> position");
    //        if(not table_prfxs.empty()) save::attr(h5file, table_prfxs, "common/table_prfxs", state_prefix, "Maps state_prefix -> one or more table
    //        prefixes"); save_log.insert({std::string(state_prefix), status});
    //    }

    void save::timer(h5pp::File &h5file, const StorageInfo &sinfo) {
        if(not settings::timer::on) return;
        if(settings::storage::storage_level_timers == StorageLevel::NONE) return;
        if(sinfo.storage_level == StorageLevel::NONE) return;
        auto t_timers   = tid::tic_token("timers", tid::level::extra);
        auto table_path = fmt::format("{}/timers", sinfo.get_state_prefix());

        auto h5_save_point = save::get_last_save_point(h5file, table_path);
        auto save_point    = std::make_pair(sinfo.iter, sinfo.step);
        if(h5_save_point and h5_save_point.value() == save_point) return;

        // Make a table entry
        auto tid_tree    = tid::search(sinfo.algo_name);
        auto table_items = std::vector<h5pp_ur::item>();
        table_items.reserve(tid_tree.size());
        for(const auto &[i, t] : iter::enumerate(tid_tree))
            table_items.emplace_back(h5pp_ur::item{t.key, t->get_time(), t.sum, t.frac * 100, t->get_time_avg(), t->get_level(), t->get_tic_count()});

        h5pp_ur::register_table_type();
        if(not h5file.linkExists(table_path)) h5file.createTable(h5pp_ur::h5_type, table_path, fmt::format("{} Timings", sinfo.algo_name));
        h5file.writeTableRecords(table_items, table_path, 0);

        h5file.writeAttribute(sinfo.iter, table_path, "iter");
        h5file.writeAttribute(sinfo.step, table_path, "step");
        h5file.writeAttribute(sinfo.bond_lim, table_path, "bond_lim");
        h5file.writeAttribute(sinfo.bond_max, table_path, "bond_max");
    }

}
