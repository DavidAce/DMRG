
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

    void save::status(h5pp::File &h5file, const StorageInfo &sinfo, const AlgorithmStatus &status) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        auto t_hdf = tid::tic_scope("status", tid::level::highest);
        sinfo.assert_well_defined();

        // Define the table
        std::string table_path = fmt::format("{}/status", sinfo.get_state_prefix());

        // Check if the current entry has already been appended
        auto attrs = tools::common::h5::save::get_save_attrs(h5file, table_path);
        if(not attrs.link_exists) h5file.createTable(h5pp_table_algorithm_status::get_h5t(), table_path, "Algorithm Status");
        if(attrs == sinfo) return;
        auto offset = tools::common::h5::save::get_table_offset(h5file, table_path, sinfo, attrs);
        tools::log->trace("Writing to table: {} | event {} | level {} | offset {}", table_path, enum2sv(sinfo.storage_event), enum2sv(sinfo.storage_level),
                          offset);
        h5file.writeTableRecords(status, table_path, offset);
        save::set_save_attrs(h5file, table_path, sinfo);
    }

    void save::mem(h5pp::File &h5file, const StorageInfo &sinfo) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        auto t_hdf = tid::tic_scope("mem", tid::level::highest);
        sinfo.assert_well_defined();

        // Define the table
        std::string table_path = fmt::format("{}/mem_usage", sinfo.get_state_prefix());

        // Check if the current entry has already been appended
        auto attrs = tools::common::h5::save::get_save_attrs(h5file, table_path);
        if(not attrs.link_exists) h5file.createTable(h5pp_table_memory_usage::get_h5t(), table_path, "Memory usage");
        if(attrs == sinfo) return;
        auto offset = tools::common::h5::save::get_table_offset(h5file, table_path, sinfo, attrs);

        // Define the table entry
        tools::log->trace("Appending to table: {}", table_path);
        h5pp_table_memory_usage::table mem_usage_entry{};
        mem_usage_entry.iter     = sinfo.iter;
        mem_usage_entry.step     = sinfo.step;
        mem_usage_entry.position = sinfo.position;
        mem_usage_entry.event    = sinfo.storage_event;
        mem_usage_entry.bond_lim = sinfo.bond_lim;
        mem_usage_entry.rss      = debug::mem_rss_in_mb();
        mem_usage_entry.hwm      = debug::mem_hwm_in_mb();
        mem_usage_entry.vm       = debug::mem_vm_in_mb();
        tools::log->trace("Writing to table: {} | offset {}", table_path, offset);
        h5file.writeTableRecords(mem_usage_entry, table_path, offset);
        save::set_save_attrs(h5file, table_path, sinfo);
    }

    void save::timer(h5pp::File &h5file, const StorageInfo &sinfo) {
        if(settings::storage::storage_level_timers == StorageLevel::NONE) return;
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        auto t_timers = tid::tic_token("timers", tid::level::higher);

        // Define the table
        auto table_path = fmt::format("{}/timers", sinfo.get_state_prefix());
        auto attrs      = tools::common::h5::save::get_save_attrs(h5file, table_path);
        if(not attrs.link_exists) h5file.createTable(h5pp_ur::get_h5t(), table_path, fmt::format("{} Timings", sinfo.algo_name));
        if(attrs == sinfo) return;

        // Make a table entry
        auto tid_tree    = tid::search(sinfo.algo_name);
        auto table_items = std::vector<h5pp_ur::item>{};
        table_items.reserve(tid_tree.size());
        for(const auto &[i, t] : iter::enumerate(tid_tree)) {
            if(settings::storage::storage_level_timers == StorageLevel::LIGHT and t->get_level() > tid::level::normal) continue;
            if(settings::storage::storage_level_timers == StorageLevel::NORMAL and t->get_level() > tid::level::higher) continue;
            if(settings::storage::storage_level_timers == StorageLevel::FULL and t->get_level() > tid::level::highest) continue;
            table_items.emplace_back(h5pp_ur::item{t.key, t->get_time(), t.sum, t.frac * 100, t->get_time_avg(), t->get_level(), t->get_tic_count()});
        }
        tools::log->trace("Appending to table: {}", table_path);
        h5file.writeTableRecords(table_items, table_path, 0);
        save::set_save_attrs(h5file, table_path, sinfo);
    }
    std::optional<std::string> find::find_duplicate_save(const h5pp::File &h5file, std::string_view state_prefix, const AlgorithmStatus &status) {
        if(not h5file.linkExists("common/status")) return std::nullopt;
        if(not h5file.linkExists("common/state_root")) return std::nullopt;
        h5pp::Options opt_status;
        opt_status.h5Type   = h5pp_table_algorithm_status::get_h5t();
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

    StorageAttrs save::get_save_attrs(const h5pp::File &h5file, std::string_view link_path) {
        StorageAttrs attrs;
        attrs.link_exists = h5file.linkExists(link_path);
        if(attrs.link_exists) {
            auto iter = h5file.readAttribute<std::optional<uint64_t>>(link_path, "iter");
            auto step = h5file.readAttribute<std::optional<uint64_t>>(link_path, "step");
            auto blim = h5file.readAttribute<std::optional<int64_t>>(link_path, "bond_lim");
            auto bmax = h5file.readAttribute<std::optional<int64_t>>(link_path, "bond_max");
            auto trnc = h5file.readAttribute<std::optional<double>>(link_path, "trnc_lim");
            auto evnt = h5file.readAttribute<std::optional<StorageEvent>>(link_path, "storage_event");
            auto levl = h5file.readAttribute<std::optional<StorageLevel>>(link_path, "storage_level");
            if(iter) attrs.iter = iter.value();
            if(step) attrs.step = step.value();
            if(blim) attrs.bond_lim = blim.value();
            if(bmax) attrs.bond_max = bmax.value();
            if(trnc) attrs.trnc_lim = trnc.value();
            if(evnt) attrs.storage_event = evnt.value();
            if(levl) attrs.storage_level = levl.value();
        }
        return attrs;
    }
    void save::set_save_attrs(h5pp::File &h5file, std::string_view link_path, const StorageAttrs &info) {
        bool link_exists = info.link_exists or h5file.linkExists(link_path);
        if(link_exists) {
            h5file.writeAttribute(info.iter, link_path, "iter");
            h5file.writeAttribute(info.step, link_path, "step");
            h5file.writeAttribute(info.bond_lim, link_path, "bond_lim");
            h5file.writeAttribute(info.bond_max, link_path, "bond_max");
            h5file.writeAttribute(info.trnc_lim, link_path, "trnc_lim");
            h5file.writeAttribute(info.storage_event, link_path, "storage_event", std::nullopt, h5_enum_storage_event::get_h5t());
            h5file.writeAttribute(info.storage_level, link_path, "storage_level", std::nullopt, h5_enum_storage_level::get_h5t());
        }
    }
    void save::set_save_attrs(h5pp::File &h5file, std::string_view link_path, const StorageInfo &sinfo) {
        bool link_exists = h5file.linkExists(link_path);
        if(link_exists) {
            h5file.writeAttribute(sinfo.iter, link_path, "iter");
            h5file.writeAttribute(sinfo.step, link_path, "step");
            h5file.writeAttribute(sinfo.bond_lim, link_path, "bond_lim");
            h5file.writeAttribute(sinfo.bond_max, link_path, "bond_max");
            h5file.writeAttribute(sinfo.trnc_lim, link_path, "trnc_lim");
            h5file.writeAttribute(sinfo.storage_event, link_path, "storage_event", std::nullopt, h5_enum_storage_event::get_h5t());
            h5file.writeAttribute(sinfo.storage_level, link_path, "storage_level", std::nullopt, h5_enum_storage_level::get_h5t());
        }
    }
    void save::initial_state_attrs(h5pp::File &h5file, const StorageInfo &sinfo) {
        if(sinfo.storage_event != StorageEvent::INIT_STATE) return;
        if(settings::strategy::initial_pattern.empty()) {
            tools::log->warn("Could not save the initial state pattern: the pattern is empty");
            return;
        }
        auto state_prefix = sinfo.get_state_prefix();
        h5file.createGroup(state_prefix);
        // Save the initial state pattern (rather than the MPS itself) and the initial state type
        h5file.writeAttribute(enum2sv(settings::strategy::initial_state), state_prefix, "initial_state");
        h5file.writeAttribute(settings::strategy::initial_pattern, state_prefix, "initial_pattern");
        h5file.writeAttribute(enum2sv(settings::strategy::initial_type), state_prefix, "initial_type");
        h5file.writeAttribute(settings::strategy::initial_axis, state_prefix, "initial_axis");
        h5file.writeAttribute(settings::strategy::target_axis, state_prefix, "target_state_axis");
    }

    hsize_t save::get_table_offset(const h5pp::File &h5file, std::string_view table_path, const StorageInfo &sinfo, const StorageAttrs &attrs) {
        // Get the table index where we should write the current entry
        if(sinfo.storage_level == StorageLevel::LIGHT) {
            // Try to find the last occurrence of this event type, to replace it
            auto events = h5file.readTableField<std::vector<StorageEvent>>(table_path, "event", h5pp::TableSelection::ALL);
            auto offset = events.size(); // This should point one past the latest table entry, so that we append
            for(const auto &[off, evn] : iter::enumerate_reverse(events)) {
                if(evn == sinfo.storage_event) {
                    offset = off;
                    break;
                }
            }
            return offset;
        } else {
            // In this case we are saving everything, unless it is a duplicate
            if(sinfo.iter > attrs.iter) {
                // The current iter is later than any of the entries in the table, so we just append
                auto dset = h5pp::hdf5::openLink<h5pp::hid::h5d>(h5file.openFileHandle(), table_path, std::nullopt, h5file.plists.dsetAccess);
                auto dims = h5pp::hdf5::getDimensions(dset);
                return dims.front();
            } else {
                // The current iter may already have been written, so we should compare with existing
                auto events = h5file.readTableField<std::vector<StorageEvent>>(table_path, "event", h5pp::TableSelection::ALL);
                auto iters  = h5file.readTableField<std::vector<uint64_t>>(table_path, "iter", h5pp::TableSelection::ALL);
                auto offset = events.size(); // This should point one past the latest table entry, so that we append
                for(const auto &[off, evn] : iter::enumerate_reverse(events)) {
                    auto iter = iters.at(off);
                    if(evn == sinfo.storage_event and iter == sinfo.iter) {
                        // We have a match. We can overwrite this entry!
                        offset = off;
                        break;
                    }
                }
                return offset;
            }
        }
    }
}
