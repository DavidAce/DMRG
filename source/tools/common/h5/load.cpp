#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "io/hdf5_types.h"
#include "tid/tid.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include <h5pp/h5pp.h>

namespace tools::common::h5 {
    void load::status(const h5pp::File &h5file, std::string_view state_prefix, AlgorithmStatus &status, const MpsInfo &info) {
        auto tic        = tid::tic_scope("status", tid::level::highest);
        auto table_path = fmt::format("{}/status", state_prefix);
        if(h5file.linkExists(table_path)) {
            tools::log->info("Loading status from table: [{}]", table_path);
            if(info.iter != -1ul or info.step != -1ul or info.event != StorageEvent::NONE) {
                // Load the whole table and find the right entry
                auto statusRecords = h5file.readTableRecords<std::vector<AlgorithmStatus>>(table_path, h5pp::TableSelection::ALL);
                for(const auto &statusRecord : iter::reverse(statusRecords)) {
                    if(info.iter != 1ul and statusRecord.iter != info.iter) continue;
                    if(info.step != 1ul and statusRecord.step != info.step) continue;
                    if(info.event != StorageEvent::NONE and statusRecord.event != info.event) continue;
                    // We have a match!
                    status = statusRecord;
                    return;
                }
                throw except::load_error("Could not find status entry matching iter {} | step {} | event {}", info.iter, info.step, enum2sv(info.event));
            } else {
                h5file.readTableRecords(status, table_path, h5pp::TableSelection::LAST); // Reads the last entry by default
            }
        } else {
            throw except::runtime_error("Could not find table [status] for state [{}] in file [{}] at table path [{}]", state_prefix, h5file.getFilePath(),
                                        table_path);
        }
    }

    void load::timer(const h5pp::File &h5file, std::string_view state_prefix, [[maybe_unused]] const AlgorithmStatus &status) {
        if(not settings::timer::on) return;
        auto table_path = fmt::format("{}/timers", state_prefix);
        if(h5file.linkExists(table_path)) {
            tools::log->info("Loading timers from table: [{}]", table_path);
            auto table = h5file.readTableRecords<std::vector<h5pp_ur::item>>(table_path, h5pp::TableSelection::ALL);
            for(const auto &t : table) {
                tools::log->trace("Loading {}", t.name);
                auto &t_ur = tid::get(t.name);
                t_ur.set_time(t.time);
                t_ur.set_level(static_cast<tid::level>(t.level));
                t_ur.set_count(t.count);
            }
        } else {
            tools::log->info("Could not load table [{}]", table_path);
        }
    }
}
