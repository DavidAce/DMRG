#include "config/settings.h"
#include "debug/exceptions.h"
#include "io/hdf5_types.h"
#include "tid/tid.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include <h5pp/h5pp.h>

namespace tools::common::h5 {

    void load::status(const h5pp::File &h5file, std::string_view state_prefix, AlgorithmStatus &status) {
        auto tic          = tid::tic_scope("status", tid::level::detailed);
        auto table_prefix = h5file.readAttribute<std::vector<std::string>>("common/table_prfxs", state_prefix).front();

        std::string table_path = fmt::format("{}/status", table_prefix);
        if(h5file.linkExists(table_path)) {
            tools::log->info("Loading status from table: [{}]", table_path);
            h5file.readTableRecords(status, table_path); // Reads the last entry by default
        } else {
            throw except::runtime_error("Could not find table [status] for state [{}] in file [{}] at table path [{}]", state_prefix, h5file.getFilePath(),
                                        table_path);
        }
    }

    void load::timer(const h5pp::File &h5file, std::string_view state_prefix, [[maybe_unused]] const AlgorithmStatus &status) {
        if(not settings::timer::on) return;
        auto state_root = h5file.readAttribute<std::string>("common/state_root", state_prefix);
        auto table_path = fmt::format("{}/tables/timers", state_root);
        if(h5file.linkExists(table_path))
            tools::log->info("Loading timers from table: [{}]", table_path);
        else
            throw except::runtime_error("Could not find table [{}] for state [{}] in file [{}] ", table_path, state_prefix, h5file.getFilePath());

        auto table = h5file.readTableRecords<std::vector<h5pp_ur::item>>(table_path, h5pp::TableSelection::ALL);
        for(const auto &t : table) {
            tools::log->trace("Loading {}", t.name);
            auto &t_ur = tid::get(t.name);
            t_ur.set_time(t.time);
            t_ur.set_level(static_cast<tid::level>(t.level));
            t_ur.set_count(t.count);
        }
    }

}
