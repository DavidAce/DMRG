#include <config/settings.h>
#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <tid/tid.h>
#include <tools/common/h5.h>
#include <tools/common/log.h>

namespace tools::common::h5 {

    void load::status(const h5pp::File &h5file, std::string_view state_prefix, AlgorithmStatus &status) {
        auto tic          = tid::tic_scope("status", tid::level::detailed);
        auto table_prefix = h5file.readAttribute<std::vector<std::string>>(state_prefix, "common/table_prfxs").front();

        std::string table_path = fmt::format("{}/status", table_prefix);
        if(h5file.linkExists(table_path)) {
            tools::log->info("Loading status from table: [{}]", table_path);
            h5file.readTableRecords(status, table_path); // Reads the last entry by default
        } else {
            throw std::runtime_error(
                fmt::format("Could not find table [status] for state [{}] in file [{}] at table path [{}]", state_prefix, h5file.getFilePath(), table_path));
        }
    }

    void load::timer(const h5pp::File &h5file, std::string_view state_prefix, [[maybe_unused]] AlgorithmStatus &status) {
        if(not settings::timer::on) return;
        auto state_root   = h5file.readAttribute<std::string>(state_prefix, "common/state_root");
        auto timer_prefix = fmt::format("{}/timers", state_root);
        if(h5file.linkExists(timer_prefix))
            tools::log->info("Loading timers from group: [{}]", timer_prefix);
        else
            throw std::runtime_error(
                fmt::format("Could not find group [timers] for state [{}] in file [{}] at prefix [{}]", state_prefix, h5file.getFilePath(), timer_prefix));
        for(const auto &dset : h5file.findDatasets("", timer_prefix)) {
            auto &t_ur = tid::get(dset);
            tools::log->trace("Loading {}", dset);
            auto data = h5file.readDataset<h5pp_ur::item>(fmt::format("{}/{}", timer_prefix, dset));
            t_ur.set_time(data.time);
            t_ur.set_count(data.count);
        }
        //        status.algo_time = tid::get_unscoped(status.algo_type_str());
    }

}
