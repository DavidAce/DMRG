#include <config/settings.h>
#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <tid/tid.h>
#include <tools/common/h5.h>
#include <tools/common/log.h>

namespace tools::common::h5 {

    void load::status(const h5pp::File &h5file, std::string_view state_prefix, AlgorithmStatus &status) {
        auto        tic        = tid::tic_scope("status", tid::level::pedant);
        std::string table_path = fmt::format("{}/status", state_prefix);
        if(h5file.linkExists(table_path)) {
            tools::log->info("Loading status from table: [{}]", table_path);
            h5file.readTableRecords(status, table_path); // Reads the last entry by default
        } else {
            throw std::runtime_error(
                fmt::format("Could not find table [status] in file [{}] at prefix [{}] at path [{}]", h5file.getFilePath(), state_prefix, table_path));
        }
    }

    void load::timer(const h5pp::File &h5file, std::string_view state_prefix, [[maybe_unused]] AlgorithmStatus &status) {
        if(not settings::profiling::on) return;
        std::string table_path = fmt::format("{}/profiling", state_prefix);
        if(h5file.linkExists(table_path))
            tools::log->info("Loading profiling from table: [{}]", table_path);
        else
            throw std::runtime_error(
                fmt::format("Could not find table [profiling] in file [{}] at prefix [{}] at path [{}]", h5file.getFilePath(), state_prefix, table_path));
        for(const auto &dset : h5file.findDatasets("", table_path)) {
            auto &t_ur = tid::get(dset);
            tools::log->trace("Loading {}", dset);
            auto data = h5file.readDataset<h5pp_ur::item>(fmt::format("{}/{}", table_path, dset));
            t_ur.set_time(data.time);
            t_ur.set_count(data.count);
        }
        //        status.algo_time = tid::get_unscoped(status.algo_type_str());
    }

}
