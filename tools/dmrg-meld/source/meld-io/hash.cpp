//
// Created by david on 2020-10-09.
//
#include "tid/tid.h"
#include <fmt/core.h>
#include <fmt/format.h>
#include <fstream>
#include <functional>
#include <h5pp/details/h5ppFilesystem.h>
#include <string>

namespace tools::hash {
    //    std::string std_hash(const std::string &str) { return std::to_string(std::hash<std::string>{}(str)); }
    size_t std_hash(const std::string &str) { return std::hash<std::string>{}(str); }

    size_t hash_file_meta(const h5pp::fs::path &fpath, const std::string &more_meta) {
        auto t_scope = tid::tic_scope(__FUNCTION__);
        auto meta =
            fmt::format("{}\n{}\n{}\n{}", fpath.string(), h5pp::fs::file_size(fpath), h5pp::fs::last_write_time(fpath).time_since_epoch().count(), more_meta);
        return std_hash(meta);
    }

}
