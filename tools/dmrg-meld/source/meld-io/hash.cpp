//
// Created by david on 2020-10-09.
//
#include "tid/tid.h"
#include <fstream>
#include <functional>
#include <h5pp/details/h5ppFilesystem.h>
#include <string>
namespace tools::hash {
    std::string std_hash(const std::string &str) { return std::to_string(std::hash<std::string>{}(str)); }

    std::string hash_file_meta(const h5pp::fs::path &fpath, const std::string &more_meta) {
        auto        t_scope = tid::tic_scope(__FUNCTION__);
        std::string meta;
        meta.reserve(512);
        meta += fpath.string() + '\n';
        meta += std::to_string(h5pp::fs::last_write_time(fpath).time_since_epoch().count()) + '\n';
        if(not more_meta.empty()) meta += more_meta + '\n';
        return std_hash(meta);
        //        return md5_string(meta);
    }

}
