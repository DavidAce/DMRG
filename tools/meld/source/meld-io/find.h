#pragma once

#include <h5pp/details/h5ppFilesystem.h>
#include <h5pp/details/h5ppFormat.h>
#include <vector>

namespace tools::io {
    template<bool RECURSIVE>
    std::vector<h5pp::fs::path> find_file(const h5pp::fs::path &base, const std::string &pattern);

    template<bool RECURSIVE>
    std::vector<h5pp::fs::path> find_dir(const h5pp::fs::path &base, const std::string &pattern, const std::string &subdir = "/output");

    std::vector<h5pp::fs::path> find_h5_dirs(const std::vector<h5pp::fs::path> &base, std::string_view sub, size_t max_dirs,
                                             const std::vector<std::string> &inc, const std::vector<std::string> &exc);

}
