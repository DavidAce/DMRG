#pragma once

#include <string>

namespace tools::hash {
    std::string    sha256_file(const std::string &fn);
    std::string    md5_file(const std::string &fn);
    std::string    md5_string(const std::string &str);
    size_t         std_hash(const std::string &str);
    size_t         hash_file_meta(const h5pp::fs::path &fpath, const std::string &more_meta = "");
    constexpr long hash_length();

}
