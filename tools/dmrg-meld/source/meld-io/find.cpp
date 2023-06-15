#include "find.h"
#include "meld-io/logger.h"
#include "mpi/mpi-tools.h"
#include "tid/tid.h"
#include <regex>
namespace tools::io {
    template<bool RECURSIVE>
    std::vector<h5pp::fs::path> find_file(const h5pp::fs::path &base, const std::string &pattern) {
        auto                        t_scope = tid::tic_scope(__FUNCTION__);
        std::vector<h5pp::fs::path> result;
        std::regex                  reg(pattern);
        using dir_iterator = typename std::conditional<RECURSIVE, h5pp::fs::recursive_directory_iterator, h5pp::fs::directory_iterator>::type;
        for(auto &obj : dir_iterator(base))
            if(h5pp::fs::is_regular_file(obj) and std::regex_match(obj.path().filename().string(), reg)) result.emplace_back(obj);
        return result;
    }
    template std::vector<h5pp::fs::path> find_file<true>(const h5pp::fs::path &base, const std::string &pattern);
    template std::vector<h5pp::fs::path> find_file<false>(const h5pp::fs::path &base, const std::string &pattern);

    template<bool RECURSIVE>
    std::vector<h5pp::fs::path> find_dir(const h5pp::fs::path &base, const std::string &pattern, const std::string &subdir) {
        auto                        t_scope = tid::tic_scope(__FUNCTION__);
        std::vector<h5pp::fs::path> result;
        const std::string           regex_suffix      = ".*";
        const std::string           pattern_wo_subdir = pattern.substr(0, pattern.find(subdir)); // Peel off 'output' if it's there
        const bool                  endswith_regex = pattern_wo_subdir.find(regex_suffix, pattern_wo_subdir.size() - regex_suffix.size()) != std::string::npos;

        using dir_iterator = typename std::conditional<RECURSIVE, h5pp::fs::recursive_directory_iterator, h5pp::fs::directory_iterator>::type;
        std::regex reg(pattern_wo_subdir);
        if(not endswith_regex) reg = std::regex(pattern_wo_subdir + ".*");

        for(const auto &obj : dir_iterator(base)) {
            auto match = std::regex_match(obj.path().filename().string(), reg);
            // A valid directory has an output subdirectory
            if(h5pp::fs::is_directory(obj.path()) and match) {
                if(h5pp::fs::exists(obj.path() / subdir)) result.emplace_back(h5pp::fs::canonical(obj.path())); // Do not append subdir here
            }
        }

        return result;
    }

    template std::vector<h5pp::fs::path> find_dir<true>(const h5pp::fs::path &base, const std::string &pattern, const std::string &subdir);
    template std::vector<h5pp::fs::path> find_dir<false>(const h5pp::fs::path &base, const std::string &pattern, const std::string &subdir);

    std::vector<h5pp::fs::path> find_h5_dirs(const std::vector<h5pp::fs::path> &base, std::string_view sub, size_t max_dirs,
                                             const std::vector<std::string> &inc, const std::vector<std::string> &exc) {
        auto                        t_scope = tid::tic_scope(__FUNCTION__);
        std::vector<h5pp::fs::path> result;
        if(mpi::world.id == 0) {
            for(const auto &src_dir : base) {
                for(const auto &dir : h5pp::fs::recursive_directory_iterator(src_dir / sub, h5pp::fs::directory_options::follow_directory_symlink)) {
                    if(dir.is_directory()) {
                        // Check that the keys in inc are present in dir
                        bool include = inc.empty() or std::all_of(inc.begin(), inc.end(), [&dir](const auto &s) { return dir.path().string().find(s) != std::string::npos; });
                        bool exclude = exc.empty() or std::any_of(exc.begin(), exc.end(), [&dir](const auto &s) { return dir.path().string().find(s) != std::string::npos; });
                        if(exclude or not include) continue;
                        logger::log->info("including: {}", dir.path().string());

                        // Check if this directory has any .h5 files
                        for(const auto &obj : h5pp::fs::directory_iterator(dir, h5pp::fs::directory_options::follow_directory_symlink)) {
                            if(obj.is_regular_file() and obj.path().extension() == ".h5") {
                                result.emplace_back(dir.path());
                                break;
                            }
                        }
                        if(max_dirs > 0 and result.size() >= max_dirs) break;
                    }
                    if(max_dirs > 0 and result.size() >= max_dirs) break;
                }
            }
            std::sort(result.begin(), result.end());
            tools::logger::log->info("Found {} directories with .h5 files. Splitting into {} groups", result.size(), mpi::world.size);
        }
        mpi::scatter_r(result, 0);
        return result;
    }

}
