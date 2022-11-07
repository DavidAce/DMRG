#pragma once
#include "config/enums.h"
#include "h5/storage_info.h"
#include <h5pp/details/h5ppHid.h>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
class AlgorithmStatus;

namespace h5pp {
    class File;
}

namespace tools::common::h5 {
    /* clang-format off */
    struct MpsInfo{
        std::string pfx;
        size_t iter = -1ul;
        size_t step = -1ul;
        StorageEvent event = StorageEvent::NONE;
    };


    namespace find{
        extern std::optional<std::string> find_duplicate_save(const h5pp::File &h5file, std::string_view state_prefix, const AlgorithmStatus & status);
    }
    namespace resume{
        extern std::optional<size_t>    extract_state_number    (std::string_view  state_prefix);
        extern std::string              extract_state_name      (std::string_view  state_prefix);
        extern std::vector<std::string> find_state_prefixes     (const h5pp::File &h5file, AlgorithmType algo_type, std::string_view name);
        extern std::vector<std::string> find_resumable_states   (const h5pp::File & h5file, AlgorithmType algo_type, std::string_view name = "", size_t iter = -1ul);
        extern std::vector<MpsInfo>     find_fully_stored_mps   (const h5pp::File & h5file, std::string_view state_prefix);
    }
    namespace save{
//        extern void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5file, std::string_view link);
//        extern void bootstrap_meta_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5file, std::string_view state_prefix);
        extern void bootstrap_meta_log(std::unordered_map<std::string, AlgorithmStatus> &save_log, const h5pp::File &h5file, std::string_view state_prefix);
        extern StorageAttrs get_save_attrs(const h5pp::File &h5file, std::string_view link_path);
        extern void         set_save_attrs(h5pp::File &h5file, std::string_view link_path, const StorageAttrs & info);
        extern void         set_save_attrs(h5pp::File &h5file, std::string_view link_path, const StorageInfo & info);
//        extern long         has_same_attrs(const h5pp::File &h5file, std::string_view link_path, const StorageInfo & info);
        extern hsize_t      get_table_offset(const h5pp::File &h5file, std::string_view table_path, const StorageInfo & sinfo, const StorageAttrs & attrs);
        template<typename AttrType>
        extern void attr     (h5pp::File & h5file, const AttrType &attrData, std::string_view linkPath, std::string_view attrName, std::string_view linkText, std::optional<h5pp::hid::h5t> h5type = std::nullopt);
        extern void status   (h5pp::File & h5file, const StorageInfo & sinfo, const AlgorithmStatus &status);
        extern void mem      (h5pp::File & h5file, const StorageInfo & sinfo);
        extern void timer    (h5pp::File & h5file, const StorageInfo & sinfo);
        extern void meta     (h5pp::File & h5file, const StorageInfo & sinfo);

    }
    namespace load {
        extern void status   (const h5pp::File & h5file, std::string_view state_prefix, AlgorithmStatus & status, const MpsInfo & info = MpsInfo());
        extern void timer    (const h5pp::File & h5file, std::string_view state_prefix, const AlgorithmStatus & status);
    }

    namespace tmp{
        namespace internal{
            struct pathpair {
             std::string original_path;
             std::string temporary_path;
            };

            inline std::map<std::string,pathpair> file_register;
            extern const pathpair & register_paths(std::string_view  filepath);
            extern const pathpair & get_paths(std::string_view  filepath);
            extern std::string get_tmp_dir();
        }
        extern std::string_view  get_temporary_filepath(std::string_view  filepath);
        extern std::string_view  get_original_filepath(std::string_view  filepath);
        extern void register_new_file (std::string_view  filepath);
        extern void copy_from_tmp(const h5pp::File &h5file, size_t iter, size_t step, StorageEvent storage_event, CopyPolicy copy_policy);
        extern void copy_from_tmp(std::string_view  filepath);
        extern void copy_into_tmp(std::string_view  filepath);
        extern void copy_file(std::string_view  src, std::string_view  tgt);
        extern void create_directory(std::string_view path);
        extern void remove_from_tmp(std::string_view filepath);
    }

    /* clang-format off */
}
