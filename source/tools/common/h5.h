#pragma once
#include <map>
#include <string>
class AlgorithmStatus;
enum class StorageLevel;
enum class ModelType;
enum class AlgorithmType;
enum class StorageReason;
enum class CopyPolicy;

namespace h5pp {
    class File;
}

namespace tools::common::h5 {
    /* clang-format off */
    namespace h5find{
        extern std::string find_latest_root(const h5pp::File & h5ppFile, std::string_view  sim_name);
    }
    namespace resume{
        extern std::optional<size_t>   extract_state_number     (std::string_view  state_prefix);
        extern std::string             extract_state_name       (std::string_view  state_prefix);
        extern std::string             find_resumable_state     (const h5pp::File & h5ppFile, AlgorithmType algo_type, std::string_view  search = "");
    }
    namespace save{
        extern void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile, std::string_view link);
        extern void bootstrap_meta_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile, std::string_view state_prefix);

        template<typename AttrType>
        extern void attr     (h5pp::File & h5ppFile, const AttrType &attrData, std::string_view attrName, std::string_view linkPath, std::string_view linkText);
        extern void status   (h5pp::File & h5ppFile, std::string_view  table_prefix, const StorageLevel & storage_level, const AlgorithmStatus &status);
        extern void mem      (h5pp::File & h5ppFile, std::string_view  table_prefix, const StorageLevel & storage_level, const AlgorithmStatus &status);
        extern void timer    (h5pp::File & h5ppFile, std::string_view  state_prefix, const StorageLevel & storage_level, const AlgorithmStatus &status);
        extern void meta (h5pp::File &h5ppFile, const StorageLevel &storage_level, const StorageReason &storage_reason,const ModelType & model_type,
                          size_t model_size, std::string_view  state_name, std::string_view state_prefix, std::string_view model_prefix,
                               const AlgorithmStatus &status);
    }
    namespace load {
        extern void status   (const h5pp::File & h5ppFile, std::string_view  state_prefix, AlgorithmStatus & status);
        extern void timer    (const h5pp::File & h5ppFile, std::string_view  state_prefix, AlgorithmStatus & status);

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
        extern void copy_from_tmp(const AlgorithmStatus &status, const h5pp::File &h5ppfile, StorageReason storage_reason, std::optional<CopyPolicy> copy_policy);
        extern void copy_from_tmp(std::string_view  filepath);
        extern void copy_into_tmp(std::string_view  filepath);
        extern void copy_file(std::string_view  src, std::string_view  tgt);
        extern void create_directory(std::string_view path);
        extern void remove_from_tmp(std::string_view filepath);
    }

    /* clang-format off */
}