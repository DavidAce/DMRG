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
        extern std::string find_latest_root(const h5pp::File & h5ppFile, const std::string & sim_name);
    }
    namespace resume{
        extern std::optional<size_t>   extract_state_number     (const std::string & state_prefix);
        extern std::string             extract_state_name       (const std::string & state_prefix);
        extern std::string             find_resumable_state     (const h5pp::File & h5ppFile, AlgorithmType algo_type, const std::string & search = "");
    }
    namespace save{
        extern void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile, const std::string &link);
        extern void bootstrap_meta_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile, const std::string &state_prefix);

        template<typename AttrType>
        extern void attr     (h5pp::File & h5ppFile, const AttrType &attrData, const std::string &attrName, const std::string &linkPath, const std::string &linkText);
        extern void status   (h5pp::File & h5ppFile, const std::string & table_prefix, const StorageLevel & storage_level, const AlgorithmStatus &status);
        extern void mem      (h5pp::File & h5ppFile, const std::string & table_prefix, const StorageLevel & storage_level, const AlgorithmStatus &status);
        extern void timer    (h5pp::File & h5ppFile, const std::string & state_prefix, const StorageLevel & storage_level, const AlgorithmStatus &status);
        extern void meta (h5pp::File &h5ppFile, const StorageLevel &storage_level, const StorageReason &storage_reason,const ModelType & model_type, size_t model_size,
                               const AlgorithmType &algo_type, const std::string & state_name, const std::string &state_prefix, const std::string &model_prefix,
                               const AlgorithmStatus &status);
    }
    namespace load {
        extern void status   (const h5pp::File & h5ppFile, const std::string & state_prefix, AlgorithmStatus & status);
        extern void timer    (const h5pp::File & h5ppFile, const std::string & state_prefix, AlgorithmStatus & status);

    }

    namespace tmp{
        namespace internal{
            struct pathpair {
             std::string original_path;
             std::string temporary_path;
            };

            inline std::map<std::string,pathpair> file_register;
            extern const pathpair & register_paths(const std::string & filepath);
            extern const pathpair & get_paths(const std::string & filepath);
            extern std::string get_tmp_dir();
        }
        extern const std::string & get_temporary_filepath(const std::string & filepath);
        extern const std::string & get_original_filepath(const std::string & filepath);
        extern void register_new_file (const std::string & filepath);
        extern void copy_from_tmp(const AlgorithmStatus &status, const h5pp::File &h5ppfile, StorageReason storage_reason, std::optional<CopyPolicy> copy_policy);
        extern void copy_from_tmp(const std::string & filepath);
        extern void copy_into_tmp(const std::string & filepath);
        extern void copy_file(const std::string & src, const std::string & tgt);
        extern void create_directory(const std::string & path);
        extern void remove_from_tmp(const std::string & filepath);
    }

    /* clang-format off */
}
