#pragma once
#include <map>
#include <string>
class class_simulation_status;
enum class StorageLevel;

namespace h5pp {
    class File;
}

namespace tools::common::io {
    /* clang-format off */
    namespace h5find{
        extern std::string find_latest_root(const h5pp::File & h5ppFile, const std::string & sim_name);
    }

    namespace h5table{
        extern void write_sim_status   (h5pp::File & h5ppFile, const std::string & state_prefix, const StorageLevel & storage_level, const class_simulation_status &sim_status);
        extern void write_profiling    (h5pp::File & h5ppFile, const std::string & state_prefix, const StorageLevel & storage_level, const class_simulation_status &sim_status);
        extern void write_mem_usage    (h5pp::File & h5ppFile, const std::string & state_prefix, const StorageLevel & storage_level, const class_simulation_status &sim_status);
    }

    namespace h5attr{
        extern void write_meta (h5pp::File & h5ppFile,
                                const std::string & sim_name, const std::string &state_prefix,const std::string &model_prefix, const std::string & model_type,
                                const StorageLevel & storage_level, const class_simulation_status &sim_status);
    }


    namespace h5resume{
        extern std::string             find_resumable_state     (const h5pp::File & h5ppFile, const std::string & sim_name, const std::string & search = "");
        extern void                    load_sim_status_from_hdf5(const h5pp::File & h5ppFile, const std::string & state_prefix, class_simulation_status & sim_status);
        extern void                    load_profiling_from_hdf5 (const h5pp::File & h5ppFile, const std::string & state_prefix);
    }

    namespace h5tmp{
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
        extern void copy_from_tmp(const std::string & filepath);
        extern void copy_into_tmp(const std::string & filepath);
        extern void copy_file(const std::string & src, const std::string & tgt);
        extern void create_directory(const std::string & path);
        extern void remove_from_tmp(const std::string & filepath);
    }
    /* clang-format off */
}
