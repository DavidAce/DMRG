#pragma once

#include <string>
class class_tensors_infinite;
class class_state_infinite;
class class_model_infinite;
class class_algorithm_status;
enum class StorageLevel;
enum class SimulationTask;
namespace h5pp {
    class File;
}

namespace tools::infinite::io {
    /* clang-format off */
    namespace h5dset{
        extern int decide_layout(std::string_view prefix_path);
        extern void save_state            (h5pp::File & h5ppFile, const std::string & state_prefix, const StorageLevel & storage_level, const class_state_infinite & state, const class_algorithm_status &status);
        extern void save_edges            (h5pp::File & h5ppFile, const std::string & edges_prefix, const StorageLevel & storage_level, const class_edges_infinite & edges);
        extern void save_model            (h5pp::File & h5ppFile, const std::string & mpo_path, const StorageLevel & storage_level, const class_model_infinite & model);
    }

    namespace h5table{
        extern void save_measurements   (h5pp::File & h5ppFile, const std::string & table_path, const StorageLevel & storage_level, const class_tensors_infinite & tensors, const class_algorithm_status & status);
        extern void save_model          (h5pp::File & h5ppFile, const std::string & table_path, const StorageLevel & storage_level, const class_model_infinite & model);
        extern void save_sim_status     (h5pp::File & h5ppFile, const std::string & table_path, const StorageLevel & storage_level, const class_algorithm_status & status);
        extern void save_profiling      (h5pp::File & h5ppFile, const std::string & table_path, const StorageLevel & storage_level, const class_algorithm_status & status);
        extern void save_mem_usage      (h5pp::File & h5ppFile, const std::string & table_path, const StorageLevel & storage_level, const class_algorithm_status & status);
    }

    namespace h5resume{
        extern void load_tensors (const h5pp::File & h5ppFile, const std::string & state_prefix, class_tensors_infinite & state, class_algorithm_status & status);
        extern void load_state (const h5pp::File & h5ppFile, const std::string & state_prefix, class_state_infinite & state, const class_algorithm_status & status);
        extern void load_model (const h5pp::File & h5ppFile, const std::string & state_prefix, class_model_infinite & state, const class_algorithm_status & status);
        extern void validate (const h5pp::File & h5ppFile, const std::string & state_prefix, const class_tensors_infinite & state, const class_algorithm_status & status);
        extern std::vector<SimulationTask>
            getTaskList(const h5pp::File &h5ppFile, const std::string &sim_name, const std::string &state_prefix, const class_state_infinite & state, const class_algorithm_status & status);

    }
    /* clang-format on */
}
