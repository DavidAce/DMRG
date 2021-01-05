#pragma once
#include <string>
class class_state_finite;
class class_model_finite;
class class_edges_finite;
class class_tensors_finite;
class class_algorithm_status;
enum class StorageLevel;
enum class AlgorithmType;
namespace h5pp {
    class File;
}

namespace tools::finite::io {
    /* clang-format off */
    namespace h5dset{
        extern int decide_layout(std::string_view prefix_path);
        extern void save_state            (h5pp::File & h5ppFile, const std::string & state_prefix, const StorageLevel & storage_level, const class_state_finite & state, const class_algorithm_status & status);
        extern void save_model            (h5pp::File & h5ppFile, const std::string & model_prefix, const StorageLevel & storage_level, const class_model_finite & model);
        extern void save_entgm            (h5pp::File & h5ppFile, const std::string & state_prefix, const StorageLevel & storage_level, const class_state_finite & state, const class_algorithm_status & status);

    }

    namespace h5table{
        extern void save_measurements   (h5pp::File & h5ppFile, const std::string & table_path, const StorageLevel & storage_level, const class_tensors_finite & tensors, const class_algorithm_status & status, AlgorithmType algo_type);
        extern void save_model          (h5pp::File & h5ppFile, const std::string & table_path, const StorageLevel & storage_level, const class_model_finite & model);
        extern void save_sim_status     (h5pp::File & h5ppFile, const std::string & table_path, const StorageLevel & storage_level, const class_algorithm_status & status);
        extern void save_profiling      (h5pp::File & h5ppFile, const std::string & table_path, const StorageLevel & storage_level, const class_algorithm_status & status);
        extern void save_mem_usage      (h5pp::File & h5ppFile, const std::string & table_path, const StorageLevel & storage_level, const class_algorithm_status & status);
    }

    namespace h5resume{
        extern void load_simulation (const h5pp::File & h5ppFile, const std::string & state_prefix, class_tensors_finite & tensors, class_algorithm_status & status);
        extern void load_state   (const h5pp::File & h5ppFile, const std::string & state_prefix, class_state_finite & state, const class_algorithm_status & status);
        extern void load_model   (const h5pp::File & h5ppFile, const std::string & state_prefix, class_model_finite & model);
        extern void validate (const h5pp::File & h5ppFile, const std::string & state_prefix, class_tensors_finite & tensors);
//        extern std::list<SimulationTask>
//            getTaskList(const h5pp::File &h5ppFile, const std::string &sim_name, const std::string &state_prefix, const class_tensors_finite & tensors, const class_algorithm_status & status);
    }
    /* clang-format on */
}
