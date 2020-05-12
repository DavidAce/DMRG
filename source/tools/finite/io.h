#pragma once
#include <string>

class class_state_finite;
class class_simulation_status;
enum class StorageLevel;
namespace h5pp {
    class File;
}

namespace tools::finite::io {
    /* clang-format off */
    namespace h5dset{
        extern int decide_layout(std::string_view prefix_path);
        extern void write_mps            (h5pp::File & h5ppFile, const std::string & state_prefix, const StorageLevel & storage_level, const class_state_finite & state);
        extern void write_mpo            (h5pp::File & h5ppFile, const std::string & model_prefix, const StorageLevel & storage_level, const class_state_finite & state);
        extern void write_ent            (h5pp::File & h5ppFile, const std::string & state_prefix, const StorageLevel & storage_level, const class_state_finite & state);
    }

    namespace h5table{
        extern void write_model          (h5pp::File & h5ppFile, const std::string & model_prefix, const StorageLevel & storage_level, const class_state_finite & state);
        extern void write_measurements   (h5pp::File & h5ppFile, const std::string & table_prefix, const StorageLevel & storage_level, const class_state_finite & state, const class_simulation_status & sim_status);
        extern void write_sim_status     (h5pp::File & h5ppFile, const std::string & table_prefix, const StorageLevel & storage_level, const class_simulation_status & sim_status);
        extern void write_profiling      (h5pp::File & h5ppFile, const std::string & table_prefix, const StorageLevel & storage_level, const class_simulation_status & sim_status);
        extern void write_mem_usage      (h5pp::File & h5ppFile, const std::string & table_prefix, const StorageLevel & storage_level, const class_simulation_status & sim_status);

    }

    namespace h5resume{
        extern void load_all (const h5pp::File & h5ppFile, const std::string & state_prefix, class_state_finite & state, class_simulation_status & sim_status);
        extern void load_mps (const h5pp::File & h5ppFile, const std::string & state_prefix, class_state_finite & state, const class_simulation_status & sim_status);
        extern void load_mpo (const h5pp::File & h5ppFile, const std::string & state_prefix, class_state_finite & state, const class_simulation_status & sim_status);
        extern void validate (const h5pp::File & h5ppFile, const std::string & state_prefix, const class_state_finite & state, const class_simulation_status & sim_status);
        extern std::list<SimulationTask>
            getTaskList(const h5pp::File &h5ppFile, const std::string &sim_name, const std::string &state_prefix, const class_state_finite & state, const class_simulation_status & sim_status);

    }

    /* clang-format on */
}
