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
        extern void write_all_state            (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_finite & state);
        extern void write_bond_matrices        (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_finite & state);
        extern void write_bond_matrix          (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_finite & state);
        extern void write_full_mps             (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_finite & state);
        extern void write_full_mpo             (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_finite & state);
        extern void write_model                (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_finite & state);
        extern void write_array_measurements   (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_finite & state);
    }

    namespace h5table{
        extern void write_measurements         (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_simulation_status & sim_status, const class_state_finite & state);
        extern void write_sim_status           (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_simulation_status & sim_status);
        extern void write_profiling            (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_simulation_status & sim_status);
    }

    namespace h5restore{
        extern void load_from_hdf5                      (const h5pp::File & h5ppFile, const std::string & path, class_simulation_status & sim_status,  class_state_finite & state);
        extern class_state_finite load_state_from_hdf5  (const h5pp::File & h5ppFile, const std::string & path);
    }
    /* clang-format on */
}
