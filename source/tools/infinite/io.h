#pragma once

#include <string>
class class_state_infinite;
class class_simulation_status;
enum class StorageLevel;

namespace h5pp{
    class File;
}

namespace tools::infinite::io{
    /* clang-format off */
    namespace h5dset{
        extern void write_all_state                    (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_infinite & state);
        extern void write_2site_mps                    (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_infinite & state);
        extern void write_2site_mpo                    (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_infinite & state);
        extern void write_2site_env                    (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_infinite & state);
        extern void write_2site_env2                   (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_infinite & state);
        extern void write_hamiltonian_params           (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_state_infinite & state);
        extern void write_all_measurements             (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_simulation_status & sim_status, const class_state_infinite & state);
    }

    namespace h5table{
        extern void write_measurements                 (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_simulation_status & sim_status, const class_state_infinite &state);
        extern void write_sim_status                   (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_simulation_status & sim_status);
        extern void write_profiling                    (h5pp::File & h5ppFile, const std::string & path, const StorageLevel & storage_level, const class_simulation_status & sim_status);
    }

    namespace h5restore{
        extern void load_from_hdf5                     (const h5pp::File & h5ppFile, const std::string & path, class_simulation_status & sim_status, class_state_infinite & state);
        extern void load_sim_status_from_hdf5          (const h5pp::File & h5ppFile, const std::string & path, class_simulation_status & sim_status);
        extern void load_superblock_from_hdf5          (const h5pp::File & h5ppFile, const std::string & path, class_state_infinite & state);
    }
    /* clang-format on */
}
