#pragma once

#include <string>
class class_state_infinite;
class class_simulation_status;

namespace h5pp{
    class File;
}



namespace tools::infinite::io{
    namespace h5dset{
        extern void write_all_state(const class_state_infinite &state, h5pp::File &h5ppFile, std::string sim_name);
        extern void write_2site_mps                    (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name);
        extern void write_2site_mpo                    (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name);
        extern void write_2site_env                    (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name);
        extern void write_2site_env2                   (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name);
        extern void write_hamiltonian_params           (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name);
        extern void write_all_measurements             (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name);
    }

    namespace h5table{
        extern void write_measurements                       (const class_state_infinite &state, const class_simulation_status &sim_status, h5pp::File & h5ppFile, const std::string & table_path);
        extern void write_sim_status                         (const class_simulation_status &sim_status, h5pp::File & h5ppFile, const std::string & table_path);
        extern void write_profiling                          (const class_simulation_status &sim_status, h5pp::File & h5ppFile, const std::string & table_path);
    }

    namespace h5restore{
        extern void load_from_hdf5                     (const h5pp::File & h5ppFile, class_state_infinite & state, class_simulation_status &sim_status, std::string sim_name);
        extern void load_superblock_from_hdf5          (const h5pp::File & h5ppFile, class_state_infinite & state, std::string sim_name);
        extern void load_sim_status_from_hdf5          (const h5pp::File & h5ppFile, class_simulation_status & sim_status, std::string sim_name);
    }
}
