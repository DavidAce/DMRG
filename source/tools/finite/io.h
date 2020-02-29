#pragma once

#include <string>
class class_state_finite;
class class_simulation_status;
namespace h5pp{
    class File;
}

namespace tools::finite::io{
    namespace h5dset{
        extern void write_all_state                              (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
        extern void write_bond_matrices                          (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
        extern void write_bond_matrix                            (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
        extern void write_full_mps                               (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
        extern void write_full_mpo                               (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
        extern void write_model                                  (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
        extern void write_array_measurements                     (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);

    }


    namespace h5table{
        extern void write_measurements                       (const class_state_finite &state, const class_simulation_status &sim_status, h5pp::File & h5ppFile, const std::string & table_path);
        extern void write_sim_status                         (const class_simulation_status &sim_status, h5pp::File & h5ppFile , const std::string & table_path);
        extern void write_profiling                          (const class_simulation_status &sim_status, h5pp::File & h5ppFile , const std::string & table_path);
    }

    namespace h5restore{
        extern void load_from_hdf5                               (const h5pp::File & h5ppFile, class_state_finite & state    , class_simulation_status & sim_status, const std::string & prefix_path);
        extern class_state_finite load_state_from_hdf5           (const h5pp::File & h5ppFile, const std::string & prefix_path);
    }


}
