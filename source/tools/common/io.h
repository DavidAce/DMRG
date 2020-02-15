#pragma once
#include <string>
namespace h5pp{
    class File;
}

template <typename table_type> class class_h5table_buffer;
class class_h5table_profiling;
class class_h5table_simulation_status;
class class_simulation_status;


namespace tools::common::io {
    namespace h5dset{
        extern void write_simulation_status(const class_simulation_status &sim_status, h5pp::File &h5ppFile, const std::string &sim_name);
    }
    namespace h5table{
        extern void write_sim_status                         (const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_simulation_status> &h5tbuf);
        extern void write_profiling                          (const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_profiling> &h5tbuf);
    }

    namespace h5restore{
        extern class_simulation_status load_sim_status_from_hdf5(const h5pp::File &h5ppFile, std::string sim_name);
    }

    namespace h5tmp{
        extern std::string set_tmp_prefix(const std::string &output_filename);
        extern std::string unset_tmp_prefix(const std::string &output_filename);
        extern void copy_from_tmp(const std::string & output_filename);
        extern void create_directory(const std::string & dir);
        extern void remove_from_temp(const std::string output_filename);

    }

}
