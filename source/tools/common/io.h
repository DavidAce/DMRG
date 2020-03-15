#pragma once
#include <string>

class class_simulation_status;
namespace h5pp{
    class File;
}



namespace tools::common::io {
    namespace h5find{
        extern std::string find_table(const h5pp::File &h5ppFile, const std::string & sim_name, const std::string &table_name, const std::string & from = "last");
    }

    namespace h5dset{
        extern void write_simulation_status(const class_simulation_status &sim_status, h5pp::File &h5ppFile, const std::string &sim_name);
    }
    namespace h5table{
        extern void write_sim_status                         (const class_simulation_status &sim_status, h5pp::File & h5ppFile , const std::string & table_path);
        extern void write_profiling                          (const class_simulation_status &sim_status, h5pp::File & h5ppFile , const std::string & table_path);
    }

    namespace h5restore{
        extern class_simulation_status load_sim_status_from_hdf5(const h5pp::File &h5ppFile, const std::string & sim_name);
        extern void load_profiling_from_hdf5(const h5pp::File &h5ppFile, const std::string & sim_name);
    }

    namespace h5tmp{
        extern std::string set_tmp_prefix(const std::string &output_filename);
        extern std::string unset_tmp_prefix(const std::string &output_filename);
        extern void copy_from_tmp(const std::string & output_filename);
        extern void copy_file(const std::string & src, const std::string & tgt);
        extern void create_directory(const std::string & dir);
        extern void remove_from_temp(const std::string output_filename);

    }

}
