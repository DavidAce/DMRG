//
// Created by david on 2018-01-26.
//

#ifndef DMRG_CLASS_HDF5_TABLE_BUFFER_H
#define DMRG_CLASS_HDF5_TABLE_BUFFER_H
#include <iostream>
#include <vector>
#include <memory>
#include <H5Cpp.h>


class class_hdf5_file;

struct table_entry{
    long   chi;
    long   chi_max;
    double energy;
    double entropy;
    double variance1;
    double variance2;
    double truncation_error;
    int    iteration;
    int    chain_length;
    double time_step;
    double wall_time;
    double phys_time;

    table_entry(
    long   chi,
    long   chi_max,
    double energy,
    double entropy,
    double variance1,
    double variance2,
    double truncation_error,
    int    iteration,
    int    chain_length,
    double time_step,
    double wall_time,
    double phys_time
    ) :
            chi(chi),
            chi_max(chi_max),
            energy(energy),
            entropy(entropy),
            variance1(variance1),
            variance2(variance2),
            truncation_error(truncation_error),
            iteration(iteration),
            chain_length(chain_length),
            time_step(time_step),
            wall_time(wall_time),
            phys_time(phys_time)
            {}
};



class class_table_entry_meta{
public:
    constexpr static hsize_t NFIELDS = 12;
    size_t dst_size = sizeof( table_entry );
    size_t dst_offset[NFIELDS] = { HOFFSET(table_entry, chi),
                                   HOFFSET(table_entry, chi_max),
                                   HOFFSET(table_entry, energy),
                                   HOFFSET(table_entry, entropy),
                                   HOFFSET(table_entry, variance1),
                                   HOFFSET(table_entry, variance2),
                                   HOFFSET(table_entry, truncation_error),
                                   HOFFSET(table_entry, iteration),
                                   HOFFSET(table_entry, chain_length),
                                   HOFFSET(table_entry, time_step),
                                   HOFFSET(table_entry, wall_time),
                                   HOFFSET(table_entry, phys_time)};
    size_t dst_sizes[NFIELDS] = {  sizeof(table_entry::chi),
                                   sizeof(table_entry::chi_max),
                                   sizeof(table_entry::energy),
                                   sizeof(table_entry::entropy),
                                   sizeof(table_entry::variance1),
                                   sizeof(table_entry::variance2),
                                   sizeof(table_entry::truncation_error),
                                   sizeof(table_entry::iteration),
                                   sizeof(table_entry::chain_length),
                                   sizeof(table_entry::time_step),
                                   sizeof(table_entry::wall_time),
                                   sizeof(table_entry::phys_time)};

    const char *field_names[NFIELDS] =
            { "chi",
              "chi_max",
              "energy",
              "entropy",
              "variance1",
              "variance2",
              "truncation_error",
              "position",
              "chain_length",
              "time_step",
              "wall_time",
              "phys_time"};
    /* Define field information */
    hid_t      field_type[NFIELDS] = {H5T_NATIVE_LONG,
                                      H5T_NATIVE_LONG,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_INT,
                                      H5T_NATIVE_INT,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE};
    hid_t      string_type = H5Tcopy( H5T_C_S1 );
    hsize_t    chunk_size  = NFIELDS;
    int        *fill_data  = NULL;
    int        compress    = 0;
    int        i;
    class_table_entry_meta(){
        H5Tset_size( string_type, 16 );
    }
};



class class_hdf5_table_buffer : public std::vector<table_entry> {
public:
    std::shared_ptr<class_hdf5_file> hdf5_out;
    class_table_entry_meta meta;
    std::string group_name      = "default_group";
    std::string table_name      = "default_table";

    int max_elements            = 5000;
    std::string table_relative_name;
    bool buffer_is_empty = false;
    bool table_is_empty  = true;
    explicit class_hdf5_table_buffer()=default;
    class_hdf5_table_buffer(std::shared_ptr<class_hdf5_file> hdf5_out_,
                            std::string group_name_,
                            std::string table_name);
    class_hdf5_table_buffer(std::nullptr_t nullp,
                            std::string group_name_,
                            std::string table_name);

    class_hdf5_table_buffer(std::string group_name_,
                            std::string table_name_);

    ~class_hdf5_table_buffer(){
        if (hdf5_out){
            write_buffer_to_file();
        }else if (!buffer_is_empty){
            std::cerr << "Warning: Output data has not been saved to file, yet it is being discarded!\n" << std::endl;
        }
    }
    void write_buffer_to_file();
};


#endif //DMRG_CLASS_HDF5_TABLE_BUFFER_H
