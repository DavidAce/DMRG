//
// Created by david on 2018-01-26.
//

#ifndef DMRG_CLASS_HDF5_TABLE_BUFFER_H
#define DMRG_CLASS_HDF5_TABLE_BUFFER_H
#include <iostream>
#include <vector>
#include <memory>
#include <hdf5_hl.h>


class class_hdf5_file;

struct table_entry{
    long   chi;
    long   chi_max;
    double energy1;
    double energy2;
    double energy3;
    double energy4;
    double energy5;
    double energy6;
    double entropy;
    double variance1;
    double variance2;
    double variance3;
    double variance4;
    double variance5;
    double variance6;
    double truncation_error;
    double parity;
    int    iteration;
    int    chain_length;
    int    sweep;
    int    position;
    double time_step;
    double wall_time;
    double phys_time;

    table_entry(
    long   chi_,
    long   chi_max_,
    double energy1_,
    double energy2_,
    double energy3_,
    double energy4_,
    double energy5_,
    double energy6_,
    double entropy_,
    double variance1_,
    double variance2_,
    double variance3_,
    double variance4_,
    double variance5_,
    double variance6_,
    double truncation_error_,
    double parity_,
    int    iteration_,
    int    chain_length_,
    int    sweep_,
    int    position_,
    double time_step_,
    double wall_time_,
    double phys_time_
    ) :
            chi(chi_),
            chi_max(chi_max_),
            energy1(energy1_),
            energy2(energy2_),
            energy3(energy3_),
            energy4(energy4_),
            energy5(energy5_),
            energy6(energy6_),
            entropy(entropy_),
            variance1(variance1_),
            variance2(variance2_),
            variance3(variance3_),
            variance4(variance4_),
            variance5(variance5_),
            variance6(variance6_),
            truncation_error(truncation_error_),
            parity(parity_),
            iteration(iteration_),
            chain_length(chain_length_),
            sweep(sweep_),
            position(position_),
            time_step(time_step_),
            wall_time(wall_time_),
            phys_time(phys_time_)
            {}
};



class class_table_entry_meta{
public:
    constexpr static hsize_t NFIELDS = 24;
    size_t dst_size = sizeof( table_entry );
    size_t dst_offset[NFIELDS] = { HOFFSET(table_entry, chi),
                                   HOFFSET(table_entry, chi_max),
                                   HOFFSET(table_entry, energy1),
                                   HOFFSET(table_entry, energy2),
                                   HOFFSET(table_entry, energy3),
                                   HOFFSET(table_entry, energy4),
                                   HOFFSET(table_entry, energy5),
                                   HOFFSET(table_entry, energy6),
                                   HOFFSET(table_entry, entropy),
                                   HOFFSET(table_entry, variance1),
                                   HOFFSET(table_entry, variance2),
                                   HOFFSET(table_entry, variance3),
                                   HOFFSET(table_entry, variance4),
                                   HOFFSET(table_entry, variance5),
                                   HOFFSET(table_entry, variance6),
                                   HOFFSET(table_entry, truncation_error),
                                   HOFFSET(table_entry, parity),
                                   HOFFSET(table_entry, iteration),
                                   HOFFSET(table_entry, chain_length),
                                   HOFFSET(table_entry, sweep),
                                   HOFFSET(table_entry, position),
                                   HOFFSET(table_entry, time_step),
                                   HOFFSET(table_entry, wall_time),
                                   HOFFSET(table_entry, phys_time)};
    size_t dst_sizes[NFIELDS] = {  sizeof(table_entry::chi),
                                   sizeof(table_entry::chi_max),
                                   sizeof(table_entry::energy1),
                                   sizeof(table_entry::energy2),
                                   sizeof(table_entry::energy3),
                                   sizeof(table_entry::energy4),
                                   sizeof(table_entry::energy5),
                                   sizeof(table_entry::energy6),
                                   sizeof(table_entry::entropy),
                                   sizeof(table_entry::variance1),
                                   sizeof(table_entry::variance2),
                                   sizeof(table_entry::variance3),
                                   sizeof(table_entry::variance4),
                                   sizeof(table_entry::variance5),
                                   sizeof(table_entry::variance6),
                                   sizeof(table_entry::truncation_error),
                                   sizeof(table_entry::parity),
                                   sizeof(table_entry::iteration),
                                   sizeof(table_entry::chain_length),
                                   sizeof(table_entry::sweep),
                                   sizeof(table_entry::position),
                                   sizeof(table_entry::time_step),
                                   sizeof(table_entry::wall_time),
                                   sizeof(table_entry::phys_time)};

    const char *field_names[NFIELDS] =
            { "chi",
              "chi_max",
              "energy1",
              "energy2",
              "energy3",
              "energy4",
              "energy5",
              "energy6",
              "entropy",
              "variance1",
              "variance2",
              "variance3",
              "variance4",
              "variance5",
              "variance6",
              "truncation_error",
              "parity",
              "iteration",
              "environment_size",
              "sweep",
              "position",
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
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_DOUBLE,
                                      H5T_NATIVE_INT,
                                      H5T_NATIVE_INT,
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
        H5Tset_size( string_type, 24 );
    }
};



class class_hdf5_table_buffer : public std::vector<table_entry> {
public:
    std::shared_ptr<class_hdf5_file> hdf5_out;
    class_table_entry_meta meta;
    std::string group_name      = "default_group";
    std::string table_name      = "default_table";
    std::string table_path;
    hsize_t recorded_elements       = 0;
    bool buffer_is_empty = false;
    bool table_is_ready  = false;
    bool mpi_on          = false;
    explicit class_hdf5_table_buffer()=default;
    class_hdf5_table_buffer(std::shared_ptr<class_hdf5_file> hdf5_out_,
                            std::string group_name_,
                            std::string table_name,
                            bool mpi_on_ = false  );
    class_hdf5_table_buffer(std::nullptr_t nullp,
                            std::string group_name_,
                            std::string table_name,
                            bool mpi_on_ = false);

    class_hdf5_table_buffer(std::string group_name_,
                            std::string table_name_,
                            bool mpi_on_ = false );
    void initialize_table();
    void write_buffer_to_file();



    ~class_hdf5_table_buffer(){
        if (hdf5_out){
            if(mpi_on){
//                write_buffer_to_file_mpi();
            }
            else {
                write_buffer_to_file();
            }
        }else if (!buffer_is_empty){
            std::cerr << "Warning: Output data has not been saved to file, yet it is being discarded!\n" << std::endl;
        }
    }

    // MPI Functions
//    int mpi_rank, mpi_size;
//    void write_buffer_to_file_mpi();
//    void initialize_table_mpi();



};


#endif //DMRG_CLASS_HDF5_TABLE_BUFFER_H
