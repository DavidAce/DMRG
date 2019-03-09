//
// Created by david on 2018-05-24.
//

#ifndef DMRG_CLASS_HDF5_TABLE_BUFFER2_H
#define DMRG_CLASS_HDF5_TABLE_BUFFER2_H


#include <hdf5.h>
#include <hdf5_hl.h>
#include <memory>
#include <algorithms/table_types.h>

class class_hdf5_file;
namespace h5pp{
    class File;
}


template<typename table_type>
class class_hdf5_table{
private:
//    std::shared_ptr<class_hdf5_file> hdf5_file;
    std::shared_ptr<h5pp::File> h5ppFile;
    std::unique_ptr<table_type> table_entries;
    std::string group_name              = "default_group";
    std::string table_name              = "default_table";
    std::string table_path;
    hsize_t recorded_elements       = 0;
    bool buffer_is_empty = false;
    bool table_is_ready  = false;
    bool mpi_on          = false;
    void initialize_table();
    void write_buffer_to_file();

public:
    class_hdf5_table() = default;
    ~class_hdf5_table();
    class_hdf5_table(std::shared_ptr<h5pp::File> h5ppFile_,
                     std::string group_name_,
                     std::string table_name_,
                     bool mpi_on_ = false  );

    template<typename ...Args>
    void append_record(Args&& ... args){
        table_entries->buffer.emplace_back(std::forward<Args> (args)...);
        if (table_entries->buffer.size() >= table_entries->meta.chunk_size){
            write_buffer_to_file();
        }
    }
};





#endif //DMRG_CLASS_HDF5_TABLE_BUFFER2_H
