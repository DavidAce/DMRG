//
// Created by david on 2018-05-24.
//
#include <iostream>
#include <iterator>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <io/nmspc_logger.h>
#include <io/class_h5table_buffer.h>
#include <io/table_types.h>
#include <h5pp/h5pp.h>


template<typename table_type>
class_h5table_buffer<table_type>::class_h5table_buffer()
{
    log = Logger::setLogger(logName,logLevel);
}



template<typename table_type>
class_h5table_buffer<table_type>::class_h5table_buffer(
        std::shared_ptr<h5pp::File> h5ppFile_,
        std::string table_path_,
        bool mpi_on_)
        :
        h5ppFile(std::move(h5ppFile_)),
        table_path(table_path_),
        mpi_on(mpi_on_)
{
    group_name = h5pp::fs::path(table_path).parent_path();
    table_name = h5pp::fs::path(table_path).filename();
    log = Logger::setLogger(logName,logLevel);
    table_entries      = std::make_unique<table_type>();
    initialize_table();
    set_buffer_size(table_entries->meta.chunk_size);
}


template<typename table_type>
class_h5table_buffer<table_type>::~class_h5table_buffer(){
    if (h5ppFile){
        write_buffer_to_file();
    }else if (!buffer_is_empty){
        log->warn("Output data has not been saved to file, yet it is being discarded!");
    }
}


template<typename table_type>
void class_h5table_buffer<table_type>::initialize_table(){
    if (table_entries->buffer.empty() and not table_is_initialized) {
        log->trace("Initializing output table: {}", table_name);
        hsize_t NRECORDS = get_num_buffered_entries();
        h5ppFile->createGroup(group_name);
        if (not h5ppFile->linkExists(table_path)){
            auto file = h5ppFile->openFileHandle();
            H5TBmake_table(table_name.c_str(), file, table_path.c_str(), table_entries->meta.NFIELDS, NRECORDS,
                           table_entries->meta.dst_size, table_entries->meta.field_names.data(), table_entries->meta.dst_offsets.data(), table_entries->meta.field_types.data(),
                           table_entries->meta.chunk_size, table_entries->meta.fill_data, table_entries->meta.compress, table_entries->buffer.data());
//            H5Fflush(file,H5F_SCOPE_GLOBAL);
        }
        table_is_initialized = true;
    }
}





template<typename table_type>
void class_h5table_buffer<table_type>:: write_buffer_to_file() {
    if (!table_entries->buffer.empty() and table_is_initialized) {
        log->trace("Writing buffer to output table: {}", table_name);
        hsize_t NRECORDS = get_num_buffered_entries();
        h5ppFile->createGroup(group_name);
        auto file = h5ppFile->openFileHandle();
        H5TBappend_records(file, table_path.c_str(), NRECORDS, table_entries->meta.dst_size,
                           table_entries->meta.dst_offsets.data(), table_entries->meta.dst_sizes.data(), table_entries->buffer.data());

        table_entries->buffer.clear();
        recorded_entries += NRECORDS;
//        H5Fflush(file,H5F_SCOPE_GLOBAL);
    }
    buffer_is_empty = true;
}


template<typename table_type> std::string class_h5table_buffer<table_type>::get_table_name()           const{return table_name;}
template<typename table_type> std::string class_h5table_buffer<table_type>::get_table_path()           const{return table_path;}
template<typename table_type> size_t      class_h5table_buffer<table_type>::get_num_recorded_entries() const{return recorded_entries;}
template<typename table_type> size_t      class_h5table_buffer<table_type>::get_num_buffered_entries() const{return table_entries->buffer.size();}
template<typename table_type> void        class_h5table_buffer<table_type>::set_buffer_size(size_t new_size) { buffer_max_size = new_size;}


//Explicit instantiations

template class class_h5table_buffer<class_h5table_measurements_finite>;
template class class_h5table_buffer<class_h5table_measurements_infinite>;
template class class_h5table_buffer<class_h5table_profiling>;
template class class_h5table_buffer<class_h5table_simulation_status>;
