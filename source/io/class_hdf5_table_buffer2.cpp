//
// Created by david on 2018-05-24.
//
#include <iostream>
#include <iterator>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <io/nmspc_logger.h>
#include <io/class_hdf5_table_buffer2.h>
#include <algorithms/table_types.h>
#include <h5pp/h5pp.h>


/*! \brief Prints the content of a vector nicely */
template<typename T, auto N>
std::ostream &operator<<(std::ostream &out, const std::array<T,N> &v) {
    if (!v.empty()) {
        out << "[ ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
        out << "]";
    }
    return out;
}

template<typename table_type>
class_hdf5_table<table_type>::class_hdf5_table(std::string logName_)
:logName(logName_)
{
    log = Logger::setLogger(logName,logLevel);
}



template<typename table_type>
class_hdf5_table<table_type>::class_hdf5_table(
        std::shared_ptr<h5pp::File> h5ppFile_,
        std::string group_name_,
        std::string table_name_,
        std::string logName_,
        bool mpi_on_)
        :
        h5ppFile(std::move(h5ppFile_)),
        group_name(group_name_),
        table_name(table_name_),
        logName(logName_),
        mpi_on(mpi_on_)
{
    log = Logger::setLogger(logName,logLevel);
    table_entries      = std::make_unique<table_type>();
    table_path = group_name + "/" + table_name;
    initialize_table();
}


template<typename table_type>
class_hdf5_table<table_type>::~class_hdf5_table(){
    if (h5ppFile){
        write_buffer_to_file();
    }else if (!buffer_is_empty){
        log->warn("Output data has not been saved to file, yet it is being discarded!");
    }
}


template<typename table_type>
void class_hdf5_table<table_type>::initialize_table(){
    if (table_entries->buffer.empty() and !table_is_ready) {
        log->trace("Initializing hdf5 table: {}", table_name);
        hsize_t NRECORDS = table_entries->buffer.size();
        h5ppFile->create_group_link(group_name);
        if (not h5ppFile->linkExists(table_path)){
            hid_t file = h5ppFile->openFileHandle();
            H5TBmake_table(table_name.c_str(), file, table_path.c_str(), table_entries->meta.NFIELDS, NRECORDS,
                           table_entries->meta.dst_size, table_entries->meta.field_names.data(), table_entries->meta.dst_offsets.data(), table_entries->meta.field_types.data(),
                           table_entries->meta.chunk_size, table_entries->meta.fill_data, table_entries->meta.compress, table_entries->buffer.data());
            H5Fflush(file,H5F_SCOPE_GLOBAL);
            h5ppFile->closeFileHandle(file);
        }
        table_is_ready = true;
    }
}

template<typename table_type>
void class_hdf5_table<table_type>:: write_buffer_to_file() {
    if (!table_entries->buffer.empty() and table_is_ready) {
        log->trace("Writing buffer to hdf5 table: {}", table_name);
        hsize_t NRECORDS = table_entries->buffer.size();
        h5ppFile->create_group_link(group_name);
        hid_t file = h5ppFile->openFileHandle();
        H5TBappend_records(file, table_path.c_str(), NRECORDS, table_entries->meta.dst_size,
                           table_entries->meta.dst_offsets.data(), table_entries->meta.dst_sizes.data(), table_entries->buffer.data());

        table_entries->buffer.clear();
        recorded_elements += NRECORDS;
        H5Fflush(file,H5F_SCOPE_GLOBAL);
        h5ppFile->closeFileHandle(file);
    }
    buffer_is_empty = true;
}



//Explicit instantiations

template class class_hdf5_table<class_table_dmrg>;
template class class_hdf5_table<class_table_tebd>;
template class class_hdf5_table<class_table_finite_chain>;
template class class_hdf5_table<class_table_profiling>;
