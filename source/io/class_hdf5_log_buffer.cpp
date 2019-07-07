//
// Created by david on 2018-05-24.
//
#include <iostream>
#include <iterator>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <io/nmspc_logger.h>
#include <io/class_hdf5_log_buffer.h>
#include <io/log_types.h>
#include <h5pp/h5pp.h>



template<typename log_type>
class_hdf5_log<log_type>::class_hdf5_log(std::string logName_)
:logName(logName_)
{
    log = Logger::setLogger(logName,logLevel);
}



template<typename log_type>
class_hdf5_log<log_type>::class_hdf5_log(
        std::shared_ptr<h5pp::File> h5ppFile_,
        std::string group_name_,
        std::string log_name_,
        std::string logName_,
        bool mpi_on_)
        :
        h5ppFile(std::move(h5ppFile_)),
        group_name(group_name_),
        log_name(log_name_),
        logName(logName_),
        mpi_on(mpi_on_)
{
    log = Logger::setLogger(logName,logLevel);
    log_entries      = std::make_unique<log_type>();
    log_path = group_name + "/" + log_name;
    initialize_table();
}


template<typename log_type>
class_hdf5_log<log_type>::~class_hdf5_log(){
    if (h5ppFile){
        write_buffer_to_file();
    }else if (!buffer_is_empty){
        log->warn("Output data has not been saved to file, yet it is being discarded!");
    }
}


template<typename log_type>
void class_hdf5_log<log_type>::initialize_table(){
    if (log_entries->buffer.empty() and not log_is_ready) {
        log->trace("Initializing hdf5 table: {}", log_name);
        hsize_t NRECORDS = log_entries->buffer.size();
        h5ppFile->create_group_link(group_name);
        if (not h5ppFile->linkExists(log_path)){
            hid_t file = h5ppFile->openFileHandle();
            H5TBmake_table(log_name.c_str(), file, log_path.c_str(), log_entries->meta.NFIELDS, NRECORDS,
                           log_entries->meta.dst_size, log_entries->meta.field_names.data(), log_entries->meta.dst_offsets.data(), log_entries->meta.field_types.data(),
                           log_entries->meta.chunk_size, log_entries->meta.fill_data, log_entries->meta.compress, log_entries->buffer.data());
            H5Fflush(file,H5F_SCOPE_GLOBAL);
            h5ppFile->closeFileHandle(file);
        }
        log_is_ready = true;
    }
}





template<typename log_type>
void class_hdf5_log<log_type>:: write_buffer_to_file() {
    if (!log_entries->buffer.empty() and log_is_ready) {
        log->trace("Writing buffer to hdf5 table: {}", log_name);
        hsize_t NRECORDS = log_entries->buffer.size();
        h5ppFile->create_group_link(group_name);
        hid_t file = h5ppFile->openFileHandle();
        H5TBappend_records(file, log_path.c_str(), NRECORDS, log_entries->meta.dst_size,
                           log_entries->meta.dst_offsets.data(), log_entries->meta.dst_sizes.data(), log_entries->buffer.data());

        log_entries->buffer.clear();
        recorded_elements += NRECORDS;
        H5Fflush(file,H5F_SCOPE_GLOBAL);
        h5ppFile->closeFileHandle(file);
    }
    buffer_is_empty = true;
}



//Explicit instantiations

template class class_hdf5_log<class_log_dmrg>;
template class class_hdf5_log<class_log_tebd>;
template class class_hdf5_log<class_log_profiling>;
template class class_hdf5_log<class_log_simulation_status>;
