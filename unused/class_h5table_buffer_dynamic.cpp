//
// Created by david on 2019-11-08.
//

#include "class_h5table_buffer_dynamic.h"
#include <iostream>
#include <iterator>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <io/nmspc_logger.h>
#include <h5pp/h5pp.h>
#include <filesystem>


namespace fs = std::filesystem;

//class dummy_table{
//    struct data_struct {
//        int dummy;
//    };
//    struct meta_struct {
//        constexpr static hsize_t NFIELDS = 1;
//        size_t dst_size = sizeof(data_struct);
//        std::array<size_t, NFIELDS>         dst_offsets     = {HOFFSET(data_struct, dummy)};
//        std::array<size_t, NFIELDS>         dst_sizes       = {sizeof(data_struct::dummy)};
//        std::array<const char *, NFIELDS>   field_names     = {"dummy"};
//        std::array<hid_t, NFIELDS>          field_types     = {H5T_NATIVE_INT};
//        hsize_t chunk_size = 4;
//        void *fill_data = nullptr;
//        int compress = 0;
//    };
//public:
//    dummy_table() = default;
//    meta_struct meta;
//    std::vector<data_struct> buffer;
//};



class_h5table_buffer_dynamic::class_h5table_buffer_dynamic()
{
    log = Logger::setLogger(logName,logLevel);
}



class_h5table_buffer_dynamic::class_h5table_buffer_dynamic(
        std::shared_ptr<h5pp::File> h5ppFile_,
        fs::path table_path_,
        bool mpi_on_)
        :
        h5ppFile(std::move(h5ppFile_)),
        table_path(table_path_.string()),
        group_name(table_path_.parent_path()),
        table_name(table_path_.filename()),
        mpi_on(mpi_on_)
{
    log = Logger::setLogger(logName,logLevel);
//    log->set_level(logLevel);
//    table_entries      = std::make_unique<table_type>();

    initialize_empty_table();
    std::cout << "GOT HERE" << std::endl;
    log->warn("Trying to initialize empty table {}", table_name);

}



class_h5table_buffer_dynamic::~class_h5table_buffer_dynamic(){
    if (h5ppFile){
        write_buffer_to_file();
    }else if (!buffer_is_empty){
        log->warn("Output data has not been saved to file, yet it is being discarded!");
    }
}



void class_h5table_buffer_dynamic::initialize_empty_table(){
    if (not table_is_initialized) {
        log->trace("Initializing output table: {}", table_name);
        hsize_t NRECORDS = 0;
        hsize_t NFIELDS = 1;
        struct data_struct {int dummy;};
        size_t dst_size = sizeof(data_struct);
        std::vector<size_t>         dst_offsets     = {HOFFSET(data_struct, dummy)};
        std::vector<size_t>         dst_sizes       = {sizeof(data_struct::dummy)};
        std::vector<const char *>   field_names     = {"dummy"};
        std::vector<hid_t>          field_types     = {H5T_NATIVE_INT};
        hsize_t chunk_size = 4;
        void *fill_data = nullptr;
        int compress = 0;

        h5ppFile->create_group_link(group_name);
        if (not h5ppFile->linkExists(table_path)){
            hid_t file = h5ppFile->openFileHandle();
            H5TBmake_table(table_name.c_str(), file, table_path.c_str(), NFIELDS, NRECORDS,
                           dst_size, field_names.data(), dst_offsets.data(), field_types.data(),
                           chunk_size, fill_data, compress, nullptr);

            h5ppFile->closeFileHandle(file);

        }
        table_is_initialized = true;

    }
}


std::vector<std::string> class_h5table_buffer_dynamic::get_field_names()const{
    if (table_is_initialized) {
        hid_t file = h5ppFile->openFileHandle();
        hsize_t    nfields;
        hsize_t    nrecords;
        H5TBget_table_info (file,table_path.c_str(), &nfields, &nrecords );

        std::vector<size_t> field_sizes(nfields);
        std::vector<size_t> field_offsets(nfields);
        std::vector<size_t> type_sizes(nfields);
        char** field_chars = new char*[nfields];
        for(hsize_t i = 0; i < nfields; ++i) field_chars[i] = new char[255];
        H5TBget_field_info(file,table_path.c_str(), field_chars, field_sizes.data(),field_offsets.data(),type_sizes.data());
        std::vector<std::string> field_names;
        for(hsize_t i = 0; i < nfields; ++i) field_names.emplace_back(field_chars[i]);
        for(hsize_t i = 0; i < nfields; ++i) delete[] field_chars[i];
        delete[] field_chars;
        h5ppFile->closeFileHandle(file);
        return field_names;
    }
    else{
        log->warn("Table has not been initialized - can't read file names");
        return std::vector<std::string>();
    }
}

template<typename T>
void class_h5table_buffer_dynamic::insert_field(const std::string &field_name){
    if (table_is_initialized) {
        auto current_field_names = get_field_names();
        auto field_iter = std::find(current_field_names.begin(), current_field_names.end(), field_name);
        if (field_iter != current_field_names.end()){
            auto pos = std::distance(current_field_names.begin(), field_iter);
            log->warn("Field is already in table: Field [{}] at position {}", field_name,pos);
            return;
        }

        hsize_t position = current_field_names.size();
        hid_t dataType = h5pp::Type::getDataType<T>();
        hid_t file = h5ppFile->openFileHandle();
        H5TBinsert_field(file,table_path.c_str(), field_name.c_str(),dataType, position, nullptr, nullptr );

        auto dummy_iter = std::find(current_field_names.begin(), current_field_names.end(), "dummy");
        if(dummy_iter != current_field_names.end()) {
            log->trace("Deleting dummy in table: {}", table_name);
            H5TBdelete_field(file, table_path.c_str(),dummy_iter->c_str());
        }
        h5ppFile->closeFileHandle(file);
    }
}

template void class_h5table_buffer_dynamic::insert_field<double>(const std::string &field_name);
template void class_h5table_buffer_dynamic::insert_field<int>(const std::string &field_name);





void class_h5table_buffer_dynamic::write_buffer_to_file() {
//    if (!table_entries->buffer.empty() and table_is_initialized) {
//        log->trace("Writing buffer to output table: {}", table_name);
//        hsize_t NRECORDS = get_num_buffered_entries();
//        h5ppFile->create_group_link(group_name);
//        hid_t file = h5ppFile->openFileHandle();
//        H5TBappend_records(file, table_path.c_str(), NRECORDS, table_entries->meta.dst_size,
//                           table_entries->meta.dst_offsets.data(), table_entries->meta.dst_sizes.data(), table_entries->buffer.data());
//
//        table_entries->buffer.clear();
//        recorded_entries += NRECORDS;
//        H5Fflush(file,H5F_SCOPE_GLOBAL);
//        h5ppFile->closeFileHandle(file);
//    }
//    buffer_is_empty = true;
}

std::string class_h5table_buffer_dynamic::get_table_name()           const{return table_name;}
std::string class_h5table_buffer_dynamic::get_table_path()           const{return table_path;}
size_t      class_h5table_buffer_dynamic::get_num_recorded_entries() const{return recorded_entries;}
size_t      class_h5table_buffer_dynamic::get_num_buffered_entries() const{return 0 ;/*return table_entries->buffer.size();*/}
void        class_h5table_buffer_dynamic::set_buffer_size(size_t new_size) { buffer_max_size = new_size;}
