//
// Created by david on 2019-11-08.
//

#pragma once
#include <memory>
#include <spdlog/spdlog.h>
#include <io/table_types.h>


//class class_hdf5_file;
namespace h5pp{
    class File;
}


class class_h5table_buffer_dynamic{
private:
    std::shared_ptr<h5pp::File> h5ppFile;
//    std::unique_ptr<table_type> table_entries;
    std::string table_path;
    std::string group_name            = "default_group";
    std::string table_name            = "default_table";
    size_t recorded_entries       = 0;
    size_t buffer_max_size        = 10;
    bool buffer_is_empty        = false;
    bool table_is_initialized   = false;
    bool mpi_on                 = false;

    // For console logging
    std::shared_ptr<spdlog::logger> log;
    size_t logLevel      = 0;
    std::string logName  = "h5t_dynamic";


    void initialize_empty_table();
    void write_buffer_to_file();
public:
    explicit class_h5table_buffer_dynamic();
    ~class_h5table_buffer_dynamic();
    class_h5table_buffer_dynamic(std::shared_ptr<h5pp::File> h5ppFile_,
                         std::string table_path_,
                         bool mpi_on_ = false  );


//    template<typename ...Args>
//    void append_record(Args&& ...args){
//        log->trace("Appending record to output table: {}", table_name);
//        table_entries->buffer.emplace_back(std::forward<Args> (args)...);
//        buffer_is_empty = false;
//        if (get_num_buffered_entries() >= table_entries->meta.chunk_size){
//            write_buffer_to_file();
//        }
//    }
    template<typename T>
    void insert_field(const std::string &field_name);
    std::string get_table_name() const;
    std::string get_table_path() const;
    std::vector<std::string> get_field_names() const;
    size_t get_num_recorded_entries() const;
    size_t get_num_buffered_entries() const;

    void set_buffer_size(size_t new_size);

};




