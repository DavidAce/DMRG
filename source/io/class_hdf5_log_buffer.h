//
// Created by david on 2018-05-24.
//

#ifndef DMRG_class_hdf5_log_BUFFER_H
#define DMRG_class_hdf5_log_BUFFER_H


#include <hdf5.h>
#include <hdf5_hl.h>
#include <memory>
#include <spdlog/spdlog.h>
#include <io/log_types.h>


class class_hdf5_file;
namespace h5pp{
    class File;
}


template<typename log_type>
class class_hdf5_log{
private:
    std::shared_ptr<h5pp::File> h5ppFile;
    std::unique_ptr<log_type> log_entries;
    std::string group_name              = "default_group";
    std::string log_name              = "default_table";
    std::string log_path;
    hsize_t recorded_elements       = 0;
    std::shared_ptr<spdlog::logger> log;
    size_t logLevel      = 3;
    std::string logName  = "DMRG";

    bool buffer_is_empty = false;
    bool log_is_ready  = false;
    bool mpi_on          = false;




    void initialize_table();
    void write_buffer_to_file();
public:
    explicit class_hdf5_log(std::string logName_ = "DMRG");
    ~class_hdf5_log();
    class_hdf5_log(std::shared_ptr<h5pp::File> h5ppFile_,
                     std::string group_name_,
                     std::string log_name_,
                     std::string logName_ = "DMRG" ,
                     bool mpi_on_ = false  );


    template<typename ...Args>
    void append_record(Args&& ...args){
        log->trace("Appending record to hdf5 table: {}", log_name);
        log_entries->buffer.emplace_back(std::forward<Args> (args)...);
        if (log_entries->buffer.size() >= log_entries->meta.chunk_size){
            write_buffer_to_file();
        }
    }
};





#endif //DMRG_class_hdf5_log_BUFFER_H
