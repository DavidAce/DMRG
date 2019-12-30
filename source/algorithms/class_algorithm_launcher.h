//
// Created by david on 7/30/17.
//
#pragma once

#include <memory>
#include <spdlog/spdlog.h>
namespace h5pp{class File;}

class class_algorithm_launcher  {
private:
    std::shared_ptr<spdlog::logger> log;
    void setLogger(std::string name);
public:

//    std::shared_ptr <class_hdf5_file> output;
    std::shared_ptr<h5pp::File> h5ppFile;

    class_algorithm_launcher(std::shared_ptr<h5pp::File> h5ppFile_);
    class_algorithm_launcher();

    void run_algorithms();
    void run_iDMRG();
    void run_fDMRG();
    void run_xDMRG();
    void run_iTEBD();

//    static std::string hdf5_temp_path;
//    static std::string hdf5_final_path;
    static void remove_temp_file();


};

