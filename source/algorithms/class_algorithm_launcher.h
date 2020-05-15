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
    void setLogger(const std::string& name);
public:

    std::shared_ptr<h5pp::File> h5pp_file;

    explicit class_algorithm_launcher(std::shared_ptr<h5pp::File> h5ppFile_);
    class_algorithm_launcher();

    void start_h5pp_file();
    void setup_temp_path();
    static void clean_up();

    void run_algorithms();
    void run_idmrg();
    void run_fdmrg();
    void run_xdmrg();
    void run_itebd();


};

