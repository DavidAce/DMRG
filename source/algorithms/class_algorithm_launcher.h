//
// Created by david on 7/30/17.
//
#pragma once

#include <memory>
namespace h5pp{class File;}

class class_algorithm_launcher  {
public:
    std::shared_ptr<h5pp::File> h5pp_file;

    explicit class_algorithm_launcher(std::shared_ptr<h5pp::File> h5ppFile_);
    class_algorithm_launcher();
    void start_h5pp_file();
    void setup_temp_path();

    void run_algorithms();
    void run_idmrg() const;
    void run_fdmrg() const;
    void run_flbit() const;
    void run_xdmrg() const;
    void run_itebd() const;


};

