#pragma once

#include <memory>
namespace h5pp {
    class File;
}

class AlgorithmLauncher {
    public:
    std::shared_ptr<h5pp::File> h5file;

    explicit AlgorithmLauncher(std::shared_ptr<h5pp::File> h5ppFile_);
    AlgorithmLauncher();
    void start_h5file();
    void setup_temp_path();

    void run_algorithms();
    void run_idmrg();
    void run_fdmrg();
    void run_flbit();
    void run_xdmrg();
    void run_itebd();
};
