//
// Created by david on 2018-01-18.
//

#pragma once

#include "class_algorithm_infinite.h"
class class_log_tebd;

/*!
 * \brief Class that runs the imaginary TEBD algorithm.
 */
class class_iTEBD :public class_algorithm_infinite {
public:
    using class_algorithm_infinite::class_algorithm_infinite;
    explicit class_iTEBD(std::shared_ptr<h5pp::File> h5ppFile_);

    std::unique_ptr<class_hdf5_log<class_log_tebd>> log_tebd;
    std::vector<Eigen::Tensor<Scalar,4>> unitary_time_evolving_operators;
    Eigen::MatrixXcd h_evn;
    Eigen::MatrixXcd h_odd;



    void run_simulation()                                   final;
    void run_preprocessing()                                final;
    void run_postprocessing()                               final;
    void single_TEBD_step();
    void check_convergence_time_step();
    void check_convergence()                                final;
    void write_logs(bool force = false)                     final;
    bool   sim_on ()                                        final;
    long   chi_max()                                        final;
    size_t num_sites()                                      final;
    size_t write_freq()                                     final;
    size_t print_freq()                                     final;
    bool   chi_grow()                                       final;
};

