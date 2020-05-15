
#pragma once
#include <general/nmspc_tensor_omp.h>
#include "class_algorithm_infinite.h"
class class_h5table_tebd;

/*!
 * \brief Class that runs the imaginary TEBD algorithm.
 */
class class_itebd :public class_algorithm_infinite {
public:
    using class_algorithm_infinite::class_algorithm_infinite;
    explicit class_itebd(std::shared_ptr<h5pp::File> h5ppFile_);

    std::vector<Eigen::Tensor<Scalar,2>> unitary_time_evolving_operators;
    Eigen::MatrixXcd h_evn;
    Eigen::MatrixXcd h_odd;



    void run_simulation()                                   final;
    void run_preprocessing()                                final;
    void run_postprocessing()                               final;
    void single_TEBD_step();
    void check_convergence_time_step();
    void check_convergence()                                final;
    bool   algo_on()                                        final;
    long   chi_max()                                        final;
//    size_t write_freq()                                     final;
    size_t print_freq()                                     final;
    bool   chi_grow()                                       final;
    long   chi_init()                                       final;
};

