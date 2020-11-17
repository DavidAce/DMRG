
#pragma once
#include "class_algorithm_infinite.h"
class class_h5table_tebd;

/*!
 * \brief Class that runs the imaginary TEBD algorithm.
 */
class class_itebd : public class_algorithm_infinite {
    public:
    using class_algorithm_infinite::class_algorithm_infinite;
    explicit class_itebd(std::shared_ptr<h5pp::File> h5ppFile_);

    std::vector<Eigen::Tensor<Scalar, 2>> unitary_time_evolving_operators;
    Eigen::MatrixXcd                      h_evn;
    Eigen::MatrixXcd                      h_odd;

    void   single_TEBD_step();
    void   run_simulation() final;
    void   run_preprocessing() final;
    void   run_postprocessing() final;
    void   check_convergence_time_step();
    void   check_convergence() final;
    bool   cfg_algorithm_is_on() final;
    long   cfg_chi_lim_max() final;
    size_t cfg_print_freq() final;
    bool   cfg_chi_lim_grow() final;
    long   cfg_chi_lim_init() final;
};
