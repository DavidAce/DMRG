
#pragma once
#include "AlgorithmInfinite.h"
#include <unsupported/Eigen/CXX11/Tensor>

/*!
 * \brief Class that runs the imaginary TEBD algorithm.
 */
class itebd : public AlgorithmInfinite {
    public:
    using AlgorithmInfinite::AlgorithmInfinite;
    explicit itebd(std::shared_ptr<h5pp::File> h5ppFile_);

    std::vector<Eigen::Tensor<Scalar, 2>> unitary_time_evolving_operators;
    Eigen::Tensor<Scalar, 2>              h_evn, h_odd;

    void single_TEBD_step();
    void run_simulation() final;
    void run_preprocessing() final;
    void run_postprocessing() final;
    void check_convergence_time_step();
    void check_convergence() final;
};
