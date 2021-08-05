#pragma once

#include "AlgorithmInfinite.h"

/*!
 * \brief Class that runs the infinite DMRG algorithm.
 */
class idmrg : public AlgorithmInfinite {
    public:
    // Inherit the constructor of class_algorithm_base
    using AlgorithmInfinite::AlgorithmInfinite;
    explicit idmrg(std::shared_ptr<h5pp::File> h5ppFile_);
    OptRitz   ritz = OptRitz::SR;
    void      single_iDMRG_step();
    void      run_simulation() final;
    void      check_convergence() final;
};
