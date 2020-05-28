//
// Created by david on 2018-02-09.
//

#pragma once

#include "class_algorithm_finite.h"
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

/*!
 * \brief Class that runs the excited-state DMRG algorithm.
 */

class class_state_finite;
class class_xdmrg : public class_algorithm_finite {
    private:
    double energy_window_growth_factor = 1.0;

    public:
    // Inherit the constructor of class_algorithm_base
    using class_algorithm_finite::class_algorithm_finite;
    explicit class_xdmrg(std::shared_ptr<h5pp::File> h5ppFile_);
    void   find_energy_range();
    void   single_xDMRG_step();
    void   single_xDMRG_step_old();
    void   reset_to_random_state_in_energy_window(ResetReason reason, std::optional<std::string> sector = std::nullopt);
    void   run_preprocessing() final;
    void   run_algorithm() final;
    void   check_convergence() final;
    bool   algo_on() final;
    long   chi_lim_max() final;
    size_t print_freq() final;
    bool   chi_lim_grow() final;
    long   chi_lim_init() final;
    bool   store_wave_function() final;
};
