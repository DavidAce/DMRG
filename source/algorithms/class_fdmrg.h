//
// Created by david on 2018-01-31.
//

#pragma once
#include "class_algorithm_finite.h"
class class_h5table_measurements_finite;

/*!
// * \brief Class that runs the finite DMRG algorithm.
 */

class class_state_finite;
class class_fdmrg : public class_algorithm_finite {
    public:
    // Inherit the constructor of class_algorithm_base
    using class_algorithm_finite::class_algorithm_finite;
    explicit class_fdmrg(std::shared_ptr<h5pp::File> h5pp_file_);
    StateRitz ritz = StateRitz::SR;
    void      single_fdmrg_step();
    void      resume() final;
    void      run_task_list(std::deque<fdmrg_task> &task_list);
    void      run_default_task_list() final;
    void      run_preprocessing() final;
    void      run_algorithm() final;
    void      check_convergence() final;
    bool      cfg_algorithm_is_on() final;
    long      cfg_chi_lim_max() final;
    size_t    cfg_print_freq() final;
    bool      cfg_chi_lim_grow() final;
    long      cfg_chi_lim_init() final;
    bool      cfg_store_wave_function() final;
};
