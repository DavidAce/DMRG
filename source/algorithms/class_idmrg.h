//
// Created by david on 2018-01-18.
//

#pragma once

#include "class_algorithm_infinite.h"
class class_log_infinite_dmrg_measurements;


/*!
 * \brief Class that runs the infinite DMRG algorithm.
 */
class class_idmrg : public class_algorithm_infinite {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_infinite::class_algorithm_infinite;
    explicit class_idmrg(std::shared_ptr<h5pp::File> h5ppFile_);
    StateRitz ritz = StateRitz::SR;
    void single_iDMRG_step();
    void run_simulation()                                   final;
    void check_convergence()                                final;
    bool    algo_on()                                       final;
    long      chi_lim_max()                                       final;
    bool      chi_lim_grow()                                      final;
    long      chi_lim_init()                                      final;
    size_t  print_freq()                                    final;
};

