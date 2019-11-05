//
// Created by david on 2018-01-18.
//

#pragma once

#include "class_algorithm_infinite.h"
class class_log_infinite_dmrg_measurements;


/*!
 * \brief Class that runs the infinite DMRG algorithm.
 */
class class_iDMRG : public class_algorithm_infinite {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_infinite::class_algorithm_infinite;
    explicit class_iDMRG(std::shared_ptr<h5pp::File> h5ppFile_);

    void single_DMRG_step(std::string ritz);
    void run_simulation()                                   final;
    void check_convergence()                                final;
    void write_logs(bool force = false)                     final;
    bool    sim_on ()                                       final;
    long    chi_max()                                       final;
    size_t  num_sites()                                     final;
    size_t  write_freq()                                    final;
    size_t  print_freq()                                    final;
    bool    chi_grow()                                      final;
    long    chi_init()                                      final;
};

