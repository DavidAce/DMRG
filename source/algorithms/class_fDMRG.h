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
class class_fDMRG : public class_algorithm_finite {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_finite::class_algorithm_finite;
    explicit class_fDMRG(std::shared_ptr<h5pp::File> h5ppFile_);

    void run_simulation()                                        final;
    void check_convergence()                                     final;
    bool   sim_on()                                              final;
    long   chi_max()                                             final;
    size_t num_sites()                                           final;
    size_t write_freq()                                          final;
    size_t print_freq()                                          final;
    bool   chi_grow()                                            final;
    long   chi_init()                                            final;
    bool   store_wave_function()                                 final;

};


