//
// Created by david on 2018-01-31.
//

#ifndef DMRG_CLASS_FINITE_DMRG_H
#define DMRG_CLASS_FINITE_DMRG_H

#include "class_algorithm_finite.h"


/*!
// * \brief Class that runs the finite DMRG algorithm.
 */

class class_finite_state;
class class_fDMRG : public class_algorithm_finite {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_finite::class_algorithm_finite;
    explicit class_fDMRG(std::shared_ptr<h5pp::File> h5ppFile_);
    bool   projected_during_saturation  = false;

    void run_simulation()                                        final;
    void check_convergence()                                     final;
    bool   sim_on()                                              final;
    long   chi_max()                                             final;
    size_t num_sites()                                           final;
    size_t store_freq()                                          final;
    size_t print_freq()                                          final;
    bool   chi_grow()                                            final;
    bool   store_wave_function()                                 final;

};



#endif //DMRG_CLASS_FINITE_DMRG_H
