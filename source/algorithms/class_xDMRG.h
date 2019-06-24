//
// Created by david on 2018-02-09.
//

#ifndef DMRG_CLASS_EXITED_DMRG_H
#define DMRG_CLASS_EXITED_DMRG_H

#include "class_algorithm_finite.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Core>





/*!
 * \brief Class that runs the excited-state DMRG algorithm.
 */

class class_finite_chain_state;
class class_xDMRG : public class_algorithm_finite {
private:

public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_finite::class_algorithm_finite;
    explicit class_xDMRG(std::shared_ptr<h5pp::File> h5ppFile_);

    bool   projected_during_saturation  = false;


    void run_simulation()                                                   final;
    void run_preprocessing()                                                final;
    void single_DMRG_step();                                                final;
    void find_energy_range();
    void check_convergence()                                                final;

    bool   sim_on()                                     final;

    long   chi_max()                                                        final;
    size_t num_sites()                                                      final;
    size_t store_freq()                                                     final;
    size_t print_freq()                                                     final;
    bool   chi_grow()                                                       final;

};



#endif //DMRG_CLASS_EXITED_DMRG_H
