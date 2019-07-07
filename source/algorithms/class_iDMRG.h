//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_INFINITE_DMRG_H
#define DMRG_CLASS_INFINITE_DMRG_H


#include "class_algorithm_infinite.h"
class class_log_dmrg;

/*!
 * \brief Class that runs the infinite DMRG algorithm.
 */
class class_iDMRG : public class_algorithm_infinite {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_infinite::class_algorithm_infinite;
    explicit class_iDMRG(std::shared_ptr<h5pp::File> h5ppFile_);
    std::unique_ptr<class_hdf5_log<class_log_dmrg>> log_dmrg;

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
};


#endif //DMRG_CLASS_INFINITE_DMRG_H
