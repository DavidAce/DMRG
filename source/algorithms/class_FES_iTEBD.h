//
// Created by david on 2018-01-31.
//

#ifndef DMRG_FES_TEBD_H
#define DMRG_FES_TEBD_H


#include "class_base_algorithm.h"
/*!
 * \brief Class that studies Finite-entanglement scaling using the imaginary TEBD algorithm.
 */
class class_FES_iTEBD : public class_base_algorithm {
public:
    //Inherit the constructor of class_base_algorithm
    using class_base_algorithm::class_base_algorithm;
    explicit class_FES_iTEBD(std::shared_ptr<class_hdf5_file> hdf5_);
    void run() override;
    void run2();
    void print_profiling()     override;
    void print_profiling_sim(class_tic_toc &t_parent) override;

};
#endif //DMRG_FES_TEBD_H
