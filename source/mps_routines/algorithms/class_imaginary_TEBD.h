//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_IMAGINARY_TEBD_H
#define DMRG_CLASS_IMAGINARY_TEBD_H
#include "class_base.h"

/*!
// * \fn iTEBD(class_superblock &superblock, class_hdf5 &hdf5)
// * \brief infinite Time evolving block decimation.
// * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
// * \param max_iter Maximum number of iterations.
// */

class class_imaginary_TEBD :public class_base {
public:
    using class_base::class_base;
    explicit class_imaginary_TEBD(std::shared_ptr<class_hdf5_file> hdf5_);
    long chi_max    = settings::itebd::chi_max;
    long max_steps  = settings::itebd::max_steps;
    int  print_freq = settings::itebd::print_freq;
    int  iteration  = 0;
    double phys_time = 0;

    void run() override;

    void print_status_full() override;
    void print_status_update() override;

    void store_table_entry() override;
    void print_profiling()   override;
};


#endif //DMRG_CLASS_IMAGINARY_TEBD_H
