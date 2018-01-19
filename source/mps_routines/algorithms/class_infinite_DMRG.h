//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_INFINITE_DMRG_H
#define DMRG_CLASS_INFINITE_DMRG_H

/*! \brief  Class that runs the infinite DMRG algorithm.*/


#include <general/class_tic_toc.h>
#include "class_base.h"

/*!
 * # infinite DMRG
 *
 * \fn \void iDMRG(class_superblock &superblock, class_storage &S, int max_length)
 * \brief Infinite DMRG, grows the chain from 2 up to `max_idmrg::length` particles.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param storage A class that stores current MPS and environments at each iteration.
 * \param max_length Maximum chain length after which the algorithm stops.
 */

class class_infinite_DMRG : public class_base {
private:
    long chi_max     = settings::idmrg::chi_max;
    int  max_length  = settings::idmrg::max_length;

public:
    class_infinite_DMRG() = default;
    explicit class_infinite_DMRG(class_hdf5 &hdf5_):class_base(&hdf5_){
    }

    void run() override;
    void run(class_superblock &superblock) override;

    void set(long chi_max_, int max_length_){
        chi_max     = chi_max_;
        max_length = max_length_;
    }

    void set(long chi_max_, int max_length_, class_hdf5 &hdf5_){
        chi_max     = chi_max_;
        max_length = max_length_;
        hdf5 = &hdf5_;

    }
    void set(long chi_max_, class_hdf5 &hdf5_){
        chi_max     = chi_max_;
        hdf5 = &hdf5_;
    }

    void set(int max_length_, class_hdf5 &hdf5_){
        max_length = max_length_;
        hdf5 = &hdf5_;
    }

};


#endif //DMRG_CLASS_INFINITE_DMRG_H
