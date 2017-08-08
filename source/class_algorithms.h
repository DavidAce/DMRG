//
// Created by david on 7/30/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_ALGORITHMS_H
#define FINITE_DMRG_EIGEN_CLASS_ALGORITHMS_H


#include <class_superblock.h>
#include <class_storage.h>
#include <class_tic_toc.h>
#include <class_hdf5.h>
#include <n_settings.h>


namespace s = settings;


class class_algorithms {
public:
    class_algorithms(){};

//    settings params;
    void iDMRG(class_superblock &superblock, class_storage &storage, class_hdf5 &hdf5);
    void fDMRG(class_superblock &superblock, class_storage &storage, class_hdf5 &hdf5);
    void iTEBD(class_superblock &superblock, class_hdf5 &hdf5);
    void FES(class_superblock &superblock, class_hdf5 &hdf5);


};


#endif //FINITE_DMRG_EIGEN_CLASS_ALGORITHMS_H
