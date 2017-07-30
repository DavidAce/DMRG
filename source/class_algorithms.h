//
// Created by david on 7/30/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_ALGORITHMS_H
#define FINITE_DMRG_EIGEN_CLASS_ALGORITHMS_H


#include <class_superblock.h>
#include <class_storage.h>
#include <class_tic_toc.h>


class class_algorithms {
private:

// Profiling objects
    class_profiling t_svd = class_profiling::class_profiling(1,5, std::string("SVD           ")) ;
    class_profiling t_eig = class_profiling::class_profiling(1,5, std::string("Diagonalize   ")) ;
    class_profiling t_env = class_profiling::class_profiling(1,5, std::string("Update Env.   ")) ;
    class_profiling t_tmp = class_profiling::class_profiling(1,5, std::string("Temporary     ")) ;
    class_profiling t_tot = class_profiling::class_profiling(1,5, std::string("Total         ")) ;

public:
    class_algorithms(){};


    void iDMRG(class_superblock &superblock, class_storage &S, int max_length);
    void fDMRG(class_superblock &superblock, class_storage &S, int sweeps);
    void iTEBD(class_superblock &superblock, int max_iter);


};


#endif //FINITE_DMRG_EIGEN_CLASS_ALGORITHMS_H
