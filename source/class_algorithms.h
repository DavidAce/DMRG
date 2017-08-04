//
// Created by david on 7/30/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_ALGORITHMS_H
#define FINITE_DMRG_EIGEN_CLASS_ALGORITHMS_H


#include <class_superblock.h>
#include <class_storage.h>
#include <class_tic_toc.h>

struct parameters{
    //Parmaters that control eigensolver and SVD precision
    int             chi                 = 20;           /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
    bool            increasing_chi      = false;
    int             chi_max             = 50;
    int             eigSteps            = 1000;         /*!< Maximum number of steps for eigenvalue solver. */
    double          eigThreshold        = 1e-10;        /*!< Minimum threshold for halting eigenvalue solver. */
    double          SVDThreshold        = 1e-12;        /*!< Minimum threshold value for keeping singular values. */
    double          delta_t             = 0.005;        /*!< Time step for iTEBD time evolution.*/
    //Parameters controlling simulation size and length.
    int             max_idmrg_length    = 50;           /*!< Final length of 1D quantum chain. */
    int             max_fdmrg_sweeps    = 4;            /*!< Number sweeps along the 1D quantum chain. */
    int             max_itebd_steps     = 5000;         /*!< Number of iTEBD iterations. */
    //Profiling
    bool            profiling           = false;
    int             time_prec           = 5;
    //Verbosity level
    int             verbosity           = 0;            /*!< Three-level verbosity */

};

class class_algorithms {
public:
    class_algorithms(){};

    parameters params;
    void iDMRG(class_superblock &superblock, class_storage &storage);
    void fDMRG(class_superblock &superblock, class_storage &storage);
    void iTEBD(class_superblock &superblock);


};


#endif //FINITE_DMRG_EIGEN_CLASS_ALGORITHMS_H
