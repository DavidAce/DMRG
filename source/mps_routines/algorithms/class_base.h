//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_DMRG_BASE_H
#define DMRG_CLASS_DMRG_BASE_H


#include <IO/class_hdf5.h>
#include <general/class_tic_toc.h>
#include <mps_routines/class_observables.h>

class class_superblock;

class class_base {
protected:
    class_base(class_hdf5 *hdf5_):hdf5(hdf5_){};
public:
    class_base(){};
    bool hdf5_is_valid;
    class_tic_toc t_tot;
    class_tic_toc t_eig;
    class_tic_toc t_sim;
    class_tic_toc t_svd;
    class_tic_toc t_env;
    class_tic_toc t_evo;
    class_tic_toc t_upd;
    class_tic_toc t_sto;
    class_tic_toc t_mps;
    class_hdf5  *hdf5;
    class_custom_cout ccout;
    virtual void run() = 0;
    virtual void run(class_superblock &superblock) = 0;

    void single_DMRG_step(class_superblock &superblock, long chi_max);
    void single_TEBD_step(class_superblock &superblock, long chi_max);
};


#endif //DMRG_CLASS_DMRG_BASE_H
