//
// Created by david on 7/30/17.
//

#ifndef DMRG_CLASS_ALGORITHMS_H
#define DMRG_CLASS_ALGORITHMS_H


#include <class_tic_toc.h>
#include <class_hdf5.h>
#include <n_settings.h>

class class_superblock;

class class_algorithms  {
public:
    class_tic_toc t_tot;
    class_tic_toc t_eig;
    class_tic_toc t_svd;
    class_tic_toc t_env;
    class_tic_toc t_evo;
    class_tic_toc t_upd;
    class_tic_toc t_sto;
    class_tic_toc t_mps;

    class_hdf5 hdf5;


public:
    class_algorithms():
    hdf5(class_hdf5(settings::hdf5::filename, settings::hdf5::path, true)){};
    void single_DMRG_step(class_superblock &superblock, long chi_max);
    void single_TEBD_step(class_superblock &superblock, long chi_max);

    void iDMRG();
    void fDMRG();
    void iTEBD();
    void FES_iTEBD();
    void FES_iDMRG();

private:



};


#endif //DMRG_CLASS_ALGORITHMS_H
