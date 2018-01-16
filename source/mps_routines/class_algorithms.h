//
// Created by david on 7/30/17.
//

#ifndef DMRG_CLASS_ALGORITHMS_H
#define DMRG_CLASS_ALGORITHMS_H
//#define EIGEN_USE_MKL_ALL


#include <general/class_tic_toc.h>
#include <IO/class_hdf5.h>
#include <sim_parameters/n_sim_settings.h>

class class_superblock;

class class_algorithms  {
private:
    class_tic_toc t_tot;
    class_tic_toc t_eig;
    class_tic_toc t_sim;
    class_tic_toc t_svd;
    class_tic_toc t_env;
    class_tic_toc t_evo;
    class_tic_toc t_upd;
    class_tic_toc t_sto;
    class_tic_toc t_mps;
    class_hdf5 hdf5;
    class_custom_cout ccout;

public:
    class_algorithms(){};
    void single_DMRG_step(class_superblock &superblock, long chi_max);
    void single_TEBD_step(class_superblock &superblock, long chi_max);

    void iDMRG();
    void fDMRG();
    void iTEBD();
    void FES_iTEBD();
    void FES_iDMRG();
};


#endif //DMRG_CLASS_ALGORITHMS_H
