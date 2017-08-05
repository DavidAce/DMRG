//
// Created by david on 7/30/17.
//

#include "class_algorithms.h"



void class_algorithms::iDMRG(class_superblock &superblock, class_storage &storage){
/*!
 * \fn void iDMRG(class_superblock &superblock, class_storage &S, int max_length)
 * \brief Infinite DMRG, grows the chain from 2 up to `max_idmrg_length` particles.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param storage A class that stores current MPS and environments at each iteration.
 * \param max_length Maximum chain length after which the algorithm stops.
 */
    class_tic_toc t_tot(params.profiling, params.time_prec, "iDMRG Total Time           ");
    class_tic_toc t_eig(params.profiling, params.time_prec, "iDMRG Eigenvalue solver    ");
    class_tic_toc t_svd(params.profiling, params.time_prec, "iDMRG SVD Truncation       ");
    class_tic_toc t_env(params.profiling, params.time_prec, "iDMRG Enlarge environment  ");
    class_tic_toc t_sto(params.profiling, params.time_prec, "iDMRG Store MPS            ");
    class_tic_toc t_mps(params.profiling, params.time_prec, "iDMRG Update MPS           ");

    storage.set_length(params.max_idmrg_length);
    t_tot.tic();
    int length = 0;
    while(length < params.max_idmrg_length){
        length += 2;
                        superblock.update_bond_dimensions();
        t_eig.tic();    superblock.find_ground_state(params.eigSteps, params.eigThreshold);     t_eig.toc();
        t_svd.tic();    superblock.truncate         (params.chi     , params.SVDThreshold);     t_svd.toc();
        t_mps.tic();    superblock.update_MPS();                                                t_mps.toc();
        t_sto.tic();    storage.store_insert(superblock);                                       t_sto.toc();
        t_env.tic();    superblock.enlarge_environment();                                       t_env.toc();
                        superblock.print_picture(params.graphics);
        superblock.print_state(params.verbosity);
                        superblock.swap_AB();
    }
    t_tot.toc();
    superblock.print_state(params.verbosity + 1);
    t_eig.print_time_w_percent();
    t_svd.print_time_w_percent();
    t_env.print_time_w_percent();
    t_sto.print_time_w_percent();
    t_mps.print_time_w_percent();
    t_tot.print_time();
}



void class_algorithms::fDMRG(class_superblock &superblock, class_storage &storage){
/*!
 * \fn void fDMRG(class_superblock &superblock, class_storage &S, int sweeps)
 * \brief Finite DMRG sweeps across the chain built during iDMRG.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param storage A class that stores current MPS and environments at each iteration.
 * \param sweeps Maximum number of sweeps.
 */
    class_tic_toc t_tot(params.profiling, params.time_prec, "fDMRG Total Time           ");
    class_tic_toc t_eig(params.profiling, params.time_prec, "fDMRG Eigenvalue solver    ");
    class_tic_toc t_svd(params.profiling, params.time_prec, "fDMRG SVD Truncation       ");
    class_tic_toc t_env(params.profiling, params.time_prec, "fDMRG Enlarge environment  ");
    class_tic_toc t_sto(params.profiling, params.time_prec, "fDMRG Store MPS            ");
    class_tic_toc t_mps(params.profiling, params.time_prec, "fDMRG Update MPS           ");
    int direction  = 1;
    int sweep = 0;
    t_tot.tic();
    Eigen::ArrayXi chi_list = Eigen::ArrayXi::LinSpaced(params.max_fdmrg_sweeps,
                                                        params.chi,
                                                        params.increasing_chi?
                                                        std::max(params.chi_max, params.chi+10)
                                                        : params.chi);

    while(sweep < params.max_fdmrg_sweeps) {
                        int chi = chi_list(sweep);
                        storage.load(superblock);
                        superblock.update_bond_dimensions();
        t_eig.tic();    superblock.find_ground_state(params.eigSteps, params.eigThreshold);     t_eig.toc();
        t_svd.tic();    superblock.truncate         (chi            , params.SVDThreshold);     t_svd.toc();
        t_mps.tic();    superblock.update_MPS();                                                t_mps.toc();
        t_sto.tic();    storage.overwrite_MPS(superblock);                                      t_sto.toc();
                        superblock.print_picture(params.graphics);
        superblock.print_state(params.verbosity);
        t_env.tic();    superblock.enlarge_environment(direction);                              t_env.toc();
        t_sto.tic();    storage.move(superblock, direction, sweep);                             t_sto.toc();
    }

    t_tot.toc();
    superblock.print_state(params.verbosity + 1);
    t_eig.print_time_w_percent();
    t_svd.print_time_w_percent();
    t_env.print_time_w_percent();
    t_sto.print_time_w_percent();
    t_mps.print_time_w_percent();
    t_tot.print_time();
}



void class_algorithms::iTEBD(class_superblock &superblock){
/*!
 * \fn iTEBD(class_superblock &superblock, int max_iter)
 * \brief infinite Time evolving block decimation.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param max_iter Maximum number of iterations.
 */
    class_tic_toc t_evo(params.profiling, params.time_prec, "iTEBD Time evolution       ");
    class_tic_toc t_svd(params.profiling, params.time_prec, "iTEBD SVD Truncation       ");
    class_tic_toc t_mps(params.profiling, params.time_prec, "iTEBD Update MPS           ");
    class_tic_toc t_tot(params.profiling, params.time_prec, "iTEBD Total Time           ");
    superblock.reset();
    t_tot.tic();
    for(auto steps = 0; steps < params.max_itebd_steps ; steps++){
                        superblock.update_bond_dimensions();
        t_evo.tic();    superblock.time_evolve()                               ; t_evo.toc();
        t_svd.tic();    superblock.truncate   (params.chi, params.SVDThreshold); t_svd.toc();
        t_mps.tic();    superblock.update_MPS();                                 t_mps.toc();
                        superblock.swap_AB();
    }
    t_tot.toc();
    superblock.print_state(params.verbosity+1);
    t_evo.print_time_w_percent();
    t_svd.print_time_w_percent();
    t_mps.print_time_w_percent();
    t_tot.print_time();
}



