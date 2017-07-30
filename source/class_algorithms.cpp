//
// Created by david on 7/30/17.
//

#include "class_algorithms.h"



void class_algorithms::iDMRG(class_superblock &superblock, class_storage &S, int max_length){
/*!
 * \fn void iDMRG(class_superblock &superblock, class_storage &S, int max_length)
 * \brief Infinite DMRG, grows the chain from 2 up to `max_length` particles.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param S A class that stores current MPS and environments at each iteration.
 * \param max_length Maximum chain length after which the algorithm stops.
 */
    int length = 0;
    while(length < max_length){
        superblock.print_picture();
        superblock.update_bond_dimensions();

        superblock.find_ground_state();
        superblock.truncate();
        superblock.update_MPS();
        S.store_insert(superblock);

        superblock.enlarge_environment();
        superblock.print_error_DMRG();

        superblock.swap_AB();
        length += 2;
    }
//    S.print_storage();
//    t_tot.toc();
//    t_svd.print_total(t_tot.total_time);cout<<endl;
//    t_eig.print_total(t_tot.total_time);cout<<endl;
//    t_env.print_total(t_tot.total_time);cout<<endl;
//    t_tmp.print_total(t_tot.total_time);cout<<endl;
//    t_tot.print_total(); cout << endl;
}



void class_algorithms::fDMRG(class_superblock &superblock, class_storage &S, int sweeps){
/*!
 * \fn void fDMRG(class_superblock &superblock, class_storage &S, int sweeps)
 * \brief Finite DMRG sweeps across the chain built during iDMRG.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param S A class that stores current MPS and environments at each iteration.
 * \param sweeps Maximum number of sweeps.
 */
//    t_tot.tic();
    int direction  = 1;
    int sweep = 0;

    while(sweep < sweeps) {
        S.load(superblock);
        superblock.update_bond_dimensions();

        superblock.print_picture();
        superblock.find_ground_state();
        superblock.truncate();
        superblock.update_MPS();
        superblock.print_error_DMRG();
        S.overwrite_MPS(superblock);

        superblock.enlarge_environment(direction);
        S.move(superblock, direction);

        if (S.position_L <= 1 || S.position_R >= S.max_length - 1) {
            direction *= -1;
        }
        if(S.position_L == S.max_length/2 -1 && S.position_R == S.max_length/2){
            sweep++;
        }
//        }
    }

//    S.print_storage();
//    t_tot.toc();
//    t_svd.print_total(t_tot.total_time);cout<<endl;
//    t_eig.print_total(t_tot.total_time);cout<<endl;
//    t_env.print_total(t_tot.total_time);cout<<endl;
//    t_tmp.print_total(t_tot.total_time);cout<<endl;
//    t_tot.print_total(); cout << endl;
}



void class_algorithms::iTEBD(class_superblock &superblock, int max_iter){
/*!
 * \fn iTEBD(class_superblock &superblock, int max_iter)
 * \brief infinite Time evolving block decimation.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param max_iter Maximum number of iterations.
 */
    for(auto iter = 0; iter < max_iter ; iter++){
        superblock.update_bond_dimensions();
        superblock.time_evolve();
        superblock.truncate();
        superblock.update_MPS();
        superblock.swap_AB();
    }
    superblock.print_error_TEBD();
}



