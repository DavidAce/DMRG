//
// Created by david on 2018-01-18.
//

#include "class_base.h"
#include <mps_routines/class_superblock.h>
namespace s = settings;
using namespace std;
using namespace Textra;


void class_base::single_DMRG_step(class_superblock &superblock, long chi_max){
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */

    superblock.update_bond_dimensions();
    t_eig.tic();    superblock.find_ground_state(s::precision::eigSteps, s::precision::eigThreshold);    t_eig.toc();
    t_svd.tic();    superblock.truncate         (chi_max,                s::precision::SVDThreshold);    t_svd.toc();
    t_mps.tic();    superblock.update_MPS();                                                             t_mps.toc();

}

void class_base::single_TEBD_step(class_superblock &superblock, long chi_max){
/*!
 * \fn single_iTEBD_step(class_superblock &superblock)
 * \brief infinite Time evolving block decimation.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 */
    superblock.update_bond_dimensions();
    t_evo.tic();    superblock.time_evolve();                                       t_evo.toc();
    t_svd.tic();    superblock.truncate   (chi_max, s::precision::SVDThreshold);    t_svd.toc();
    t_mps.tic();    superblock.update_MPS();                                        t_mps.toc();
}