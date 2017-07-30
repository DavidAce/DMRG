/*! \file */

#include <class_tic_toc.h>
#include <n_tensor_extra.h>
#include <class_superblock.h>
#include <class_storage.h>

using namespace std;
using namespace Eigen;
using namespace Textra;




void iDMRG(class_superblock &superblock, class_storage &S, int max_length);
void fDMRG(class_superblock &superblock, class_storage &S, int sweeps);
void iTEBD(class_superblock &superblock, int max_iter);





int main() {
    /*!
        \brief  Main function
        \return an integer 0 upon exit success
    */

    int     chi          = 200;         /*! \param int chi  Keep \f$ \chi \f$  singular values in SVD truncation. */
    int     L            = 200;         /*! \param int L    Final length of 1D quantum chain. */
    int     N            = 10000;       /*! \param int N    Number of iTEBD iterations. */
    int     sweeps       = 4;
    double  SVDThreshold = 1e-12;
    double  eigThreshold = 1e-12;
    int     eigSteps     = 1000;

    //Initialize the superblock and storage for each chain length
    class_superblock superblock(eigSteps,eigThreshold,SVDThreshold,chi);
    class_storage S(L);


    iDMRG(superblock,S,L);
    fDMRG(superblock, S, sweeps);
    superblock.reset();
    iTEBD(superblock, N);
    superblock.print_error_TEBD();
    return 0;
}



void iDMRG(class_superblock &superblock, class_storage &S, int max_length){
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



void fDMRG(class_superblock &superblock, class_storage &S, int sweeps){
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



void iTEBD(class_superblock &superblock, int max_iter){
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




