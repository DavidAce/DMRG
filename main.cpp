/* iDMRG code to find the ground state of the infinite 1D XXZ Ising chain.
 * This implementation follows the steps in the article
 *      “Phase Diagram of the Anisotropic Spin-2 XXZ Model: Infinite-System Density Matrix Renormalization Group Study.”
 *      Kjäll, Zaletel, Mong, Bardarson, and Pollmann.
 *      Physical Review B 87 (23): 235106. doi:10.1103/PhysRevB.87.235106.
 *
 *
 *
 */


#include <class_tic_toc.h>
#include <n_tensor_extra.h>
#include <class_TwoSiteHamiltonian.h>
#include <class_SVD.h>
#include <class_storage.h>
#include <class_MPS.h>
#include <class_superblock.h>

using namespace std;
using namespace Eigen;
using namespace Textra;



// First define the parameters of the simulation

long chi = 100;
int N = 50;
double SVDThreshold = 1e-12;
double eigThreshold = 1e-12;
int    eigSteps = 1000;

// Profiling objects
class_profiling t_svd   (1,5, string("SVD           ")) ;
class_profiling t_eig   (1,5, string("Diagonalize   ")) ;
class_profiling t_env   (1,5, string("Update Env.   ")) ;
class_profiling t_tmp   (1,5, string("Temporary     ")) ;
class_profiling t_tot   (1,5, string("Total         ")) ;



void infinite_DMRG(class_superblock &superblock, class_storage &S, long max_length){
//    t_tot.tic();
    long length = 0;
    while(length < max_length){
        superblock.print_picture();
        superblock.update_bond_dimensions();


        superblock.find_ground_state();
        superblock.truncate();
        superblock.update_MPS();
        S.store_insert(superblock);

        superblock.enlarge_environment();
        superblock.print_error();

        superblock.swap_AB();
        length += 2;
    }
    S.print_storage();
//    t_tot.toc();
//    t_svd.print_total(t_tot.total_time);cout<<endl;
//    t_eig.print_total(t_tot.total_time);cout<<endl;
//    t_env.print_total(t_tot.total_time);cout<<endl;
//    t_tmp.print_total(t_tot.total_time);cout<<endl;
//    t_tot.print_total(); cout << endl;
}



void finite_DMRG(class_superblock &superblock, class_storage &S, long sweeps){
//    t_tot.tic();
    long direction  = 1;
    long sweep = 0;

    while(sweep < sweeps) {
//    for (long sweep = 0; sweep < sweeps; sweep++){
//        for (size_t step = 0; step < S.max_length - 2; step++) {

            S.load(superblock);
            superblock.update_bond_dimensions();

            superblock.print_picture();
            superblock.find_ground_state();
            superblock.truncate();
            superblock.update_MPS();
            superblock.print_error();

            S.overwrite(superblock);
//            superblock.enlarge_environment(direction);
            S.move(superblock, direction);

            if (S.position_L <= 1 || S.position_R >= S.max_length - 1) {
                direction *= -1;
            }
            if(S.position_L == S.max_length/2 -1 && S.position_R == S.max_length/2){
                sweep++;
            }
//        }
    }

    S.print_storage();
//    t_tot.toc();
//    t_svd.print_total(t_tot.total_time);cout<<endl;
//    t_eig.print_total(t_tot.total_time);cout<<endl;
//    t_env.print_total(t_tot.total_time);cout<<endl;
//    t_tmp.print_total(t_tot.total_time);cout<<endl;
//    t_tot.print_total(); cout << endl;
}


int main() {
    class_superblock superblock(eigSteps,eigThreshold,SVDThreshold,chi);
    size_t L = 10;
    class_storage S(L);
    infinite_DMRG(superblock,S,L);
    finite_DMRG(superblock,S, 4);
    return 0;
}

