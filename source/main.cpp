/*! \mainpage

  # DMRG
  This program finds the ground state of a 1D quantum Ising chain in a transverse field.
  The implementation follows the steps in the article

  > [Phase Diagram of the Anisotropic Spin-2 XXZ Model: Infinite-System Density Matrix Renormalization Group Study.”](https://arxiv.org/abs/1212.6255)<br>
  > by Kjäll, Zaletel, Mong, Bardarson, and Pollmann. Physical Review B 87 (23): 235106. <br>


 */


#include <class_tic_toc.h>
#include <n_tensor_extra.h>
#include <class_superblock.h>
#include <class_storage.h>

using namespace std;
using namespace Eigen;
using namespace Textra;




// Profiling objects
class_profiling t_svd   (1,5, string("SVD           ")) ;
class_profiling t_eig   (1,5, string("Diagonalize   ")) ;
class_profiling t_env   (1,5, string("Update Env.   ")) ;
class_profiling t_tmp   (1,5, string("Temporary     ")) ;
class_profiling t_tot   (1,5, string("Total         ")) ;

void infinite_DMRG  (class_superblock &superblock, class_storage &S, int max_length);
void finite_DMRG    (class_superblock &superblock, class_storage &S, int sweeps);

int main() {
    // First define the parameters of the simulation

    int    chi          = 40;
    int    L            = 50;
    double SVDThreshold = 1e-12;
    double eigThreshold = 1e-12;
    int    eigSteps     = 1000;

    class_superblock superblock(eigSteps,eigThreshold,SVDThreshold,chi);
    class_storage S(L);
    infinite_DMRG(superblock,S,L);
    finite_DMRG(superblock,S, 4);
    return 0;
}



void infinite_DMRG(class_superblock &superblock, class_storage &S, int max_length){
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



void finite_DMRG(class_superblock &superblock, class_storage &S, int sweeps){
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



