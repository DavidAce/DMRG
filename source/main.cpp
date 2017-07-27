
/*************************************************************************************************/

/*! \mainpage
 * \brief This program finds the ground state of a 1D quantum Ising chain in a transverse field using the DMRG algorithm.

  # DMRG (in development)
  [Density matrix renormalization group](https://en.wikipedia.org/wiki/Density_matrix_renormalization_group) (DMRG) is a variational numerical technique to study the low-energy physics of many-body quantum systems.

  This algorithm constructs and minimizes trial wave functions, in the shape of [Matrix Product States](https://en.wikipedia.org/wiki/Matrix_product_state) (MPS), iteratively in order to find the ground state of one-dimensional quantum systems with high precision.

  This implementation loosely follows the steps outlined in:

  > [Phase Diagram of the Anisotropic Spin-2 XXZ Model: Infinite-System Density Matrix Renormalization Group Study.](https://arxiv.org/abs/1212.6255)<br>
  > by Kjäll, Zaletel, Mong, Bardarson, and Pollmann. Physical Review B 87 (23): 235106. <br>

  > [Efficient Numerical Simulations Using Matrix-Product States](http://quantumtensor.pks.mpg.de/wp-content/uploads/2016/06/notes_1.pdf)<br>
  > by Frank Pollmann. <br>

  > [The density-matrix renormalization group in the age of matrix product states](https://arxiv.org/abs/1008.3477)<br>
  > by Ulrich Schollwöck. <br>


 ## Notation

 The *Vidal canonical form*, i.e. \f$\Gamma\Lambda\Gamma\Lambda\f$"..., is used throughout this code.
 In code we denote

 - \f$\Gamma \rightarrow\f$ `G`.
 - \f$\Lambda \rightarrow\f$ `L`.

 ## Tensor index order convention.
 The tensor index order used here follows the convention:
 - physical indices first, from left to right or for MPO's, up to down.
 - other dimensions (like bond dimensions) next, from left to right.

 #### Example:
 Consider for some position \f$n\f$ on the chain \f$\Gamma = \Gamma^{\sigma_n}_{a,b}\f$.
Here \f$\sigma_n \in [-1,1]\f$ is a particle with local (physical) dimension \f$d\f$ = 2, and \f$a,b\f$ are the remaining dimensions, in practice they are
bond dimensions of \f$\Lambda^{n-1}\f$ and \f$\Lambda^{n}\f$, respectively, which can be numbers \f$\in [1,\chi]\f$.

In diagrammatic tensor notation this is:
@verbatim
                 	    [d]          0
            G     =	[a]__|__[b] = 1__|__2
@endverbatim
where after the second equality the index order is shown. In code this corresponds to

\code{.cpp}
 Textra::Tensor3 G(d,a,b);
\endcode

Similarly, we have for \f$\Theta^{\sigma_n,\sigma_{n+1}}_{\chi_a,\chi_b}\f$:

@verbatim
                 	           	[d] [d]                0   1
            Theta     =	[chia]___|___|___[chib] = 2 ___|___|___ 3
@endverbatim

which in code reads

\code{.cpp}
 Textra::Tensor4 G(d,d,chia,chib);
\endcode

# Requirements
The following software is required and has been included:
* [Eigen](http://eigen.tuxfamily.org)
* [Spectra](https://spectralib.org/)

# Details
 \author    David Aceituno
 \version   1.0
 \date      07-2017
 \copyright MPL2

*/

/*************************************************************************************************/


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

    S.print_storage();
//    t_tot.toc();
//    t_svd.print_total(t_tot.total_time);cout<<endl;
//    t_eig.print_total(t_tot.total_time);cout<<endl;
//    t_env.print_total(t_tot.total_time);cout<<endl;
//    t_tmp.print_total(t_tot.total_time);cout<<endl;
//    t_tot.print_total(); cout << endl;
}



