#ifndef FINITE_DMRG_EIGEN_CLASS_MPS_H
#define FINITE_DMRG_EIGEN_CLASS_MPS_H

#include <n_tensor_extra.h>

using namespace Textra;
using namespace std;

/*! \brief Contains the Matrix Product State (MPS) for two sites, A and B, corresponding to sites \f$n\f$ and \f$n+1\f$, respectively.*/

/*!
 # MPS Class
 During the initial infinite-DMRG stage, sites A and B are swapped in each iteration to mimic the movement between even and odd positions
 in a translationally invariant chain.

 ## Background
 In both the infinite and finite DMRG algorithms, the MPS in this class will be used to construct  \f$\Theta\f$, a two-site MPS tensor of rank 4:

 \f[
  \Theta = \Lambda^B_{n-1} \Gamma^A_{n} \Lambda^A_{n} \Gamma^B_{n+1} \Lambda^B_{n+1}.
 \f]

 > Note that the first \f$ \Lambda^B_{n-1} \f$  is special:
 > - During infinite-DMRG \f$ \Lambda^B_{n-1} = \Lambda^B_{n+1} \f$ because we expect the MPS to be symmetric about \f$n\f$.
 > - During finite-DMRG it is simply the bond directly to the left of \f$n\f$.

 ### Diagram

 In diagrammatic tensor notation the two-site MPS is:
@verbatim
                 	           	  [d]    [d]
            Theta     =	[chia]__.__|__.__|__.__[chib]

                 	           LB  GA LA GB LB
@endverbatim
 where the values in [ ] denote the dimension of each leg of the full tensor:
 - \f$d\f$ = local (or physical) dimension. Default is 2 for spin-1/2, Ising spins or qubits.
 - \f$\chi_a\f$ = Left dimension equivalent to bond dimension of \f$\Lambda^B_{n-1}\f$
 - \f$\chi_b\f$ = Right dimension equivalent to bond dimension of \f$\Lambda^B_{n+1}\f$


 ### Tensor index order convention: See the <a href="index.html">Home page</a>.
*/


class class_MPS {
private:
    long local_dimension;                           /*!< Local (or physical) spin or qubit dimension, usually denoted \f$d\f$ elsewhere. */
    bool swap;                                      /*!< Tracks the swapped state of A and B positions. */
    Tensor3 tmp3;                                   /*!< Temporary holder for swapping*/
    Tensor1 tmp1;                                   /*!< Temporary holder for swapping*/
public:

    Tensor3 GA;                                     /*!< \f$\Gamma^A\f$*/
    Tensor3 GB;                                     /*!< \f$\Gamma^B\f$*/
    Tensor1 LA;                                     /*!< \f$\Lambda^A\f$*/
    Tensor1 LB;                                     /*!< \f$\Lambda^B\f$*/
    Tensor1 L_tail;                                 /*!< \f$\Lambda^B_{n+1}\f$ or \f$\Lambda^B_{n-1}\f$ in iDMRG or fDMRG respectively.*/

    class_MPS(){}
    void initialize(const long local_dimension_);   /*! Sets local dimension*/
    Tensor4 get_theta() const;                      /*! Returns rank 4 tensor \f$\Theta\f$.*/
    void swap_AB();                                 /*! Swaps the roles of A and B. Used in infinite DMRG.*/
    double get_energy(const Tensor4 &Hamiltonian);   /*! Computes the current energy by contracting the current MPS with the Hamiltonian MPO.*/
};


#endif //FINITE_DMRG_EIGEN_CLASS_MPS_H
