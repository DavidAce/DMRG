#ifndef DMRG_CLASS_MPS_H
#define DMRG_CLASS_MPS_H

#include "general/nmspc_tensor_extra.h"


/*!
 \class class_mps
 \brief Contains the Matrix Product State (MPS) for two sites, A and B, corresponding to sites \f$n\f$ and \f$n+1\f$, respectively.

 \section description Description
 During the initial infinite-DMRG stage, sites A and B are swapped in each iteration to mimic the movement between even and odd positions
  in a translationally invariant chain.

 \subsection background Background
 In both the infinite and finite DMRG algorithms, the MPS in this class will be used to construct  \f$\Theta\f$, a two-site MPS tensor of rank 4:

 \f[
  \Theta = \Lambda^B_{n-1} \Gamma^A_{n} \Lambda^A_{n} \Gamma^B_{n+1} \Lambda^B_{n+1}.
 \f]

 > Note that the first \f$ \Lambda^B_{n-1} \f$  is special:
 > - During infinite-DMRG \f$ \Lambda^B_{n-1} = \Lambda^B_{n+1} \f$ because we expect the MPS to be symmetric about \f$n\f$.
 > - During finite-DMRG it is simply the bond directly to the left of \f$n\f$.

 \subsection diagram Diagram

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

 Tensor index order convention: See the <a href="index.html">Home page</a>.
*/


template<typename Scalar>
class class_mps {
public:
private:
    long local_dimension;                          /*!< Local (or physical) spin or qubit dimension, usually denoted \f$d\f$ elsewhere. */
    Textra::Tensor<Scalar,3> tmp3;                         /*!< Temporary holder for swapping*/
    Textra::Tensor<Scalar,1> tmp1;                         /*!< Temporary holder for swapping*/
public:
    bool swapped;                                  /*!< Tracks the swapped state of A and B positions. */

    Textra::Tensor<Scalar,3> GA;                  /*!< \f$\Gamma^A\f$*/
    Textra::Tensor<Scalar,3> GB;                  /*!< \f$\Gamma^B\f$*/
    Textra::Tensor<Scalar,1> LA;                  /*!< \f$\Lambda^A\f$*/
    Textra::Tensor<Scalar,1> LB;                  /*!< \f$\Lambda^B\f$*/
    Textra::Tensor<Scalar,1> L_tail;              /*!< \f$\Lambda^B_{n+1}\f$ or \f$\Lambda^B_{n-1}\f$ in iDMRG or fDMRG respectively.*/

    class_mps(){};

    void initialize(long local_dimension_);         /*!< Sets local dimension*/
    void swap_AB();                                 /*!< Swaps the roles of A and B. Used in infinite DMRG.*/
    Textra::Tensor<Scalar,3> A() const;
    Textra::Tensor<Scalar,3> B() const;
    Textra::Tensor<Scalar,4> thetaL() const;
    Textra::Tensor<Scalar,4> thetaR() const;

    Textra::Tensor<Scalar,4> get_theta() const ;             /*!< Returns rank 4 tensor \f$\Theta\f$.*/
    Textra::Tensor<Scalar,4> get_theta_swapped() const ;     /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
    Textra::Tensor<Scalar,4> get_transfer_matrix_L()const;
    Textra::Tensor<Scalar,4> get_transfer_matrix_R()const;
    Textra::Tensor<Scalar,4> get_transfer_2_site_matrix_L()const;
    Textra::Tensor<Scalar,4> get_transfer_2_site_matrix_R()const;

    template<typename T>
    Textra::Tensor<T,6> get_transfer_matrix_L(const Textra::Tensor<T,4> &MPO_1site)const{
        return  A().template cast<T>()
                .contract(MPO_1site, Textra::idx<1>({0},{2}))
                .contract(A().conjugate(), Textra::idx<1>({4},{0}))
                .shuffle(Textra::array6{0,4,2,1,5,3});
    }

    template<typename T>
    Textra::Tensor<T,6> get_transfer_matrix_R(const Textra::Tensor<T,4> &MPO_1site)const{
        return B()
                .contract(MPO_1site, Textra::idx<1>({0},{2}))
                .contract(B().conjugate(), Textra::idx<1>({4},{0}))
                .shuffle(Textra::array6{1,5,3,0,4,2});
    };

    template<typename T>
    Textra::Tensor<T,6> get_transfer_2_site_matrix_L(const Textra::Tensor<T,4> &MPO_1site)const{
        return thetaL()
                .contract(MPO_1site,            Textra::idx<1>({0}  ,{2}))
                .contract(MPO_1site,            Textra::idx<2>({4,0},{0,2}))
                .contract(thetaL().conjugate(), Textra::idx<2>({3,5},{0,1}))
                .shuffle(Textra::array6{0,4,2,1,5,3});
//                .shuffle(array6{1,5,3,0,4,2});
    };
    template<typename T>
    Textra::Tensor<T,6> get_transfer_2_site_matrix_R(const Textra::Tensor<T,4> &MPO_1site)const{
        return thetaR()
                .contract(MPO_1site,            Textra::idx<1>({0}  ,{2}))
                .contract(MPO_1site,            Textra::idx<2>({4,0},{0,2}))
                .contract(thetaR().conjugate(), Textra::idx<2>({3,5},{0,1}))
                .shuffle(Textra::array6{1,5,3,0,4,2});
//                .shuffle(array6{0,4,2,1,5,3});
    };

    template<typename T>
    Textra::Tensor<T,6> get_transfer_2_site_matrix_L(const Textra::Tensor<T,6> &MPO_2site)const{
        return thetaL()
                .contract(MPO_2site,            Textra::idx<2>({0,1},{2,3}))
                .contract(thetaL().conjugate(), Textra::idx<2>({4,5},{0,1}))
                .shuffle(Textra::array6{0,4,2,1,5,3});
    };

    template<typename T>
    Textra::Tensor<T,6> get_transfer_2_site_matrix_R(const Textra::Tensor<T,6> &MPO_2site)const{
        return thetaR()
                .contract(MPO_2site,            Textra::idx<2>({0,1},{2,3}))
                .contract(thetaR().conjugate(), Textra::idx<2>({4,5},{0,1}))
                .shuffle(Textra::array6{1,5,3,0,4,2});
    };

};





#endif //DMRG_CLASS_MPS_H
