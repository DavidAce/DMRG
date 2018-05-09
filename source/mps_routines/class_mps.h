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


class class_mps {
public:
    using Scalar = std::complex<double>;
private:
    long local_dimension;                          /*!< Local (or physical) spin or qubit dimension, usually denoted \f$d\f$ elsewhere. */
    Textra::Tensor<Scalar,3> tmp3;                 /*!< Temporary holder for swapping*/
    Textra::Tensor<Scalar,1> tmp1;                 /*!< Temporary holder for swapping*/
public:
    bool swapped = false;                                  /*!< Tracks the swapped state of A and B positions. */
    double truncation_error = 1;

    Textra::Tensor<Scalar,3> GA;                  /*!< \f$\Gamma^A\f$*/
    Textra::Tensor<Scalar,3> GB;                  /*!< \f$\Gamma^B\f$*/
    Textra::Tensor<Scalar,1> LA;                  /*!< \f$\Lambda^A\f$*/
    Textra::Tensor<Scalar,1> LB;                  /*!< \f$\Lambda^B\f$*/
    Textra::Tensor<Scalar,1> LB_left;              /*!< \f$\Lambda^B_{n+1}\f$ or \f$\Lambda^B_{n-1}\f$ in iDMRG or fDMRG respectively.*/


    Textra::Tensor<Scalar,4> theta,theta_evn_normalized, theta_odd_normalized;
    Textra::Tensor<Scalar,4> theta_sw ;
    Textra::Tensor<Scalar,3> LBGA, LAGB;
    Textra::Tensor<Scalar,2> l_evn, r_evn;
    Textra::Tensor<Scalar,2> l_odd, r_odd ;

    Textra::Tensor<Scalar,4> transfer_matrix_LBGA;
    Textra::Tensor<Scalar,4> transfer_matrix_LAGB;
    Textra::Tensor<Scalar,4> transfer_matrix_evn;
    Textra::Tensor<Scalar,4> transfer_matrix_odd;
//    Textra::Tensor<T,4> transfer_matrix_thetaL;
//    Textra::Tensor<T,4> transfer_matrix_thetaR_Sw;





    class_mps() = default;

    void initialize(long local_dimension_);         /*!< Sets local dimension*/
    void swap_AB();                                 /*!< Swaps the roles of A and B. Used in infinite DMRG.*/
    void compute_mps_components();


    Textra::Tensor<Scalar,4> get_theta(Scalar norm = 1.0) const;              /*!< Returns rank 4 tensor \f$\Theta\f$.*/
    Textra::Tensor<Scalar,4> get_theta_swapped(Scalar norm = 1.0) const;      /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
    Textra::Tensor<Scalar,4> get_theta_evn(Scalar norm = 1.0) const;          /*!< Returns rank 4 tensor \f$\Theta\f$.*/
    Textra::Tensor<Scalar,4> get_theta_odd(Scalar norm = 1.0) const;          /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
    Textra::Tensor<Scalar,3> A() const;
    Textra::Tensor<Scalar,3> B() const;
    Textra::Tensor<Scalar,4> get_transfer_matrix_zero() const;
    Textra::Tensor<Scalar,4> get_transfer_matrix_LBGA(Scalar norm = 1.0)const;
    Textra::Tensor<Scalar,4> get_transfer_matrix_GALA(Scalar norm = 1.0)const;
    Textra::Tensor<Scalar,4> get_transfer_matrix_GBLB(Scalar norm = 1.0)const;
    Textra::Tensor<Scalar,4> get_transfer_matrix_LAGB(Scalar norm = 1.0)const;
    Textra::Tensor<Scalar,4> get_transfer_matrix_theta_evn(Scalar norm  = 1.0)const;
    Textra::Tensor<Scalar,4> get_transfer_matrix_theta_odd(Scalar norm  = 1.0)const;
    Textra::Tensor<Scalar,4> get_transfer_matrix_AB(int p)const;
private:

//    Textra::Tensor<T,4> get_regularization_fixpointA()const;
//    Textra::Tensor<T,4> get_regularization_fixpointB()const;
//    Textra::Tensor<T,4> get_transfer_matrix_LBGA_regularized()const;
//    Textra::Tensor<T,4> get_transfer_matrix_GALA_regularized()const;
//    Textra::Tensor<T,4> get_transfer_matrix_GBLB_regularized()const;
//    Textra::Tensor<T,4> get_transfer_matrix_LAGB_regularized()const;
//    Textra::Tensor<T,4> get_transfer_matrix_AB_regularized()const;
//    Textra::Tensor<T,4> get_transfer_matrix_BA_regularized()const;
//    Textra::Tensor<T,4> get_transfer_matrix_AB_regularized_term(int p)const;
//    Textra::Tensor<T,4> get_transfer_matrix_regularized_inverseA()const;
//    Textra::Tensor<T,4> get_transfer_matrix_regularized_inverseB()const;
//    Textra::Tensor<T,4> get_transfer_matrix_regularized_inverseAB()const;
//    Textra::Tensor<T,4> get_transfer_matrix_regularized_inverseBA()const;

};
//    Textra::Tensor<T,4> get_transfer_matrix_L()const;
//    Textra::Tensor<T,4> get_transfer_matrix_R()const;
//    Textra::Tensor<T,4> get_transfer_2_site_matrix_L()const;
//    Textra::Tensor<T,4> get_transfer_2_site_matrix_R()const;


//
//
//template<typename T>
//Textra::Tensor<T,6> get_transfer_matrix_L(const Textra::Tensor<T,4> &MPO_1site)const{
//    return  A().template cast<T>()
//            .contract(MPO_1site, Textra::idx<1>({0},{2}))
//            .contract(A().conjugate(), Textra::idx<1>({4},{0}))
//            .shuffle(Textra::array6{0,4,2,1,5,3});
//}
//
//template<typename T>
//Textra::Tensor<T,6> get_transfer_matrix_R(const Textra::Tensor<T,4> &MPO_1site)const{
//    return B()
//            .contract(MPO_1site, Textra::idx<1>({0},{2}))
//            .contract(B().conjugate(), Textra::idx<1>({4},{0}))
//            .shuffle(Textra::array6{1,5,3,0,4,2});
//};
//
//template<typename T>
//Textra::Tensor<T,6> get_transfer_2_site_matrix_L(const Textra::Tensor<T,4> &MPO_1site)const{
//    return get_thetaL()
//            .contract(MPO_1site,            Textra::idx<1>({0}  ,{2}))
//            .contract(MPO_1site,            Textra::idx<2>({4,0},{0,2}))
//            .contract(get_thetaL().conjugate(), Textra::idx<2>({3,5},{0,1}))
//            .shuffle(Textra::array6{0,4,2,1,5,3});
////                .shuffle(array6{1,5,3,0,4,2});
//};
//template<typename T>
//Textra::Tensor<T,6> get_transfer_2_site_matrix_R(const Textra::Tensor<T,4> &MPO_1site)const{
//    return get_theta_evn()
//            .contract(MPO_1site,            Textra::idx<1>({0}  ,{2}))
//            .contract(MPO_1site,            Textra::idx<2>({4,0},{0,2}))
//            .contract(get_theta_evn().conjugate(), Textra::idx<2>({3,5},{0,1}))
//            .shuffle(Textra::array6{1,5,3,0,4,2});
////                .shuffle(array6{0,4,2,1,5,3});
//};
//
//template<typename T>
//Textra::Tensor<T,6> get_transfer_2_site_matrix_L(const Textra::Tensor<T,6> &MPO_2site)const{
//    return get_thetaL()
//            .contract(MPO_2site,            Textra::idx<2>({0,1},{2,3}))
//            .contract(get_thetaL().conjugate(), Textra::idx<2>({4,5},{0,1}))
//            .shuffle(Textra::array6{0,4,2,1,5,3});
//};
//
//template<typename T>
//Textra::Tensor<T,6> get_transfer_2_site_matrix_R(const Textra::Tensor<T,6> &MPO_2site)const{
//    return get_theta_evn()
//            .contract(MPO_2site,            Textra::idx<2>({0,1},{2,3}))
//            .contract(get_theta_evn().conjugate(), Textra::idx<2>({4,5},{0,1}))
//            .shuffle(Textra::array6{1,5,3,0,4,2});
//};

#endif //DMRG_CLASS_MPS_H
