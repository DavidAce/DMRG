#ifndef DMRG_CLASS_MPS_H
#define DMRG_CLASS_MPS_H

#include <memory>
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





class class_vidal_mps {
public:
    using Scalar = std::complex<double>;
    void set_mps(const Eigen::Tensor<Scalar,3> &G_, const Eigen::Tensor<Scalar,1> &L_){G = G_; L = L_;}
    void set_L(const Eigen::Tensor<Scalar,1> &L_){L=L_;}
    void set_G(const Eigen::Tensor<Scalar,3> &G_){G=G_;}
    const auto &get_G() const {return std::as_const(G);}
    const auto &get_L() const {return std::as_const(L);}
    class_vidal_mps() = default;
private:
    Eigen::Tensor<Scalar,3> G;                  /*!< \f$\Gamma \f$*/
    Eigen::Tensor<Scalar,1> L;                  /*!< \f$\Lambda\f$*/
};




class class_mps_2site {
public:
    using Scalar = std::complex<double>;
private:

    long spin_dimension;                         /*!< Local (or physical) spin or qubit dimension, usually denoted \f$d\f$ elsewhere. */
    Eigen::Tensor<Scalar,3> tmp3;                 /*!< Temporary holder for swapping*/
    Eigen::Tensor<Scalar,1> tmp1;                 /*!< Temporary holder for swapping*/

    template< class T >
    std::unique_ptr<T> copy_unique(const std::unique_ptr<T>& source)
    {
        return source ? std::make_unique<T>(*source) : nullptr;
    }
public:
    bool swapped = false;                                  /*!< Tracks the swapped state of A and B positions. */
    double truncation_error = 1;

    std::unique_ptr<class_vidal_mps> MPS_A;
    std::unique_ptr<class_vidal_mps> MPS_B;
    Eigen::Tensor<Scalar,1> LC;

    auto chiA () {return MPS_A->get_L().dimension(0);}
    auto chiB () {return MPS_B->get_L().dimension(0);}
    auto chiC () {return LC.dimension(0);}

    auto A() const{
        using namespace Textra;
        return asDiagonal(MPS_A->get_L()).contract(MPS_A->get_G(), idx({1},{1})).shuffle(array3{1,0,2});
    };

    auto B() const{
        using namespace Textra;
        return MPS_B->get_G().contract(asDiagonal(MPS_B->get_L()), idx({2},{0}));
    };

    auto C() const{
        using namespace Textra;
        return asDiagonal(LC);
    };

    void set_mps(const Eigen::Tensor<Scalar,1> &LA,
                 const Eigen::Tensor<Scalar,3> &GA,
                 const Eigen::Tensor<Scalar,1> &LC_,
                 const Eigen::Tensor<Scalar,3> &GB,
                 const Eigen::Tensor<Scalar,1> &LB){
        MPS_A->set_mps(GA,LA);
        MPS_B->set_mps(GB,LB);
        LC = LC_;
    }


    Eigen::Tensor<Scalar,4> theta,theta_evn_normalized, theta_odd_normalized;
    Eigen::Tensor<Scalar,4> theta_sw ;
    Eigen::Tensor<Scalar,3> LBGA, LAGB;
    Eigen::Tensor<Scalar,2> l_evn, r_evn;
    Eigen::Tensor<Scalar,2> l_odd, r_odd;

    Eigen::Tensor<Scalar,4> transfer_matrix_LBGA;
    Eigen::Tensor<Scalar,4> transfer_matrix_LAGB;
    Eigen::Tensor<Scalar,4> transfer_matrix_evn;
    Eigen::Tensor<Scalar,4> transfer_matrix_odd;


    class_mps_2site();
    class_mps_2site(const class_mps_2site &other);

    void initialize(int spin_dim);                              /*!< Initializes the MPS*/
    void swap_AB();                                             /*!< Swaps the roles of A and B. Used in infinite DMRG.*/
    void compute_mps_components();


    Eigen::Tensor<Scalar,4> get_theta(Scalar norm = 1.0) const;              /*!< Returns rank 4 tensor \f$\Theta\f$.*/
    Eigen::Tensor<Scalar,4> get_theta_swapped(Scalar norm = 1.0) const;      /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
    Eigen::Tensor<Scalar,4> get_theta_evn(Scalar norm = 1.0) const;          /*!< Returns rank 4 tensor \f$\Theta\f$.*/
    Eigen::Tensor<Scalar,4> get_theta_odd(Scalar norm = 1.0) const;          /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/

    Eigen::Tensor<Scalar,4> get_transfer_matrix_zero() const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_LBGA(Scalar norm = 1.0)const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_GALC(Scalar norm = 1.0)const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_GBLB(Scalar norm = 1.0)const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_LCGB(Scalar norm = 1.0)const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_theta_evn(Scalar norm  = 1.0)const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_theta_odd(Scalar norm  = 1.0)const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_AB(int p)const;
};


#endif //DMRG_CLASS_MPS_H
