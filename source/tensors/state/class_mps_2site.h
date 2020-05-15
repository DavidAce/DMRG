#pragma once

#include "general/nmspc_tensor_extra.h"
#include <memory>

class class_mps_site;

/*!
 \class class_mps_2site
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

class class_mps_2site {
    public:
    using Scalar = std::complex<double>;

    private:
    //    long spin_dimension = 2;                      /*!< Local (or physical) spin or qubit dimension, usually denoted \f$d\f$ elsewhere. */
    Eigen::Tensor<Scalar, 3> tmp3; /*!< Temporary holder for swapping*/
    Eigen::Tensor<Scalar, 1> tmp1; /*!< Temporary holder for swapping*/

    template<class T>
    std::unique_ptr<T> copy_unique(const std::unique_ptr<T> &source) {
        return source ? std::make_unique<T>(*source) : nullptr;
    }

    public:
    //    explicit class_mps_2site(std::string model_type_str);
    class_mps_2site() = default;
    explicit class_mps_2site(const class_mps_2site &other);

    bool   swapped          = false; /*!< Tracks the swapped state of A and B positions. */
    double truncation_error = 0;

    std::unique_ptr<class_mps_site> MPS_A;
    std::unique_ptr<class_mps_site> MPS_B;
    void                            set_mps(const Eigen::Tensor<Scalar, 3> &A, const Eigen::Tensor<Scalar, 1> &LC_, const Eigen::Tensor<Scalar, 3> &B);
    void set_mps(const Eigen::Tensor<Scalar, 1> &LA, const Eigen::Tensor<Scalar, 3> &A, const Eigen::Tensor<Scalar, 1> &LC_, const Eigen::Tensor<Scalar, 3> &B,
                 const Eigen::Tensor<Scalar, 1> &LB);

    void                                          assert_validity() const;
    [[nodiscard]] bool                            is_real() const;
    [[nodiscard]] bool                            has_nan() const;
    [[nodiscard]] long                            chiA() const;
    [[nodiscard]] long                            chiB() const;
    [[nodiscard]] long                            chiC() const;
    [[nodiscard]] long                            spin_dim_A() const;
    [[nodiscard]] long                            spin_dim_B() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &A_bare() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &A() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &B() const;
    [[nodiscard]] Eigen::Tensor<Scalar, 2>        LC() const;
    [[nodiscard]] Eigen::Tensor<Scalar, 3>        GA() const;
    [[nodiscard]] Eigen::Tensor<Scalar, 3>        GB() const;
    [[nodiscard]] Eigen::Tensor<Scalar, 2>        LA() const;
    [[nodiscard]] Eigen::Tensor<Scalar, 2>        LB() const;

    Eigen::DSizes<long, 3>   dimensions() const;
    void                     initialize(int spin_dim);         /*!< Initializes the MPS*/
    void                     swap_AB();                        /*!< Swaps the roles of A and B. Used in infinite DMRG.*/
    Eigen::Tensor<Scalar, 3> get_mps(Scalar norm = 1.0) const; /*!< Returns rank 3 tensor having contracted A and B, with physical leg size d² \f$\Theta\f$.*/
};
