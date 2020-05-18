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

                           LA  GA LC GB LB
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
    std::unique_ptr<class_mps_site> MPS_A;
    std::unique_ptr<class_mps_site> MPS_B;
    bool                            swapped = false; /*!< Tracks the swapped state of A and B positions. */
    public:

    class_mps_2site();
    ~class_mps_2site();                                           // Read comment on implementation
    class_mps_2site(class_mps_2site &&other) noexcept;            // default move ctor
    class_mps_2site &operator=(class_mps_2site &&other) noexcept; // default move assign
    class_mps_2site(const class_mps_2site &other);                // copy ctor
    class_mps_2site &operator=(const class_mps_2site &other);     // copy assign


    void set_mps(const class_mps_site & mpsA, const class_mps_site & mpsB);
    void set_mps(const Eigen::Tensor<Scalar, 3> &A, const Eigen::Tensor<Scalar, 1> &LC_, const Eigen::Tensor<Scalar, 3> &B);
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
    [[nodiscard]] Eigen::DSizes<long, 3>          dimensions() const;
    void                                          initialize(long spin_dim); /*!< Initializes the MPS*/
    void                                          swap_AB();                 /*!< Swaps the roles of A and B. Used in infinite DMRG.*/
    void                                          set_positions(size_t posA, size_t posB);
    void                                          set_positionA(size_t pos);
    void                                          set_positionB(size_t pos);
    [[nodiscard]] const class_mps_site &          get_mps_siteA() const;
    [[nodiscard]] const class_mps_site &          get_mps_siteB() const;
    [[nodiscard]] class_mps_site &                get_mps_siteA();
    [[nodiscard]] class_mps_site &                get_mps_siteB();
    [[nodiscard]] std::pair<size_t, size_t>       get_positions() const ;
    [[nodiscard]] size_t                          get_positionA() const ;
    [[nodiscard]] size_t                          get_positionB() const ;
    [[nodiscard]] Eigen::Tensor<Scalar, 3>        get_2site_tensor(Scalar norm = 1.0) const;
};
