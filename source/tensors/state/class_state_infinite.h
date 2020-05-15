//
// Created by david on 7/22/17.
//

#pragma once

#include <config/nmspc_settings.h>
#include <general/class_tic_toc.h>
#include <general/nmspc_tensor_extra.h>
#include <io/nmspc_logger.h>
#include <math/nmspc_eigutils.h>
#include <measure/state_measure_infinite.h>
#include <memory>
#include <optional>
class class_mps_2site;

/*!
  \class class_infinite_state
  \brief This class contains the Hamiltonian MPO, current wave function MPS, left and right environment blocks and routines to contract, diagonalize, truncate
   and update them.
*/

class class_state_infinite {
    private:
    std::optional<long> chi_lim;
    std::optional<long> chi_max;

    public:
    using Scalar = std::complex<double>;
    explicit class_state_infinite();

    std::unique_ptr<class_mps_2site> MPS; /*!< A matrix product state for two sites , A and B, and a center bond. In Vidal Canonical Form \f$\Lambda^A\Gamma^A
                                             \Lambda^C \Gamma^B\Lambda^B\f$. */

    bool is_real() const;
    bool has_nan() const;
    void assert_validity() const;

    long                     get_chi() const;
    long                     get_chi_lim() const;
    void                     set_chi_lim(long chi_lim_);
    long                     get_chi_max() const;
    void                     set_chi_max(long chi_max_);
    double                   get_truncation_error() const;
    Eigen::Tensor<Scalar, 3> get_mps() const;
    Eigen::DSizes<long, 3>   dimensions() const;
    void                     assert_positions() const;
    void                     enlarge_environment(int direction = 0);

    void set_mps(const Eigen::Tensor<Scalar,3> &MA,
                 const Eigen::Tensor<Scalar,1> &LC,
                 const Eigen::Tensor<Scalar,3> &MB);
    void set_mps(const Eigen::Tensor<Scalar,1> &LA,
                 const Eigen::Tensor<Scalar,3> &MA,
                 const Eigen::Tensor<Scalar,1> &LC,
                 const Eigen::Tensor<Scalar,3> &MB,
                 const Eigen::Tensor<Scalar,1> &LB);
    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> get_H_local_matrix() const;
    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> get_H_local_sq_matrix() const;
    void                                             set_positions(size_t position);
    void                                             swap_AB(); /*!< Swap the roles of A and B. Used in the infinite-DMRG stage.*/

    mutable state_measure_infinite measurements;
    mutable double                 lowest_recorded_variance = 1.0;

    void do_all_measurements() const;
    void clear_measurements() const;
    void clear_cache() const;

    private:
    struct Cache {
        std::optional<Eigen::Tensor<Scalar, 3>> mps = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 4>> theta = std::nullopt;
        // Possibly more views?
    };
    mutable Cache cache;
};
