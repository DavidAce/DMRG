//
// Created by david on 7/22/17.
//

#pragma once

#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <list>
#include <measure/state_measure_infinite.h>
#include <memory>
#include <optional>

class class_mps_site;

/*!
  \class class_infinite_state
  \brief This class contains the current 2-site translationally invariant wave function in MPS form
*/

class class_state_infinite {
    public:
    using Scalar = std::complex<double>;

    private:
    struct Cache {
        std::optional<Eigen::Tensor<Scalar, 3>> twosite_mps     = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 4>> theta           = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 3>> GA              = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 3>> GB              = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 2>> LC_diag         = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 2>> LA_diag         = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 2>> LB_diag         = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 2>> LC_diag_inv     = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 2>> LA_diag_inv     = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 2>> LB_diag_inv     = std::nullopt;

    };

    std::unique_ptr<class_mps_site> MPS_A;
    std::unique_ptr<class_mps_site> MPS_B;
    bool                            swapped = false; /*!< Tracks the swapped state of A and B positions. */

    //    std::unique_ptr<class_mps_2site> mps_sites; /*!< A matrix product state for two sites , A and B, and a center bond. In Vidal Canonical Form
    //                                                   \f$\Lambda^A\Gamma^A \Lambda^C \Gamma^B\Lambda^B\f$. */
    std::optional<long> chi_lim;
    std::optional<long> chi_max;
    mutable Cache       cache;

    public:
    mutable state_measure_infinite measurements;
    mutable double                 lowest_recorded_variance = 1.0;

    public:
    class_state_infinite();
    ~class_state_infinite();                                                // Read comment on implementation
    class_state_infinite(class_state_infinite &&other);                     // default move ctor
    class_state_infinite &operator=(class_state_infinite &&other);          // default move assign
    class_state_infinite(const class_state_infinite &other);                // copy ctor
    class_state_infinite &operator=(const class_state_infinite &other);     // copy assign

    void initialize(ModelType model_type);

    void                                          assert_validity() const;
    [[nodiscard]] bool                            is_real() const;
    [[nodiscard]] bool                            has_nan() const;
    [[nodiscard]] double                          get_truncation_error() const;
    [[nodiscard]] std::pair<size_t, size_t>       get_positions();
    [[nodiscard]] size_t                          get_positionA();
    [[nodiscard]] size_t                          get_positionB();
    [[nodiscard]] long                            chiC() const;
    [[nodiscard]] long                            chiA() const;
    [[nodiscard]] long                            chiB() const;
    [[nodiscard]] long                            get_chi_lim() const;
    [[nodiscard]] long                            get_chi_max() const;
    [[nodiscard]] long                            get_spin_dimA() const;
    [[nodiscard]] long                            get_spin_dimB() const;
    [[nodiscard]] Eigen::DSizes<long, 3>          dimensions() const;
    [[nodiscard]] const class_mps_site &          get_mps_siteA() const;
    [[nodiscard]] const class_mps_site &          get_mps_siteB() const;
    [[nodiscard]] class_mps_site &                get_mps_siteA();
    [[nodiscard]] class_mps_site &                get_mps_siteB();
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &A_bare() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &A() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &B() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 2> &LC_diag() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 2> &LA_diag() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 2> &LB_diag() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 2> &LC_diag_inv() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 2> &LA_diag_inv() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 2> &LB_diag_inv() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> & GA() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> & GB() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 1> & LC() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 1> & LA() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 1> & LB() const;

    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &get_2site_mps(Scalar norm = 1.0) const;
    void                                          set_chi_lim(long chi_lim_);
    void                                          set_chi_max(long chi_max_);
    void                                          set_positions(size_t position);

    void swap_AB(); /*!< Swap the roles of A and B. Used in the infinite-DMRG stage.*/
    void set_mps(const Eigen::Tensor<Scalar, 3> &twosite_tensor);
    void set_mps(const std::list<class_mps_site> &mps_list);
    void set_mps(const class_mps_site &mpsA, const class_mps_site &mpsB);
    void set_mps(const Eigen::Tensor<Scalar, 3> &MA, const Eigen::Tensor<Scalar, 1> &LC, const Eigen::Tensor<Scalar, 3> &MB);
    void set_mps(const Eigen::Tensor<Scalar, 1> &LA, const Eigen::Tensor<Scalar, 3> &MA, const Eigen::Tensor<Scalar, 1> &LC, const Eigen::Tensor<Scalar, 3> &MB,
                 const Eigen::Tensor<Scalar, 1> &LB);

    //    template<typename T>
    //    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> get_H_local_matrix() const;
    //    template<typename T>
    //    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> get_H_local_sq_matrix() const;

    void do_all_measurements() const;
    void clear_measurements() const;
    void clear_cache() const;
};
