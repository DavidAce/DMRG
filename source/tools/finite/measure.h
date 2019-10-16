//
// Created by david on 2019-10-13.
//

#pragma once

#include <vector>
#include <list>
#include <complex>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

class class_finite_state;
class class_model_base;
class class_simulation_status;

namespace tools::finite::measure{
    using Scalar = std::complex<double>;

//            extern void do_all_measurements                           (class_finite_state & state);
    extern int length                                         (const class_finite_state & state);
    extern size_t bond_dimension_current                      (const class_finite_state & state);
    extern size_t bond_dimension_midchain                     (const class_finite_state & state);
    extern std::vector<size_t> bond_dimensions                (const class_finite_state & state);
    extern double norm                                        (const class_finite_state & state);



    extern double energy                                      (const class_finite_state & state);
    extern double energy_per_site                             (const class_finite_state & state);
    extern double energy_variance                             (const class_finite_state & state);
    extern double energy_variance_per_site                    (const class_finite_state & state);
    extern double spin_component                              (const class_finite_state & state, Eigen::Matrix2cd paulimatrix);
    extern Eigen::Tensor<Scalar,1> mps_wavefn                 (const class_finite_state & state);
    extern double entanglement_entropy_current                (const class_finite_state & state);
    extern double entanglement_entropy_midchain               (const class_finite_state & state);
    extern std::vector<double> entanglement_entropies         (const class_finite_state & state);
    extern std::vector<double> spin_components                (const class_finite_state & state);

    namespace twosite{
        extern double energy_minus_energy_reduced                 (const class_finite_state & state, const Eigen::Tensor<Scalar,4> & theta);
        extern double energy                                      (const class_finite_state & state, const Eigen::Tensor<Scalar,4> & theta);
        extern double energy_per_site                             (const class_finite_state & state, const Eigen::Tensor<Scalar,4> & theta);
        extern double energy_variance                             (const class_finite_state & state, const Eigen::Tensor<Scalar,4> & theta);
        extern double energy_variance_per_site                    (const class_finite_state & state, const Eigen::Tensor<Scalar,4> & theta);
    }

    namespace multisite{
        namespace internal{
            inline double digits;
            double significant_digits(double H2, double E2);
        }
        extern double energy_minus_energy_reduced             (const class_finite_state & state, const Eigen::Tensor<Scalar,3> & multitheta);
        extern double energy                                  (const class_finite_state & state, const Eigen::Tensor<Scalar,3> & multitheta);
        extern double energy_per_site                         (const class_finite_state & state, const Eigen::Tensor<Scalar,3> & multitheta);
        extern double energy_variance                         (const class_finite_state & state, const Eigen::Tensor<Scalar,3> & multitheta);
        extern double energy_variance_per_site                (const class_finite_state & state, const Eigen::Tensor<Scalar,3> & multitheta);
        extern double energy                                  (const class_finite_state & state);
        extern double energy_per_site                         (const class_finite_state & state);
        extern double energy_variance                         (const class_finite_state & state);
        extern double energy_variance_per_site                (const class_finite_state & state);
    }

    template<typename Derived>
    double energy_minus_energy_reduced(const class_finite_state & state, const Eigen::TensorBase<Derived,Eigen::WriteAccessors> & theta){
        constexpr int rank = Derived::NumIndices;
        if constexpr (rank == 4) return twosite::energy_minus_energy_reduced(state,theta);
        if constexpr (rank == 3) return multisite::energy_minus_energy_reduced(state,theta);
        static_assert("Wrong rank, expected 3 or 4" and (rank == 3 or rank == 4));
    }
    template<typename Derived>
    double energy(const class_finite_state & state, const Eigen::TensorBase<Derived,Eigen::WriteAccessors> & theta){
        constexpr int rank = Derived::NumIndices;
        if constexpr (rank == 4) return twosite::energy(state,theta);
        if constexpr (rank == 3) return multisite::energy(state,theta);
        static_assert("Wrong rank, expected 3 or 4" and (rank == 3 or rank == 4));
    }
    template<typename Derived>
    double energy_per_site(const class_finite_state & state, const Eigen::TensorBase<Derived,Eigen::WriteAccessors> & theta){
        constexpr int rank = Derived::NumIndices;
        if constexpr (rank == 4) return twosite::energy_per_site(state,theta);
        if constexpr (rank == 3) return multisite::energy_per_site(state,theta);
        static_assert("Wrong rank, expected 3 or 4" and (rank == 3 or rank == 4));
    }
    template<typename Derived>
    double energy_variance(const class_finite_state & state, const Eigen::TensorBase<Derived,Eigen::WriteAccessors> & theta){
        constexpr int rank = Derived::NumIndices;
        if constexpr (rank == 4) return twosite::energy_variance(state,theta);
        if constexpr (rank == 3) return multisite::energy_variance(state,theta);
        static_assert("Wrong rank, expected 3 or 4" and (rank == 3 or rank == 4));
    }
    template<typename Derived>
    double energy_variance_per_site(const class_finite_state & state, const Eigen::TensorBase<Derived,Eigen::WriteAccessors> & theta){
        constexpr int rank = Derived::NumIndices;
        if constexpr (rank == 4) return twosite::energy_variance_per_site(state,theta);
        if constexpr (rank == 3) return multisite::energy_variance_per_site(state,theta);
        static_assert("Wrong rank, expected 3 or 4" and (rank == 3 or rank == 4));

    }
}

