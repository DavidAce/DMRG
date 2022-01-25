#pragma once
#include <complex>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

class TensorsInfinite;
class StateInfinite;
class ModelInfinite;
class EdgesInfinite;
namespace tools::infinite::measure {
    using Scalar = std::complex<double>;
    /* clang-format off */
    extern void   do_all_measurements             (const TensorsInfinite & tensors);
    extern void   do_all_measurements             (const StateInfinite & state);

    extern size_t length                          (const TensorsInfinite & tensors);
    extern size_t length                          (const EdgesInfinite & edges);
    extern long   bond_dimension                  (const StateInfinite & state);
    extern double truncation_error                (const StateInfinite & state);
    extern double norm                            (const StateInfinite & state);
    extern double entanglement_entropy            (const StateInfinite & state);

    template<typename state_or_mps_type>
    double energy_minus_energy_shift     (const state_or_mps_type & state, const ModelInfinite & model, const EdgesInfinite & edges);
    template<typename state_or_mps_type>
    double energy_mpo                      (const state_or_mps_type & state, const ModelInfinite & model, const EdgesInfinite & edges);
    template<typename state_or_mps_type>
    double energy_per_site_mpo             (const state_or_mps_type & state, const ModelInfinite & model, const EdgesInfinite & edges);
    template<typename state_or_mps_type>
    double energy_variance_mpo             (const state_or_mps_type & state, const ModelInfinite & model, const EdgesInfinite & edges);
    template<typename state_or_mps_type>
    double energy_variance_per_site_mpo    (const state_or_mps_type & state, const ModelInfinite & model, const EdgesInfinite & edges);

    extern double energy_mpo                      (const TensorsInfinite & tensors);
    extern double energy_per_site_mpo             (const TensorsInfinite & tensors);
    extern double energy_variance_mpo             (const TensorsInfinite & tensors);
    extern double energy_variance_per_site_mpo    (const TensorsInfinite & tensors);

    extern double energy_mpo                      (const Eigen::Tensor<Scalar,3> &mps, const TensorsInfinite & tensors);
    extern double energy_per_site_mpo             (const Eigen::Tensor<Scalar,3> &mps, const TensorsInfinite & tensors);
    extern double energy_variance_mpo             (const Eigen::Tensor<Scalar,3> &mps, const TensorsInfinite & tensors);
    extern double energy_variance_per_site_mpo    (const Eigen::Tensor<Scalar,3> &mps, const TensorsInfinite & tensors);

    extern double energy_per_site_ham             (const TensorsInfinite & tensors);
    extern double energy_per_site_mom             (const TensorsInfinite & tensors);
    extern double energy_variance_per_site_ham    (const TensorsInfinite & tensors);
    extern double energy_variance_per_site_mom    (const TensorsInfinite & tensors);

    /* clang-format on */

}