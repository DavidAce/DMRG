#pragma once
#include <complex>
#include <general/nmspc_tensor_omp.h>
#include <vector>

class class_tensors_infinite;
class class_state_infinite;
class class_model_infinite;
class class_edges_infinite;
namespace tools::infinite::measure{
    using Scalar = std::complex<double>;
    /* clang-format off */
    extern void   do_all_measurements             (const class_tensors_infinite & tensors);
    extern void   do_all_measurements             (const class_state_infinite & state);

    extern size_t length                          (const class_tensors_infinite & tensors);
    extern size_t length                          (const class_edges_infinite & edges);
    extern long   bond_dimension                  (const class_state_infinite & state);
    extern double truncation_error                (const class_state_infinite & state);
    extern double norm                            (const class_state_infinite & state);
    extern double entanglement_entropy            (const class_state_infinite & state);

    template<typename state_or_mps_type>
    double energy_minus_energy_reduced     (const state_or_mps_type & state, const class_model_infinite & model, const class_edges_infinite & edges);
    template<typename state_or_mps_type>
    double energy_mpo                      (const state_or_mps_type & state, const class_model_infinite & model, const class_edges_infinite & edges);
    template<typename state_or_mps_type>
    double energy_per_site_mpo             (const state_or_mps_type & state, const class_model_infinite & model, const class_edges_infinite & edges);
    template<typename state_or_mps_type>
    double energy_variance_mpo             (const state_or_mps_type & state, const class_model_infinite & model, const class_edges_infinite & edges);
    template<typename state_or_mps_type>
    double energy_variance_per_site_mpo    (const state_or_mps_type & state, const class_model_infinite & model, const class_edges_infinite & edges);

    extern double energy_mpo                      (const class_tensors_infinite & tensors);
    extern double energy_per_site_mpo             (const class_tensors_infinite & tensors);
    extern double energy_variance_mpo             (const class_tensors_infinite & tensors);
    extern double energy_variance_per_site_mpo    (const class_tensors_infinite & tensors);

    extern double energy_mpo                      (const Eigen::Tensor<Scalar,3> &mps, const class_tensors_infinite & tensors);
    extern double energy_per_site_mpo             (const Eigen::Tensor<Scalar,3> &mps, const class_tensors_infinite & tensors);
    extern double energy_variance_mpo             (const Eigen::Tensor<Scalar,3> &mps, const class_tensors_infinite & tensors);
    extern double energy_variance_per_site_mpo    (const Eigen::Tensor<Scalar,3> &mps, const class_tensors_infinite & tensors);

    extern double energy_per_site_ham             (const class_tensors_infinite & tensors);
    extern double energy_per_site_mom             (const class_tensors_infinite & tensors);
    extern double energy_variance_per_site_ham    (const class_tensors_infinite & tensors);
    extern double energy_variance_per_site_mom    (const class_tensors_infinite & tensors);

    /* clang-format on */

}