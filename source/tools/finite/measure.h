//
// Created by david on 2019-10-13.
//

#pragma once
#include <complex>
#include <general/nmspc_tensor_omp.h>
#include <vector>

/* clang-format off */
class class_state_finite;
class class_model_finite;
class class_edges_finite;
class class_tensors_finite;
class class_mpo_site;
class class_algorithm_status;
namespace tools::finite::measure{
    using Scalar = std::complex<double>;
    extern void do_all_measurements(const class_tensors_finite & tensors);
    extern void do_all_measurements(const class_state_finite & state);

    extern size_t length                                      (const class_tensors_finite & tensors);
    extern size_t length                                      (const class_state_finite & state);
    extern size_t length                                      (const class_model_finite & model);
    extern size_t length                                      (const class_edges_finite & edges);
    extern long   bond_dimension_current                      (const class_state_finite & state);
    extern long   bond_dimension_midchain                     (const class_state_finite & state);
    extern std::vector<long> bond_dimensions                  (const class_state_finite & state);
    extern double norm                                        (const class_state_finite & state);
    extern double norm_fast                                   (const class_state_finite & state);


    extern double spin_component                              (const class_state_finite & state, const Eigen::Matrix2cd &paulimatrix);
    extern double spin_component                              (const class_state_finite & state, const std::string & axis);
    extern Eigen::Tensor<Scalar,1> mps_wavefn                 (const class_state_finite & state);
    extern double entanglement_entropy_current                (const class_state_finite & state);
    extern double entanglement_entropy_midchain               (const class_state_finite & state);
    extern std::vector<double> entanglement_entropies         (const class_state_finite & state);
    extern std::array<double,3> spin_components               (const class_state_finite & state);
    extern std::vector<double> truncation_errors              (const class_state_finite & state);


    template<typename state_or_mps_type>
    extern double energy_minus_energy_reduced             (const state_or_mps_type & state, const class_model_finite & model, const class_edges_finite & edges);
    template<typename state_or_mps_type>
    extern double energy                                  (const state_or_mps_type & state, const class_model_finite & model, const class_edges_finite & edges);
    template<typename state_or_mps_type>
    extern double energy_per_site                         (const state_or_mps_type & state, const class_model_finite & model, const class_edges_finite & edges);
    template<typename state_or_mps_type>
    extern double energy_variance                         (const state_or_mps_type & state, const class_model_finite & model, const class_edges_finite & edges);
    template<typename state_or_mps_type>
    extern double energy_variance_per_site                (const state_or_mps_type & state, const class_model_finite & model, const class_edges_finite & edges);

    template<typename state_or_mps_type>
    extern double energy_normalized                       (const state_or_mps_type & state, const class_model_finite & model, const class_edges_finite & edges, double energy_minimum, double energy_maximum);


    extern double energy_minus_energy_reduced(const class_tensors_finite & tensors);
    extern double energy(const class_tensors_finite & tensors);
    extern double energy_per_site(const class_tensors_finite & tensors);
    extern double energy_variance(const class_tensors_finite & tensors);
    extern double energy_variance_per_site(const class_tensors_finite & tensors);
    extern double energy_normalized(const class_tensors_finite & tensors, double energy_minimum, double energy_maximum);

    extern double energy_minus_energy_reduced(const Eigen::Tensor<Scalar,3> &mps, const class_tensors_finite & tensors);
    extern double energy(const Eigen::Tensor<Scalar,3> &mps, const class_tensors_finite & tensors);
    extern double energy_per_site(const Eigen::Tensor<Scalar,3> &mps, const class_tensors_finite & tensors);
    extern double energy_variance(const Eigen::Tensor<Scalar,3> &mps, const class_tensors_finite & tensors);
    extern double energy_variance_per_site(const Eigen::Tensor<Scalar,3> &mps, const class_tensors_finite & tensors);
    extern double energy_normalized(const Eigen::Tensor<Scalar,3> &mps, const class_tensors_finite & tensors, double energy_minimum, double energy_maximum);
}

/* clang-format on */
