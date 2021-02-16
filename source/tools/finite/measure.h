//
// Created by david on 2019-10-13.
//

#pragma once
#include <complex>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

/* clang-format off */
class class_state_finite;
class class_model_finite;
class class_edges_finite;
class class_tensors_finite;
class class_mpo_site;
class class_algorithm_status;
struct tensors_measure_finite;
namespace tools::finite::measure{
    using Scalar = std::complex<double>;
    extern void do_all_measurements(const class_tensors_finite & tensors);
    extern void do_all_measurements(const class_state_finite & state);

    [[nodiscard]] extern size_t length                                      (const class_tensors_finite & tensors);
    [[nodiscard]] extern size_t length                                      (const class_state_finite & state);
    [[nodiscard]] extern size_t length                                      (const class_model_finite & model);
    [[nodiscard]] extern size_t length                                      (const class_edges_finite & edges);
    [[nodiscard]] extern long   bond_dimension_current                      (const class_state_finite & state);
    [[nodiscard]] extern long   bond_dimension_midchain                     (const class_state_finite & state);
    [[nodiscard]] extern std::vector<long> bond_dimensions_merged           (const class_state_finite & state);
    [[nodiscard]] extern std::vector<long> bond_dimensions                  (const class_state_finite & state);
    [[nodiscard]] extern double norm                                        (const class_state_finite & state);
//  [[nodiscard]]  extern double norm_fast                                   (const class_state_finite & state);


    [[nodiscard]] extern double spin_component                              (const class_state_finite & state, const Eigen::Matrix2cd &paulimatrix);
    [[nodiscard]] extern double spin_component                              (const class_state_finite & state, const std::string & axis);
    [[nodiscard]] extern Eigen::Tensor<Scalar,1> mps_wavefn                 (const class_state_finite & state);
    [[nodiscard]] extern double entanglement_entropy_current                (const class_state_finite & state);
    [[nodiscard]] extern double entanglement_entropy_midchain               (const class_state_finite & state);
    [[nodiscard]] extern std::vector<double> entanglement_entropies         (const class_state_finite & state);
    [[nodiscard]] extern std::vector<double> renyi_entropies                (const class_state_finite & state, double q);
    [[nodiscard]] extern double number_entropy_current                      (const class_state_finite & state);
    [[nodiscard]] extern double number_entropy_midchain                     (const class_state_finite & state);
    [[nodiscard]] extern std::vector<double> number_entropies               (const class_state_finite & state);
    [[nodiscard]] extern std::array<double,3> spin_components               (const class_state_finite & state);
    [[nodiscard]] extern std::vector<double> truncation_errors              (const class_state_finite & state);
    [[nodiscard]] extern std::vector<double> truncation_errors_active       (const class_state_finite & state);


    template<typename state_or_mps_type>
    [[nodiscard]] double energy_minus_energy_reduced             (const state_or_mps_type & state, const class_model_finite & model, const class_edges_finite & edges, tensors_measure_finite * measurements = nullptr);
    template<typename state_or_mps_type>
    [[nodiscard]] double energy                                  (const state_or_mps_type & state, const class_model_finite & model, const class_edges_finite & edges, tensors_measure_finite * measurements = nullptr);
    template<typename state_or_mps_type>
    [[nodiscard]] double energy_per_site                         (const state_or_mps_type & state, const class_model_finite & model, const class_edges_finite & edges, tensors_measure_finite * measurements = nullptr);
    template<typename state_or_mps_type>
    [[nodiscard]] double energy_variance                         (const state_or_mps_type & state, const class_model_finite & model, const class_edges_finite & edges, tensors_measure_finite * measurements = nullptr);
    template<typename state_or_mps_type>
    [[nodiscard]] double energy_variance_per_site                (const state_or_mps_type & state, const class_model_finite & model, const class_edges_finite & edges, tensors_measure_finite * measurements = nullptr);

    template<typename state_or_mps_type>
    [[nodiscard]] double energy_normalized                       (const state_or_mps_type & state, const class_model_finite & model, const class_edges_finite & edges, double energy_minimum, double energy_maximum, tensors_measure_finite * measurements = nullptr);


    [[nodiscard]] extern double energy_reduced(const class_tensors_finite & tensors);
    [[nodiscard]] extern double energy_per_site_reduced(const class_tensors_finite & tensors);
    [[nodiscard]] extern double energy_minus_energy_reduced(const class_tensors_finite & tensors);
    [[nodiscard]] extern double energy(const class_tensors_finite & tensors);
    [[nodiscard]] extern double energy_per_site(const class_tensors_finite & tensors);
    [[nodiscard]] extern double energy_variance(const class_tensors_finite & tensors);
    [[nodiscard]] extern double energy_variance_per_site(const class_tensors_finite & tensors);
    [[nodiscard]] extern double energy_normalized(const class_tensors_finite & tensors, double energy_minimum, double energy_maximum);

    [[nodiscard]] extern double energy_minus_energy_reduced(const Eigen::Tensor<Scalar,3> &mps, const class_tensors_finite & tensors);
    [[nodiscard]] extern double energy(const Eigen::Tensor<Scalar,3> &mps, const class_tensors_finite & tensors);
    [[nodiscard]] extern double energy_per_site(const Eigen::Tensor<Scalar,3> &mps, const class_tensors_finite & tensors);
    [[nodiscard]] extern double energy_variance(const Eigen::Tensor<Scalar,3> &mps, const class_tensors_finite & tensors);
    [[nodiscard]] extern double energy_variance_per_site(const Eigen::Tensor<Scalar,3> &mps, const class_tensors_finite & tensors);
    [[nodiscard]] extern double energy_normalized(const Eigen::Tensor<Scalar,3> &mps, const class_tensors_finite & tensors, double energy_minimum, double energy_maximum);
}

/* clang-format on */
