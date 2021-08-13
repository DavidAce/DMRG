#pragma once
#include <complex>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

/* clang-format off */
class StateFinite;
class ModelFinite;
class EdgesFinite;
class TensorsFinite;
class MpoSite;
class AlgorithmStatus;
struct tensors_measure_finite;
namespace tools::finite::measure{
    using real = double;
    using cplx = std::complex<double>;
    extern void do_all_measurements(const TensorsFinite & tensors);
    extern void do_all_measurements(const StateFinite & state);

    [[nodiscard]] extern size_t length                                      (const TensorsFinite & tensors);
    [[nodiscard]] extern size_t length                                      (const StateFinite & state);
    [[nodiscard]] extern size_t length                                      (const ModelFinite & model);
    [[nodiscard]] extern size_t length                                      (const EdgesFinite & edges);
    [[nodiscard]] extern long   bond_dimension_current                      (const StateFinite & state);
    [[nodiscard]] extern long   bond_dimension_midchain                     (const StateFinite & state);
    [[nodiscard]] extern std::vector<long> bond_dimensions_merged           (const StateFinite & state);
    [[nodiscard]] extern std::vector<long> bond_dimensions                  (const StateFinite & state);
    [[nodiscard]] extern double norm                                        (const StateFinite & state);
//  [[nodiscard]]  extern double norm_fast                                   (const StateFinite & state);


    [[nodiscard]] extern double spin_component                              (const StateFinite & state, const Eigen::Matrix2cd &paulimatrix);
    [[nodiscard]] extern double spin_component                              (const StateFinite & state, std::string_view axis);
    [[nodiscard]] extern Eigen::Tensor<cplx,1> mps_wavefn                 (const StateFinite & state);
    [[nodiscard]] extern double entanglement_entropy_current                (const StateFinite & state);
    [[nodiscard]] extern double entanglement_entropy_midchain               (const StateFinite & state);
    [[nodiscard]] extern std::vector<double> entanglement_entropies         (const StateFinite & state);
    [[nodiscard]] extern std::vector<double> renyi_entropies                (const StateFinite & state, double q);
    [[nodiscard]] extern double number_entropy_current                      (const StateFinite & state);
    [[nodiscard]] extern double number_entropy_midchain                     (const StateFinite & state);
    [[nodiscard]] extern std::vector<double> number_entropies               (const StateFinite & state);
    [[nodiscard]] extern std::array<double,3> spin_components               (const StateFinite & state);
    [[nodiscard]] extern std::vector<double> truncation_errors              (const StateFinite & state);
    [[nodiscard]] extern std::vector<double> truncation_errors_active       (const StateFinite & state);


    template<typename state_or_mps_type>
    [[nodiscard]] double energy_minus_energy_reduced             (const state_or_mps_type & state, const ModelFinite & model, const EdgesFinite & edges, tensors_measure_finite * measurements = nullptr);
    template<typename state_or_mps_type>
    [[nodiscard]] double energy                                  (const state_or_mps_type & state, const ModelFinite & model, const EdgesFinite & edges, tensors_measure_finite * measurements = nullptr);
    template<typename state_or_mps_type>
    [[nodiscard]] double energy_per_site                         (const state_or_mps_type & state, const ModelFinite & model, const EdgesFinite & edges, tensors_measure_finite * measurements = nullptr);
    template<typename state_or_mps_type>
    [[nodiscard]] double energy_variance                         (const state_or_mps_type & state, const ModelFinite & model, const EdgesFinite & edges, tensors_measure_finite * measurements = nullptr);
    template<typename state_or_mps_type>
    [[nodiscard]] double energy_variance_per_site                (const state_or_mps_type & state, const ModelFinite & model, const EdgesFinite & edges, tensors_measure_finite * measurements = nullptr);

    template<typename state_or_mps_type>
    [[nodiscard]] double energy_normalized                       (const state_or_mps_type & state, const ModelFinite & model, const EdgesFinite & edges, double energy_minimum, double energy_maximum, tensors_measure_finite * measurements = nullptr);


    [[nodiscard]] extern double energy_reduced                  (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_per_site_reduced         (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_minus_energy_reduced     (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy                          (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_per_site                 (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_variance                 (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_variance_per_site        (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_normalized               (const TensorsFinite & tensors, double energy_minimum, double energy_maximum);

    [[nodiscard]] extern double energy_minus_energy_reduced(const StateFinite & state, const TensorsFinite & tensors);
    [[nodiscard]] extern double energy                     (const StateFinite & state, const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_per_site            (const StateFinite & state, const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_variance            (const StateFinite & state, const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_variance_per_site   (const StateFinite & state, const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_normalized          (const StateFinite & state, const TensorsFinite & tensors, double energy_minimum, double energy_maximum);


    [[nodiscard]] extern double energy_minus_energy_reduced (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors);
    [[nodiscard]] extern double energy                      (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_per_site             (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_variance             (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_variance_per_site    (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_normalized           (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors, double energy_minimum, double energy_maximum);

    [[nodiscard]] extern double max_grad_norm               (const Eigen::Tensor<cplx,3> & mps, const TensorsFinite & tensors);

    template<typename Scalar>
    [[nodiscard]] extern double max_grad_norm               (const Eigen::Tensor<Scalar, 3> &mps,
                                                             const Eigen::Tensor<Scalar, 4> &mpo1,
                                                             const Eigen::Tensor<Scalar, 3> &en1L,
                                                             const Eigen::Tensor<Scalar, 3> &en1R,
                                                             const Eigen::Tensor<Scalar, 4> &mpo2,
                                                             const Eigen::Tensor<Scalar, 3> &en2L,
                                                             const Eigen::Tensor<Scalar, 3> &en2R);

    template<typename Scalar>
    [[nodiscard]] extern double max_grad_norm               (const Eigen::Tensor<Scalar, 3> &mps,
                                                             const Eigen::Tensor<Scalar, 4> &mpo2,
                                                             const Eigen::Tensor<Scalar, 3> &en2L,
                                                             const Eigen::Tensor<Scalar, 3> &en2R);



}

/* clang-format on */
