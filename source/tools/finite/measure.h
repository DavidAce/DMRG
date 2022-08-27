#pragma once
#include <complex>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

class StateFinite;
class ModelFinite;
class EdgesFinite;
class TensorsFinite;
class MpoSite;
class AlgorithmStatus;
struct MeasurementsTensorsFinite;
namespace tools::finite::measure {
    using real = double;
    using cplx = std::complex<double>;
    extern void do_all_measurements(const TensorsFinite &tensors);
    extern void do_all_measurements(const StateFinite &state);

    struct LocalObservableOp {
        Eigen::Tensor<cplx, 2> op;
        long                   pos;
        mutable bool           used = false;
    };

    struct LocalObservableMpo {
        Eigen::Tensor<cplx, 4> mpo;
        long                   pos;
        mutable bool           used = false;
    };

    /* clang-format off */


    [[nodiscard]] extern size_t length                                      (const TensorsFinite & tensors);
    [[nodiscard]] extern size_t length                                      (const StateFinite & state);
    [[nodiscard]] extern size_t length                                      (const ModelFinite & model);
    [[nodiscard]] extern size_t length                                      (const EdgesFinite & edges);
    [[nodiscard]] extern long   bond_dimension_current                      (const StateFinite & state);
    [[nodiscard]] extern long   bond_dimension_midchain                     (const StateFinite & state);
    [[nodiscard]] extern std::vector<long> bond_dimensions_merged           (const StateFinite & state);
    [[nodiscard]] extern std::vector<long> bond_dimensions                  (const StateFinite & state);
    [[nodiscard]] extern double norm                                        (const StateFinite & state, bool full = false);

//  [[nodiscard]]  extern double norm_fast                                   (const StateFinite & state);


    [[nodiscard]] extern double spin_component                              (const StateFinite & state, const Eigen::Matrix2cd &paulimatrix);
    [[nodiscard]] extern double spin_component                              (const StateFinite & state, std::string_view axis);
    [[nodiscard]] extern double spin_alignment                              (const StateFinite & state, std::string_view axis);
    [[nodiscard]] extern int    spin_sign                                   (const StateFinite & state, std::string_view axis);
    [[nodiscard]] extern Eigen::Tensor<cplx,1> mps_wavefn                   (const StateFinite & state);
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
    [[nodiscard]] double energy_minus_energy_shift               (const state_or_mps_type & state, const ModelFinite & model, const EdgesFinite & edges, MeasurementsTensorsFinite * measurements = nullptr);
    template<typename state_or_mps_type>
    [[nodiscard]] double energy                                  (const state_or_mps_type & state, const ModelFinite & model, const EdgesFinite & edges, MeasurementsTensorsFinite * measurements = nullptr);
    template<typename state_or_mps_type>
    [[nodiscard]] double energy_per_site                         (const state_or_mps_type & state, const ModelFinite & model, const EdgesFinite & edges, MeasurementsTensorsFinite * measurements = nullptr);
    template<typename state_or_mps_type>
    [[nodiscard]] double energy_variance                         (const state_or_mps_type & state, const ModelFinite & model, const EdgesFinite & edges, MeasurementsTensorsFinite * measurements = nullptr);
    template<typename state_or_mps_type>
    [[nodiscard]] double energy_variance_per_site                (const state_or_mps_type & state, const ModelFinite & model, const EdgesFinite & edges, MeasurementsTensorsFinite * measurements = nullptr);

    template<typename state_or_mps_type>
    [[nodiscard]] double energy_normalized                       (const state_or_mps_type & state, const ModelFinite & model, const EdgesFinite & edges, double energy_min, double energy_max, MeasurementsTensorsFinite * measurements = nullptr);


    [[nodiscard]] extern double energy_shift                    (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_shift_per_site           (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_minus_energy_shift       (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy                          (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_per_site                 (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_variance                 (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_variance_per_site        (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_normalized               (const TensorsFinite & tensors, double energy_minimum, double energy_maximum);

    [[nodiscard]] extern double energy_minus_energy_shift  (const StateFinite & state, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy                     (const StateFinite & state, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy_per_site            (const StateFinite & state, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy_variance            (const StateFinite & state, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy_variance_per_site   (const StateFinite & state, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy_normalized          (const StateFinite & state, const TensorsFinite & tensors, double energy_minimum, double energy_maximum, MeasurementsTensorsFinite * measurements = nullptr);


    [[nodiscard]] extern double energy_minus_energy_shift   (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy                      (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy_per_site             (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy_variance             (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy_variance_per_site    (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy_normalized           (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors, double energy_minimum, double energy_maximum, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double residual_norm                    (const Eigen::Tensor<cplx, 3> &mps,
                                                             const Eigen::Tensor<cplx, 4> &mpo,
                                                             const Eigen::Tensor<cplx, 3> &envL,
                                                             const Eigen::Tensor<cplx, 3> &envR);
    [[nodiscard]] extern double residual                    (const TensorsFinite & tensors);



    [[nodiscard]] extern double                   expectation_value      (const StateFinite & op, const std::vector<LocalObservableOp> & ops);
    [[nodiscard]] extern double                   expectation_value      (const StateFinite & state, const std::vector<LocalObservableMpo> & mpos);
    [[nodiscard]] extern Eigen::Tensor<double, 1> expectation_values     (const StateFinite & state, const Eigen::Tensor<cplx,2> &op);
    [[nodiscard]] extern Eigen::Tensor<double, 1> expectation_values     (const StateFinite & state, const Eigen::Tensor<cplx,4> &mpo);
    [[nodiscard]] extern double                   correlation            (const StateFinite & state, const Eigen::Tensor<cplx,2> &op1, const Eigen::Tensor<cplx,2> &op2, long pos1, long pos2);
    [[nodiscard]] extern Eigen::Tensor<double, 2> correlation_matrix     (const StateFinite & state, const Eigen::Tensor<cplx,2> &op1, const Eigen::Tensor<cplx,2> &op2);
    [[nodiscard]] extern Eigen::Tensor<double, 2> kvornings_matrix       (const StateFinite & state);
                  extern void                     kvornings_marker       (const StateFinite & state);
                  extern void                     expectation_values_xyz (const StateFinite & state);
                  extern void                     correlation_matrix_xyz (const StateFinite & state);
    [[nodiscard]] extern double                   structure_factor       (const StateFinite & state, const Eigen::Tensor<double, 2> &correlation_matrix);
                  extern void                     structure_factors_xyz  (const StateFinite & state);

                  extern void                     parity_components(const StateFinite &state, const Eigen::Matrix2cd &paulimatrix);

}

/* clang-format on */
