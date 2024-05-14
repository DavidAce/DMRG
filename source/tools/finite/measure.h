#pragma once
#include "math/float.h"
#include <complex>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

class StateFinite;
class ModelFinite;
class EdgesFinite;
class TensorsFinite;
class MpoSite;
class MpsSite;
class AlgorithmStatus;
struct MeasurementsTensorsFinite;
template<typename T>
struct env_pair;
class EnvEne;
class EnvVar;
namespace svd {
    struct config;
}

namespace tools::finite::measure {
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
    [[nodiscard]] extern std::vector<long>    bond_dimensions_active        (const StateFinite & state);
    [[nodiscard]] extern std::pair<long,long> bond_dimensions               (const StateFinite & state, long pos);
    [[nodiscard]] extern std::vector<long>    bond_dimensions               (const StateFinite & state);
    [[nodiscard]] extern double norm                                        (const StateFinite & state, bool full = false);

//  [[nodiscard]]  extern double norm_fast                                   (const StateFinite & state);


    [[nodiscard]] extern double spin_component                              (const StateFinite & state, const Eigen::Matrix2cd &paulimatrix);
    [[nodiscard]] extern double spin_component                              (const StateFinite & state, std::string_view axis);
    [[nodiscard]] extern double spin_alignment                              (const StateFinite & state, std::string_view axis);
    [[nodiscard]] extern int    spin_sign                                   (const StateFinite & state, std::string_view axis);
    template<typename Scalar = cplx>
    [[nodiscard]] extern Eigen::Tensor<cplx,1> mps2tensor                           (const std::vector<std::unique_ptr<MpsSite>> & mps_sites, std::string_view name);
    [[nodiscard]] extern Eigen::Tensor<cplx,1> mps2tensor                           (const StateFinite & state);
    [[nodiscard]] extern double entanglement_entropy                                (const Eigen::Tensor<cplx,1> & bond);
    [[nodiscard]] extern double entanglement_entropy_current                        (const StateFinite & state);
    [[nodiscard]] extern double entanglement_entropy_midchain                       (const StateFinite & state);
    [[nodiscard]] extern std::vector<double> entanglement_entropies                 (const StateFinite & state);
    [[nodiscard]] extern std::vector<double> entanglement_entropies_log2            (const StateFinite & state);
    [[nodiscard]] extern Eigen::ArrayXXd subsystem_entanglement_entropies_log2      (const StateFinite & state);
    [[nodiscard]] extern Eigen::ArrayXXd subsystem_entanglement_entropies_swap_log2 (const StateFinite & state, const svd::config &svd_cfg);
    [[nodiscard]] extern double renyi_entropy_midchain                              (const StateFinite & state, double q);
    [[nodiscard]] extern std::vector<double> renyi_entropies                        (const StateFinite & state, double q);
    [[nodiscard]] extern double number_entropy_current                              (const StateFinite & state);
    [[nodiscard]] extern double number_entropy_midchain                             (const StateFinite & state);
    [[nodiscard]] extern std::vector<double> number_entropies                       (const StateFinite & state);
    [[nodiscard]] extern std::array<double,3> spin_components                       (const StateFinite & state);
    [[nodiscard]] extern std::vector<double> truncation_errors                      (const StateFinite & state);
    [[nodiscard]] extern std::vector<double> truncation_errors_active               (const StateFinite & state);

    [[nodiscard]] double energy_minus_energy_shift               (const StateFinite & state, const ModelFinite & model, const EdgesFinite & edges, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] double energy                                  (const StateFinite & state, const ModelFinite & model, const EdgesFinite & edges, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] double energy_variance                         (const StateFinite & state, const ModelFinite & model, const EdgesFinite & edges, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] double energy_normalized                       (const StateFinite & state, const ModelFinite & model, const EdgesFinite & edges, double energy_min, double energy_max, MeasurementsTensorsFinite * measurements = nullptr);

    [[nodiscard]] extern double energy_shift                    (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_minus_energy_shift       (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy                          (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_variance                 (const TensorsFinite & tensors);
    [[nodiscard]] extern double energy_normalized               (const TensorsFinite & tensors, double energy_minimum, double energy_maximum);

    [[nodiscard]] extern double energy_minus_energy_shift       (const StateFinite & state, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy                          (const StateFinite & state, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy_variance                 (const StateFinite & state, const TensorsFinite & tensors, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy_normalized               (const StateFinite & state, const TensorsFinite & tensors, double energy_minimum, double energy_maximum, MeasurementsTensorsFinite * measurements = nullptr);


    [[nodiscard]] double energy_minus_energy_shift               (const Eigen::Tensor<cplx,3> & multisite_mps, const ModelFinite & model, const EdgesFinite & edges, std::optional<svd::config> svd_cfg, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] double energy                                  (const Eigen::Tensor<cplx,3> & multisite_mps, const ModelFinite & model, const EdgesFinite & edges, std::optional<svd::config> svd_cfg, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] double energy_variance                         (const Eigen::Tensor<cplx,3> & multisite_mps, const ModelFinite & model, const EdgesFinite & edges, std::optional<svd::config> svd_cfg, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] double energy_normalized                       (const Eigen::Tensor<cplx,3> & multisite_mps, const ModelFinite & model, const EdgesFinite & edges, double energy_min, double energy_max, std::optional<svd::config> svd_cfg, MeasurementsTensorsFinite * measurements = nullptr);


    [[nodiscard]] extern double energy_minus_energy_shift   (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors, std::optional<svd::config> svd_cfg, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy                      (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors, std::optional<svd::config> svd_cfg, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy_variance             (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors, std::optional<svd::config> svd_cfg, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double energy_normalized           (const Eigen::Tensor<cplx,3> &mps, const TensorsFinite & tensors, double energy_minimum, double energy_maximum, std::optional<svd::config> svd_cfg, MeasurementsTensorsFinite * measurements = nullptr);
    [[nodiscard]] extern double residual_norm               (const Eigen::Tensor<cplx, 3> &mps,
                                                             const Eigen::Tensor<cplx, 4> &mpo,
                                                             const Eigen::Tensor<cplx, 3> &envL,
                                                             const Eigen::Tensor<cplx, 3> &envR);
    [[nodiscard]] extern double residual_norm                (const TensorsFinite & tensors);

    [[nodiscard]] extern double residual_norm_full           (const StateFinite &state, const ModelFinite &model);


    [[nodiscard]] extern cplx                     expectation_value      (const StateFinite & state, const std::vector<LocalObservableOp> & ops);
    [[nodiscard]] extern cplx                     expectation_value      (const StateFinite & state, const std::vector<LocalObservableMpo> & mpos);
    [[nodiscard]] extern cplx                     expectation_value      (const StateFinite & state1, const StateFinite & state2,
                                                                          const std::vector<Eigen::Tensor<cplx,4>> & mpos);
    [[nodiscard]] extern cplx_t                   expectation_value      (const StateFinite & state1, const StateFinite & state2,
                                                                          const std::vector<Eigen::Tensor<cplx_t,4>> & mpos_t);
    [[nodiscard]] extern cplx                     expectation_value      (const StateFinite & state1, const StateFinite & state2,
                                                                          const std::vector<Eigen::Tensor<cplx,4>> & mpos,
                                                                          const Eigen::Tensor<cplx,1> & ledge,
                                                                          const Eigen::Tensor<cplx,1> & redge);
    template<typename EnvType>
    [[nodiscard]] extern cplx                     expectation_value      (const std::vector<std::reference_wrapper<const MpsSite>> & mpsBra,
                                                                          const std::vector<std::reference_wrapper<const MpsSite>> & mpsKet,
                                                                          const std::vector<std::reference_wrapper<const MpoSite>> & mpos,
                                                                          const env_pair<EnvType> & envs);
    template<typename EnvType>
    [[nodiscard]] extern cplx                     expectation_value      (const Eigen::Tensor<cplx, 3> &mpsBra, const Eigen::Tensor<cplx, 3> &mpsKet,
                                                                          const std::vector<std::reference_wrapper<const MpoSite>> &mpos, const env_pair<EnvType> &envs,
                                                                          std::optional<svd::config> svd_cfg = std::nullopt);
    template<typename EnvType>
    [[nodiscard]] extern cplx                     expectation_value      (const Eigen::Tensor<cplx, 3> &multisite_mps,
                                                                          const std::vector<std::reference_wrapper<const MpoSite>> &mpos, const env_pair<EnvType> &envs,
                                                                          std::optional<svd::config> svd_cfg = std::nullopt);
    [[nodiscard]] extern Eigen::Tensor<cplx, 1>   expectation_values     (const StateFinite & state, const Eigen::Tensor<cplx,2> &op);
    [[nodiscard]] extern Eigen::Tensor<cplx, 1>   expectation_values     (const StateFinite & state, const Eigen::Tensor<cplx,4> &mpo);
    [[nodiscard]] extern Eigen::Tensor<cplx, 1>   expectation_values     (const StateFinite & state, const Eigen::Matrix2cd &op);
    [[nodiscard]] extern cplx                     correlation            (const StateFinite & state, const Eigen::Tensor<cplx,2> &op1, const Eigen::Tensor<cplx,2> &op2, long pos1, long pos2);
    [[nodiscard]] extern Eigen::Tensor<cplx, 2>   correlation_matrix     (const StateFinite & state, const Eigen::Tensor<cplx,2> &op1, const Eigen::Tensor<cplx,2> &op2);
    [[nodiscard]] extern Eigen::Tensor<cplx, 2>   opdm                   (const StateFinite & state);
    [[nodiscard]] extern Eigen::Tensor<double, 1> opdm_spectrum          (const StateFinite & state);
    [[nodiscard]] extern std::array<Eigen::Tensor<double, 1>, 3>
                                                  expectation_values_xyz (const StateFinite & state);
    [[nodiscard]] extern std::array<double, 3>    expectation_value_xyz  (const StateFinite & state);
    [[nodiscard]] extern std::array<Eigen::Tensor<double,2>, 3>
                                                  correlation_matrix_xyz (const StateFinite & state);
    [[nodiscard]] extern double                   structure_factor       (const StateFinite & state, const Eigen::Tensor<cplx, 2> &correlation_matrix);
    [[nodiscard]] extern std::array<double, 3>    structure_factor_xyz   (const StateFinite & state);
                  extern void                     parity_components(const StateFinite &state, const Eigen::Matrix2cd &paulimatrix);


}

/* clang-format on */
