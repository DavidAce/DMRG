#pragma once
#include "math/svd/config.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>
class StateFinite;
class MpoSite;

namespace tools::finite::ops {
    using Scalar = std::complex<double>;

    /* clang-format off */
                  extern void                  apply_mpo                        (StateFinite & state, const Eigen::Tensor<Scalar,4> & mpo, const Eigen::Tensor<Scalar,3> &Ledge, const Eigen::Tensor<Scalar,3> & Redge);
                  extern void                  apply_mpos                       (StateFinite & state, const std::vector<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,1> & Ledge, const Eigen::Tensor<Scalar,1> & Redge);
                  extern void                  apply_mpos                       (StateFinite & state, const std::vector<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
    [[nodiscard]] extern std::optional<double> get_spin_component_along_axis    (StateFinite & state, std::string_view axis);
                  extern void                  project_to_axis                  (StateFinite & state, const Eigen::MatrixXcd & paulimatrix, int sign, std::optional<svd::config> svd_cfg = std::nullopt);
    [[nodiscard]] extern int                   project_to_nearest_axis          (StateFinite & state, std::string_view axis, std::optional<svd::config> svd_cfg = std::nullopt);
    [[nodiscard]] extern StateFinite           get_projection_to_axis           (const StateFinite & state, const Eigen::MatrixXcd & paulimatrix, int sign, std::optional<svd::config> svd_cfg = std::nullopt);
    [[nodiscard]] extern StateFinite           get_projection_to_nearest_axis   (const StateFinite & state, std::string_view  axis, std::optional<svd::config> svd_cfg = std::nullopt);
    [[nodiscard]] extern double                overlap                          (const StateFinite & state1, const StateFinite & state2);
    /* clang-format on */
}