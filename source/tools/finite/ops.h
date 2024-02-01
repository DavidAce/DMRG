#pragma once
#include "math/svd/config.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>
class StateFinite;
class MpoSite;

namespace tools::finite::ops {
    using cplx = std::complex<double>;

    /* clang-format off */
                  extern void                  apply_mpo                        (StateFinite & state, const Eigen::Tensor<cplx,4> & mpo, const Eigen::Tensor<cplx,3> &Ledge, const Eigen::Tensor<cplx,3> & Redge);
                  extern void                  apply_mpos                       (StateFinite & state, const std::vector<Eigen::Tensor<cplx,4>> & mpos, const Eigen::Tensor<cplx,1> & Ledge, const Eigen::Tensor<cplx,1> & Redge, bool adjoint = false);
                  extern void                  apply_mpos                       (StateFinite & state, const std::vector<Eigen::Tensor<cplx,4>> & mpos, const Eigen::Tensor<cplx,3> & Ledge, const Eigen::Tensor<cplx,3> & Redge, bool adjoint = false);
                  extern void                  apply_mpos_general               (StateFinite & state, const std::vector<Eigen::Tensor<cplx,4>> & mpos, const Eigen::Tensor<cplx,1> & Ledge, const Eigen::Tensor<cplx,1> & Redge, const svd::config & svd_cfg);
                  extern void                  apply_mpos_general               (StateFinite & state, const std::vector<Eigen::Tensor<cplx,4>> & mpos, const svd::config & svd_cfg);
    [[nodiscard]] extern std::optional<double> get_spin_component_along_axis    (StateFinite & state, std::string_view axis);
                  extern void                  project_to_axis                  (StateFinite & state, const Eigen::MatrixXcd & paulimatrix, int sign, std::optional<svd::config> svd_cfg = std::nullopt);
    [[nodiscard]] extern int                   project_to_nearest_axis          (StateFinite & state, std::string_view axis, std::optional<svd::config> svd_cfg = std::nullopt);
    [[nodiscard]] extern StateFinite           get_projection_to_axis           (const StateFinite & state, const Eigen::MatrixXcd & paulimatrix, int sign, std::optional<svd::config> svd_cfg = std::nullopt);
    [[nodiscard]] extern StateFinite           get_projection_to_nearest_axis   (const StateFinite & state, std::string_view  axis, std::optional<svd::config> svd_cfg = std::nullopt);
    [[nodiscard]] extern cplx                  overlap                          (const StateFinite & state1, const StateFinite & state2);
    /* clang-format on */
}