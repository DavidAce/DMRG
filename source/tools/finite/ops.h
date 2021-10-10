#pragma once
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>
class StateFinite;
class MpoSite;

namespace tools::finite::ops {
    using Scalar = std::complex<double>;

    /* clang-format off */
    extern void apply_mpo                                   (StateFinite & state, const Eigen::Tensor<Scalar,4> & mpo, const Eigen::Tensor<Scalar,3> &Ledge, const Eigen::Tensor<Scalar,3> & Redge);
    extern void apply_mpos                                  (StateFinite & state, const std::vector<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,1> & Ledge, const Eigen::Tensor<Scalar,1> & Redge);
    extern void apply_mpos                                  (StateFinite & state, const std::vector<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
    extern std::optional<double> get_spin_component_in_sector(StateFinite & state, std::string_view sector);
    extern void project_to_nearest_sector                   (StateFinite & state, std::string_view sector, std::optional<long> chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void project_to_sector                           (StateFinite & state, const Eigen::MatrixXcd & paulimatrix, int sign, std::optional<long> chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    [[nodiscard]] extern
    StateFinite get_projection_to_sector                    (const StateFinite & state, const Eigen::MatrixXcd & paulimatrix, int sign, std::optional<long> chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    [[nodiscard]] extern
    StateFinite get_projection_to_nearest_sector            (const StateFinite & state, std::string_view  sector, std::optional<long> chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);


    [[nodiscard]] extern double overlap                     (const StateFinite & state1, const StateFinite & state2);
//    [[nodiscard]] extern double expectation_value           (const StateFinite & state1, const StateFinite & state2, const std::vector<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
//    [[nodiscard]] extern double exp_sq_value                (const StateFinite & state1, const StateFinite & state2, const std::vector<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,4> & Ledge, const Eigen::Tensor<Scalar,4> & Redge);
    /* clang-format on */
}