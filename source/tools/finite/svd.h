#pragma once
#include <list>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>
class class_state_finite;
namespace tools::finite::svd {
    /* clang-format off */
    using Scalar = std::complex<double>;

    extern void truncate_theta(const Eigen::Tensor<Scalar,3> & theta, class_state_finite & state, std::optional<size_t> chi_lim = std::nullopt, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_left (const Eigen::Tensor<Scalar,3> & theta, class_state_finite & state, std::optional<size_t> chi_lim = std::nullopt, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_right(const Eigen::Tensor<Scalar,3> & theta, class_state_finite & state, std::optional<size_t> chi_lim = std::nullopt, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_theta(const Eigen::Tensor<Scalar,4> & theta, class_state_finite & state, std::optional<size_t> chi_lim = std::nullopt, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_all_sites                  (class_state_finite & state, const size_t & chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_active_sites               (class_state_finite & state, const size_t & chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_next_sites                 (class_state_finite & state, const size_t & chi_lim, size_t num_sites = 4, std::optional<double> svd_threshold = std::nullopt);

    /* clang-format on */
}
