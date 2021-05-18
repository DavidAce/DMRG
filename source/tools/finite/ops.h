#pragma once
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>
class class_state_finite;
class class_mpo_site;

namespace tools::finite::ops {
    using Scalar = std::complex<double>;
    /* clang-format off */
    extern void apply_mpo                                   (class_state_finite & state, const Eigen::Tensor<Scalar,4> & mpo, const Eigen::Tensor<Scalar,3> &Ledge, const Eigen::Tensor<Scalar,3> & Redge);
    extern void apply_mpos                                  (class_state_finite & state, const std::vector<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,1> & Ledge, const Eigen::Tensor<Scalar,1> & Redge);
    extern void apply_mpos                                  (class_state_finite & state, const std::vector<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
    extern void project_to_nearest_sector                   (class_state_finite & state, const std::string &sector);
    extern void project_to_sector                           (class_state_finite & state, const Eigen::MatrixXcd & paulimatrix, int sign);
    [[nodiscard]] extern
    class_state_finite get_projection_to_sector             (const class_state_finite & state, const Eigen::MatrixXcd & paulimatrix, int sign);
    [[nodiscard]] extern
    class_state_finite get_projection_to_nearest_sector     (const class_state_finite & state, const std::string & sector);

    [[nodiscard]] extern double overlap                     (const class_state_finite & state1, const class_state_finite & state2);
    [[nodiscard]] extern double expectation_value           (const class_state_finite & state1, const class_state_finite & state2, const std::vector<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
    [[nodiscard]] extern double exp_sq_value                (const class_state_finite & state1, const class_state_finite & state2, const std::vector<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,4> & Ledge, const Eigen::Tensor<Scalar,4> & Redge);
    /* clang-format on */
}