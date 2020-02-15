#pragma once
#include <general/nmspc_tensor_omp.h>
#include <list>
class class_state_finite;
class class_model_base;

namespace tools::finite::ops {
    using Scalar = std::complex<double>;
    extern std::list<Eigen::Tensor<Scalar,4>>
    make_mpo_list                 (const std::list<std::unique_ptr<class_model_base>> & mpos_L, const std::list<std::unique_ptr<class_model_base>> & mpos_R);
    extern void apply_mpo                     (class_state_finite & state, const Eigen::Tensor<Scalar,4> & mpo, const Eigen::Tensor<Scalar,3> &Ledge, const Eigen::Tensor<Scalar,3> & Redge);
    extern void apply_mpos                    (class_state_finite & state, const std::list<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
    extern class_state_finite
    get_projection_to_parity_sector           (const class_state_finite & state, const Eigen::MatrixXcd & paulimatrix, int sign);
    extern class_state_finite
    get_projection_to_closest_parity_sector   (const class_state_finite & state, const Eigen::MatrixXcd & paulimatrix);
    extern class_state_finite
    get_projection_to_closest_parity_sector   (const class_state_finite & state, std::string parity_sector);
    extern double overlap                     (const class_state_finite & state1, const class_state_finite & state2);
    extern double expectation_value           (const class_state_finite & state1, const class_state_finite & state2, const std::list<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
    extern double exp_sq_value                (const class_state_finite & state1, const class_state_finite & state2, const std::list<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,4> & Ledge, const Eigen::Tensor<Scalar,4> & Redge);
}