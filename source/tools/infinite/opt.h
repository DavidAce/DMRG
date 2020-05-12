#pragma once
#include <string>
#include <general/nmspc_tensor_omp.h>
#include <simulation/enums.h>
class class_state_infinite;

namespace tools::infinite::opt {
    using Scalar = std::complex<double>;
    extern Eigen::Tensor<Scalar, 4> find_ground_state(const class_state_infinite &state, std::string_view ritz);
    extern Eigen::Tensor<Scalar, 4> find_ground_state(const class_state_infinite &state, StateRitz ritz);
    extern Eigen::Tensor<Scalar, 4> time_evolve_theta(const class_state_infinite &state, const Eigen::Tensor<Scalar, 4> &U);
    extern void                     truncate_theta(Eigen::Tensor<Scalar, 4> &theta, class_state_infinite &state);

}
