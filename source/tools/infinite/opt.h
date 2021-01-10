#pragma once
#include <config/enums.h>
#include <string>
#include <unsupported/Eigen/CXX11/Tensor>

class class_tensors_infinite;
namespace tools::infinite::opt {
    using Scalar = std::complex<double>;
    extern Eigen::Tensor<Scalar, 3> find_ground_state(const class_tensors_infinite &tensors, std::string_view ritz);
    extern Eigen::Tensor<Scalar, 3> find_ground_state(const class_tensors_infinite &tensors, StateRitz ritz);
    extern Eigen::Tensor<Scalar, 3> time_evolve_state(const class_state_infinite &state, const Eigen::Tensor<Scalar, 2> &U);
}
