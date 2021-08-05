#pragma once
#include <config/enums.h>
#include <string>
#include <unsupported/Eigen/CXX11/Tensor>

class TensorsInfinite;
class StateInfinite;
namespace tools::infinite::opt {
    using cplx = std::complex<double>;
    extern Eigen::Tensor<cplx, 3> find_ground_state(const TensorsInfinite &tensors, std::string_view ritz);
    extern Eigen::Tensor<cplx, 3> find_ground_state(const TensorsInfinite &tensors, OptRitz ritz);
    extern Eigen::Tensor<cplx, 3> time_evolve_state(const StateInfinite &state, const Eigen::Tensor<cplx, 2> &U);
}
