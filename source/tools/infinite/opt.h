#pragma once
#include "config/enums.h"
#include <string>
#include <unsupported/Eigen/CXX11/Tensor>

class TensorsInfinite;
class StateInfinite;
namespace tools::infinite::opt {
    using cx64 = std::complex<double>;
    extern Eigen::Tensor<cx64, 3> find_ground_state(const TensorsInfinite &tensors, std::string_view ritz);
    extern Eigen::Tensor<cx64, 3> find_ground_state(const TensorsInfinite &tensors, OptRitz ritz);
    extern Eigen::Tensor<cx64, 3> time_evolve_state(const StateInfinite &state, const Eigen::Tensor<cx64, 2> &U);
}
