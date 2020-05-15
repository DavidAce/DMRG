#pragma once
#include <complex>
#include <unsupported/Eigen/CXX11/Tensor>

class class_state_infinite;

namespace tools::infinite::svd{
    using Scalar = std::complex<double>;
    extern void                     truncate_theta(Eigen::Tensor<Scalar, 3> &mps, class_state_infinite &state);
}
