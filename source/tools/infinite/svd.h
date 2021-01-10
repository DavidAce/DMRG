#pragma once
#include <complex>
#include <general/eigen_tensor_fwd_decl.h>

class class_state_infinite;

namespace tools::infinite::svd {
    using Scalar = std::complex<double>;
    extern void truncate_theta(Eigen::Tensor<Scalar, 3> &mps, class_state_infinite &state);
}
