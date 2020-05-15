#pragma once
#include <list>
#include <complex>
#include <unsupported/Eigen/CXX11/Tensor>
class class_mps_site;

namespace tools::common::svd{
    using Scalar = std::complex<double>;
    extern std::list<class_mps_site> split_mps (Eigen::Tensor<Scalar,3> & mps, long spin_dim = 2);
}