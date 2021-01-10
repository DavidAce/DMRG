#pragma once
#include <Eigen/Core>
#include <complex>
#include <config/enums.h>
#include <general/eigen_tensor_fwd_decl.h>
#include <optional>
#include <set>
#include <string>

class class_state_infinite;
namespace tools::infinite::mps {
    using Scalar = std::complex<double>;
    extern void merge_twosite_tensor(class_state_infinite &state, const Eigen::Tensor<Scalar, 3> &twosite_tensor, long chi_lim,
                                     std::optional<double> svd_threshold = std::nullopt);
    extern void random_product_state(const class_state_infinite &state, [[maybe_unused]] const std::string &sector, [[maybe_unused]] long bitfield,
                                     bool use_eigenspinors);
}