#pragma once
#include <string>

class class_state_infinite;
namespace tools::infinite::mps {
    using Scalar = std::complex<double>;
    extern class_state_infinite set_random_state(const class_state_infinite &state, [[maybe_unused]] const std::string &parity);
//    extern void insert_twosite_tensor (class_state_infinite & state, const Eigen::Tensor<Scalar,3> & twosite_tensor, std::optional<long> chi_lim = std::nullopt, std::optional<double> svd_threshold = std::nullopt);
}