#pragma once
#include <complex>
#include <general/eigen_tensor_fwd_decl.h>
#include <optional>
#include <vector>
class class_state_finite;
class class_model_finite;

namespace tools::finite::views {
    extern const Eigen::Tensor<std::complex<double>, 3> &get_multisite_tensor(const class_state_finite &         state,
                                                                              std::optional<std::vector<size_t>> active_sites = std::nullopt);
    extern const Eigen::Tensor<std::complex<double>, 4> &get_multisite_tensor(const class_model_finite &         model,
                                                                              std::optional<std::vector<size_t>> active_sites = std::nullopt);

}