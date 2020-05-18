#pragma once
#include <complex>
#include <optional>
#include <list>
#include <unsupported/Eigen/CXX11/Tensor>

class class_state_finite;
class class_model_finite;

namespace tools::finite::views{
    extern const Eigen::Tensor<std::complex<double>,3> &get_multisite_tensor(const class_state_finite & state, std::optional<std::list<size_t>> active_sites = std::nullopt);
    extern const Eigen::Tensor<std::complex<double>,4> &get_multisite_tensor(const class_model_finite & model, std::optional<std::list<size_t>> active_sites = std::nullopt);

}