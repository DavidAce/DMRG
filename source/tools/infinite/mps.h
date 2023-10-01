#pragma once
#include "config/enums.h"
#include "math/svd/config.h"
#include "math/tenx/fwd_decl.h"
#include <complex>
#include <optional>
#include <set>
#include <string>

class StateInfinite;
namespace tools::infinite::mps {
    using Scalar = std::complex<double>;
    extern void merge_twosite_tensor(StateInfinite &state, const Eigen::Tensor<Scalar, 3> &twosite_tensor, std::optional<svd::config> svd_cfg = std::nullopt);
    extern void random_product_state(const StateInfinite &state, [[maybe_unused]] std::string_view sector, bool use_eigenspinors,
                                     [[maybe_unused]] std::string & pattern);
}