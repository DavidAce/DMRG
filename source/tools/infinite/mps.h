#pragma once
#include <complex>
#include <config/enums.h>
#include <math/svd/settings.h>
#include <math/tenx/fwd_decl.h>
#include <optional>
#include <set>
#include <string>

class StateInfinite;
namespace tools::infinite::mps {
    using Scalar = std::complex<double>;
    extern void merge_twosite_tensor(StateInfinite &state, const Eigen::Tensor<Scalar, 3> &twosite_tensor, long bond_limit,
                                     std::optional<svd::settings> svd_settings = std::nullopt);
    extern void random_product_state(const StateInfinite &state, [[maybe_unused]] std::string_view sector, [[maybe_unused]] long bitfield,
                                     bool use_eigenspinors);
}