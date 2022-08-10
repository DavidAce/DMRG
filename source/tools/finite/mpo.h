#pragma once

#include "config/enums.h"
#include "math/tenx/fwd_decl.h"
#include <complex>
#include <vector>

/* clang-format off */
class ModelFinite;
namespace tools::finite::mpo {
    using cplx = std::complex<double>;
    extern std::pair<Eigen::Tensor<cplx, 4>,Eigen::Tensor<cplx, 4>>
                swap_mpo    (const Eigen::Tensor<cplx, 4> & mpoL, const Eigen::Tensor<cplx, 4> & mpoR);
    extern void swap_sites  (ModelFinite & model, size_t posL, size_t posR, std::vector<size_t> & sites);
}

/* clang-format on */
