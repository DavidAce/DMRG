#pragma once

#include "config/enums.h"
#include "math/tenx/fwd_decl.h"
#include <complex>
#include <vector>
#include <math/float.h>

/* clang-format off */
class ModelFinite;
namespace tools::finite::mpo {
    extern std::pair<Eigen::Tensor<cx64, 4>,Eigen::Tensor<cx64, 4>>
                swap_mpo    (const Eigen::Tensor<cx64, 4> & mpoL, const Eigen::Tensor<cx64, 4> & mpoR);
    extern void swap_sites  (ModelFinite & model, size_t posL, size_t posR, std::vector<size_t> & sites);

    extern std::vector<Eigen::Tensor<cx64,4>> get_mpos_with_edges (const std::vector<Eigen::Tensor<cx64,4>> & mpos, const Eigen::Tensor<cx64,1> & Ledge, const Eigen::Tensor<cx64,1> & Redge);
    extern std::vector<Eigen::Tensor<cx128,4>> get_mpos_with_edges_t (const std::vector<Eigen::Tensor<cx128,4>> & mpos, const Eigen::Tensor<cx128,1> & Ledge, const Eigen::Tensor<cx128,1> & Redge);

    extern std::vector<Eigen::Tensor<cx64,4>> get_compressed_mpos (std::vector<Eigen::Tensor<cx64, 4>> mpos, MpoCompress mpoCompress);
    extern std::vector<Eigen::Tensor<cx64,4>> get_compressed_mpos (const std::vector<Eigen::Tensor<cx64,4>> & mpos, const Eigen::Tensor<cx64,1> & Ledge, const Eigen::Tensor<cx64,1> & Redge, MpoCompress mpoCompress);
    extern std::vector<Eigen::Tensor<cx64,4>> get_svdcompressed_mpos (std::vector<Eigen::Tensor<cx64, 4>> mpos);
    extern std::vector<Eigen::Tensor<cx64,4>> get_svdcompressed_mpos (const std::vector<Eigen::Tensor<cx64,4>> & mpos, const Eigen::Tensor<cx64,1> & Ledge, const Eigen::Tensor<cx64,1> & Redge);
    extern std::vector<Eigen::Tensor<cx64,4>> get_inverted_mpos (const std::vector<Eigen::Tensor<cx64,4>> & mpos);
    // extern std::vector<Eigen::Tensor<cx64,4>> get_inverted_mpos (std::vector<Eigen::Tensor<cx64, 4>> mpos);
    // extern std::vector<Eigen::Tensor<cx64,4>> get_inverted_mpos (const std::vector<Eigen::Tensor<cx64,4>> & mpos, const Eigen::Tensor<cx64,1> & Ledge, const Eigen::Tensor<cx64,1> & Redge);
  extern std::vector<Eigen::Tensor<cx64,4>> get_deparallelized_mpos (std::vector<Eigen::Tensor<cx64, 4>> mpos);
  extern std::vector<Eigen::Tensor<cx64,4>> get_deparallelized_mpos (const std::vector<Eigen::Tensor<cx64,4>> & mpos, const Eigen::Tensor<cx64,1> & Ledge, const Eigen::Tensor<cx64,1> & Redge);

}

/* clang-format on */
