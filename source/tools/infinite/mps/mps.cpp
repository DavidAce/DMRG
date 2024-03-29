#include "../mps.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateInfinite.h"
#include "tools/common/split.h"
using Scalar = std::complex<double>;

void tools::infinite::mps::merge_twosite_tensor(StateInfinite &state, const Eigen::Tensor<Scalar, 3> &twosite_tensor, std::optional<svd::config> svd_cfg) {
    long   dA       = state.get_spin_dimA();
    long   dB       = state.get_spin_dimB();
    size_t posA     = state.get_positionA();
    size_t posB     = state.get_positionB();
    auto   mps_list = tools::common::split::split_mps(twosite_tensor, {dA, dB}, {posA, posB}, static_cast<long>(posA), svd_cfg);
    if(mps_list.size() != 2) throw except::logic_error("Got {} MPS sites from two-site tensor.", mps_list.size());
    state.get_mps_siteA().fuse_mps(mps_list.front());
    state.get_mps_siteB().fuse_mps(mps_list.back());
}

void tools::infinite::mps::random_product_state([[maybe_unused]] const StateInfinite &state, [[maybe_unused]] std::string_view sector,
                                                [[maybe_unused]] bool use_eigenspinors, [[maybe_unused]] size_t bitfield) {
    throw except::runtime_error("random product state for infinite state not implemented yet");
}
