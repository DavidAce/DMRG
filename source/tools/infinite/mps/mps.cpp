#include "../mps.h"
#include <config/settings.h>
#include <io/fmt.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateInfinite.h>
#include <tools/common/split.h>
using Scalar = std::complex<double>;

void tools::infinite::mps::merge_twosite_tensor(StateInfinite &state, const Eigen::Tensor<Scalar, 3> &twosite_tensor, long chi_lim,
                                                std::optional<svd::settings> svd_settings) {
    long   dA       = state.get_spin_dimA();
    long   dB       = state.get_spin_dimB();
    size_t posA     = state.get_positionA();
    size_t posB     = state.get_positionB();
    auto   mps_list = tools::common::split::split_mps(twosite_tensor, {dA, dB}, {posA, posB}, static_cast<long>(posA), chi_lim, svd_settings);
    if(mps_list.size() != 2) throw std::logic_error(fmt::format("Got {} MPS sites from two-site tensor.", mps_list.size()));
    state.get_mps_siteA().fuse_mps(mps_list.front());
    state.get_mps_siteB().fuse_mps(mps_list.back());
}

void tools::infinite::mps::random_product_state([[maybe_unused]] const StateInfinite &state, [[maybe_unused]] const std::string &sector,
                                                [[maybe_unused]] long bitfield, [[maybe_unused]] bool use_eigenspinors) {
    throw std::runtime_error("random product state for infinite state not implemented yet");
}

// void tools::infinite::mps::merge_multisite_tensor(StateInfinite &state, const Eigen::Tensor<Scalar, 3> &multisite_mps, std::optional<long> chi_lim,
// std::optional<double> svd_threshold) {
//    // Some sanity checks
//    state.get_mps()
//    if(multisite_mps.dimension(1) != state.get_mps(positions.front()).get_chiL())
//        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim1 {} != chiL on left-most site {}", multisite_mps.dimension(1),
//                                             state.get_mps(positions.front()).get_chiL(), positions.front()));
//
//    if(multisite_mps.dimension(2) != state.get_mps(positions.back()).get_chiR())
//        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim2 {} != chiR on right-most site {}",
//        multisite_mps.dimension(2),
//                                             state.get_mps(positions.back()).get_chiR(), positions.back()));
//    if(not chi_lim)
//        chi_lim = state.get_chi_lim();
//
//    std::vector<long> spin_dims;
//    for(const auto &site : positions) spin_dims.emplace_back(state.get_mps(site).spin_dim());
//
//    // Split the multisite mps into single-site mps objects
//    auto mps_list = tools::common::split::split_mps(multisite_mps, spin_dims, positions, center_position, chi_lim.value(), svd_threshold);
//
//    if(positions.size() != mps_list.size())
//        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: number of sites mismatch: positions.size() {} != mps_list.size() {}",
//                                             positions.size(), mps_list.size()));
//
//    // Note that one of the positions on the split will be a center, so we need to unset
//    // the center in our current state so we don't get duplicate centers
//    state.get_mps().unset_LC();
//
//    // Copy the split up mps components into the current state
//    auto mps_tgt = std::next(state.MPS.begin(), static_cast<long>(positions.front()));
//    for(const auto &mps_src : mps_list) {
//        if(mps_tgt->get_position() != mps_src.get_position())
//            throw std::runtime_error(fmt::format("Could not merge multisite mps into state: Position mismatch: mps_tgt pos {} != mps_src pos {}",
//                                                 mps_tgt->get_position(), mps_src.get_position()));
//        mps_tgt->set_M(mps_src.get_M());
//        if(mps_src.get_L().size() > 0) // The edges have empty "L"
//            mps_tgt->set_L(mps_src.get_L());
//        if(mps_src.isCenter()) mps_tgt->set_LC(mps_src.get_LC());
//        mps_tgt++;
//    }
//    state.clear_cache();
//    state.clear_measurements();
//}
