#include <math/tenx.h>
// -- (textra first)
#include "../mps.h"
#include "config/enums.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "math/linalg/tensor.h"
#include "math/num.h"
#include "math/svd.h"
#include "qm/gate.h"
#include "qm/mpo.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"
#include "tools/common/split.h"
#include "tools/finite/measure.h"
#include "tools/finite/ops.h"

namespace settings {
    inline constexpr bool debug_merge = false;
    inline constexpr bool debug_gates = false;
    inline constexpr bool debug_moves = false;
}

bool tools::finite::mps::init::bitfield_is_valid(std::optional<long> bitfield) {
    return bitfield.has_value() and bitfield.value() > 0 and init::used_bitfields.count(bitfield.value()) == 0;
}
bool tools::finite::mps::init::is_valid_axis(std::string_view sector) {
    return std::find(valid_axis_str.begin(), valid_axis_str.end(), sector) != valid_axis_str.end();
}

size_t tools::finite::mps::move_center_point_single_site(StateFinite &state, long bond_limit, std::optional<svd::settings> svd_settings) {
    auto t_move = tid::tic_scope("move");
    if(state.position_is_outward_edge()) {
        if(state.get_direction() == -1 and state.get_mps_site(0l).get_chiL() != 1)
            throw except::logic_error("chiL at position 0 must have dimension 1, but it has dimension {}. Mps dims {}", state.get_mps_site(0l).get_chiL(),
                                      state.get_mps_site(0l).dimensions());
        if(state.get_direction() == 1 and state.get_mps_site().get_chiR() != 1)
            throw except::logic_error("chiR at position {} must have dimension 1, but it has dimension {}. Mps dims {}", state.get_position(),
                                      state.get_mps_site().get_chiR(), state.get_mps_site().dimensions());
        state.flip_direction(); // Instead of moving out of the chain, just flip the direction and return
        return 0;               // No moves this time, return 0
    } else {
        long pos  = state.get_position<long>();  // If all sites are B's, then this is -1. Otherwise this is the current "A*LC" site
        long posC = pos + state.get_direction(); // This is the site which becomes the new center position
        if(pos < -1 or pos >= state.get_length<long>()) throw std::runtime_error(fmt::format("pos out of bounds: {}", pos));
        if(posC < -1 or posC >= state.get_length<long>()) throw std::runtime_error(fmt::format("posC out of bounds: {}", posC));
        if(state.get_direction() != posC - pos) throw std::logic_error(fmt::format("Expected posC - pos == {}. Got {}", state.get_direction(), posC - pos));

        if constexpr(settings::debug_moves) {
            if(posC > pos) tools::log->info("Moving {} -> {}", pos, posC);
            if(posC < pos) tools::log->info("Moving {} <- {}", posC, pos);
        }

        Eigen::Tensor<cplx, 1> LC(1);
        LC.setConstant(1); // Store the LC bond in a temporary. It will become a regular "L" bond later
        if(pos >= 0) LC = state.get_mps_site(pos).get_LC();

        if(state.get_direction() == 1) {
            auto  posC_ul  = static_cast<size_t>(posC);                                                          // Cast to unsigned
            auto &mpsC     = state.get_mps_site(posC);                                                           // This becomes the new center position AC
            long  bond_new = std::min(bond_limit, mpsC.spin_dim() * std::min(mpsC.get_chiL(), mpsC.get_chiR())); // Bond dimensions growth limit
            // Construct a single-site tensor. This is equivalent to state.get_multisite_mps(...) but avoid normalization checks.
            auto onesite_tensor = tools::common::contraction::contract_bnd_mps_temp(LC, mpsC.get_M());
            tools::finite::mps::merge_multisite_mps(state, onesite_tensor, {posC_ul}, posC, bond_new, svd_settings, LogPolicy::QUIET);
        } else if(state.get_direction() == -1) {
            auto  pos_ul         = static_cast<size_t>(pos);                                                        // Cast to unsigned
            auto &mps            = state.get_mps_site(pos);                                                         // This AC becomes the new B
            long  bond_new       = std::min(bond_limit, mps.spin_dim() * std::min(mps.get_chiL(), mps.get_chiR())); // Bond dimensions growth limit
            auto  onesite_tensor = mps.get_M(); // No need to contract anything this time. Note that we must take a copy! Not a reference (LC is unset later)
            tools::finite::mps::merge_multisite_mps(state, onesite_tensor, {pos_ul}, posC, bond_new, svd_settings, LogPolicy::QUIET);
        }
        state.clear_cache(LogPolicy::QUIET);
        state.clear_measurements(LogPolicy::QUIET);
        return 1; // Moved once, so return 1
    }
}

size_t tools::finite::mps::move_center_point(StateFinite &state, long bond_limit, std::optional<svd::settings> svd_settings) {
    auto t_move = tid::tic_scope("move");
    if(state.position_is_outward_edge(2)) {
        state.flip_direction(); // Instead of moving out of the chain, just flip the direction and return
        return 0;               // No moves this time, return 0
    } else {
        long pos = state.get_position<long>();
        if(pos < -1 or pos >= state.get_length<long>()) throw std::runtime_error(fmt::format("pos out of bounds: {}", pos));

        long  posL    = state.get_direction() == 1 ? pos + 1 : pos - 1;
        long  posR    = state.get_direction() == 1 ? pos + 2 : pos;
        auto  posL_ul = static_cast<size_t>(posL);
        auto  posR_ul = static_cast<size_t>(posR);
        auto &mps     = state.get_mps_site();
        auto &mpsL    = state.get_mps_site(posL); // Becomes the new center position
        auto &mpsR    = state.get_mps_site(posR); // The site to the right of the new center position
        // Store the special LC bond in a temporary. It needs to be put back afterwards
        // Do the same with its truncation error
        Eigen::Tensor<cplx, 1> LC                  = mps.get_LC();
        double                 truncation_error_LC = mps.get_truncation_error_LC();
        auto                   twosite_tensor      = state.get_multisite_mps({posL_ul, posR_ul});
        tools::finite::mps::merge_multisite_mps(state, twosite_tensor, {static_cast<size_t>(posL), static_cast<size_t>(posR)}, static_cast<long>(posL),
                                                bond_limit, svd_settings, LogPolicy::QUIET);
        state.clear_cache(LogPolicy::QUIET);
        state.clear_measurements(LogPolicy::QUIET);

        // Put LC where it belongs.
        // Recall that mpsL, mpsR are on the new position, not the old one!
        if(state.get_direction() == 1)
            mpsL.set_L(LC, truncation_error_LC);
        else
            mpsR.set_L(LC, truncation_error_LC);
        return 1; // Moved once, so return 1
    }
}

size_t tools::finite::mps::move_center_point_to_pos(StateFinite &state, long pos, long bond_limit, std::optional<svd::settings> svd_settings) {
    if(pos != std::clamp<long>(pos, -1l, state.get_length<long>() - 1))
        throw std::logic_error(fmt::format("move_center_point_to_pos: Given pos [{}]. Expected range [-1,{}]", pos, state.get_length<long>() - 1));
    if((state.get_direction() < 0 and pos > state.get_position<long>()) or //
       (state.get_direction() > 0 and pos < state.get_position<long>()))   //
        state.flip_direction();                                            // Turn direction towards new position

    size_t moves = 0;
    while(not state.position_is_at(pos)) moves += move_center_point_single_site(state, bond_limit, svd_settings);
    return moves;
}

size_t tools::finite::mps::move_center_point_to_pos_dir(StateFinite &state, long pos, int dir, long bond_limit, std::optional<svd::settings> svd_settings) {
    if(pos != std::clamp<long>(pos, -1l, state.get_length<long>() - 1))
        throw std::logic_error(fmt::format("move_center_point_to_pos_dir: Given pos [{}]. Expected range [-1,{}]", pos, state.get_length<long>() - 1));
    if(std::abs(dir) != 1) throw std::logic_error("move_center_point_to_pos_dir: dir must be 1 or -1");
    if((state.get_direction() < 0 and pos > state.get_position<long>()) or //
       (state.get_direction() > 0 and pos < state.get_position<long>()))   //
        state.flip_direction();                                            // Turn direction towards new position
    size_t moves = 0;
    while(not state.position_is_at(pos)) moves += move_center_point_single_site(state, bond_limit, svd_settings);
    if(dir != state.get_direction()) state.flip_direction();
    return moves;
}

size_t tools::finite::mps::move_center_point_to_edge(StateFinite &state, long bond_limit, std::optional<svd::settings> svd_settings) {
    size_t moves = 0;
    while(not state.position_is_inward_edge()) moves += move_center_point_single_site(state, bond_limit, svd_settings);
    return moves;
}

size_t tools::finite::mps::move_center_point_to_middle(StateFinite &state, long bond_limit, std::optional<svd::settings> svd_settings) {
    size_t moves = 0;
    while(not state.position_is_the_middle_any_direction()) moves += move_center_point_single_site(state, bond_limit, svd_settings);
    return moves;
}

size_t tools::finite::mps::merge_multisite_mps(StateFinite &state, const Eigen::Tensor<cplx, 3> &multisite_mps, const std::vector<size_t> &sites,
                                               long center_position, long bond_limit, std::optional<svd::settings> svd_settings,
                                               std::optional<LogPolicy> logPolicy) {
    auto t_merge          = tid::tic_scope("merge");
    auto current_position = state.get_position<long>();
    auto moves            = static_cast<size_t>(std::abs(center_position - current_position));
    if constexpr(settings::debug_merge or settings::debug)
        if(logPolicy == LogPolicy::NORMAL)
            tools::log->trace("merge_multisite_mps: sites {} | chi limit {} | dimensions {} | center {} -> {} | {}", sites, bond_limit,
                              multisite_mps.dimensions(), current_position, center_position, state.get_labels());

    // Some sanity checks
    if(multisite_mps.dimension(1) != state.get_mps_site(sites.front()).get_chiL())
        throw std::runtime_error(fmt::format("merge_multisite_mps: mps dim1 {} != chiL {} on left-most site", multisite_mps.dimension(1),
                                             state.get_mps_site(sites.front()).get_chiL(), sites.front()));

    if(multisite_mps.dimension(2) != state.get_mps_site(sites.back()).get_chiR())
        throw std::runtime_error(fmt::format("merge_multisite_mps: mps dim2 {} != chiR {} on right-most site", multisite_mps.dimension(2),
                                             state.get_mps_site(sites.back()).get_chiR(), sites.back()));
    if constexpr(settings::debug_merge or settings::debug) {
        // Never allow nan's in the multisite_mps
        if(tenx::hasNaN(multisite_mps))
            throw except::runtime_error("merge_multisite_mps: multisite_mps has nan's:\n"
                                        "sites            :{}\n"
                                        "center_position  :{}\n"
                                        "current_position :{}\n"
                                        "bond_limit          :{}\n"
                                        "multisite_mps    :\n{}",
                                        sites, center_position, current_position, bond_limit, linalg::tensor::to_string(multisite_mps, 3, 6));

        if(state.has_nan())
            throw except::runtime_error("merge_multisite_mps: state has nan's:\n"
                                        "sites            :{}\n"
                                        "center_position  :{}\n"
                                        "current_position :{}\n"
                                        "bond_limit          :{}\n"
                                        "multisite_mps    :\n{}",
                                        sites, center_position, current_position, bond_limit, linalg::tensor::to_string(multisite_mps, 3, 6));

        // We have to allow non-normalized multisite mps! Otherwise we won't be able to make them normalized
        auto norm = tenx::VectorCast(multisite_mps).norm();
        if(std::abs(norm - 1) > 1e-8) tools::log->debug("Multisite mps for positions {} has norm far from unity: {:.16f}", sites, norm);
    }

    // Can't set center on one of sites if the current center is too far away: we would end up with interleaved A's and B sites
    bool center_in_sites = center_position == std::clamp<long>(center_position, static_cast<long>(sites.front()), static_cast<long>(sites.back()));
    bool center_in_range = current_position == std::clamp<long>(current_position, static_cast<long>(sites.front() - 1), static_cast<long>(sites.back()));
    if(center_in_sites and not center_in_range)
        throw std::runtime_error(fmt::format("merge_multisite_mps: cannot merge multisite_mps {} with new center at {}: current center {} is too far", sites,
                                             center_position, current_position));

    long              spin_prod = 1;
    std::vector<long> spin_dims;
    spin_dims.reserve(sites.size());
    for(const auto &pos : sites) {
        spin_dims.emplace_back(state.get_mps_site(pos).spin_dim());
        spin_prod *= spin_dims.back();
    }
    if(spin_prod != multisite_mps.dimension(0))
        throw std::runtime_error(fmt::format("merge_multisite_mps: multisite mps dim0 {} != spin_prod {}", multisite_mps.dimension(0), spin_prod));

    // Hold L on the edge in case we need to convert AL to A or LB to B
    //    std::optional<stash<Eigen::Tensor<cplx, 1>>> ll_edge = std::nullopt;
    //    std::optional<stash<Eigen::Tensor<cplx, 1>>> lr_edge = std::nullopt;
    //    if(num::cmp_greater_equal(center_position, sites.back())) {
    //        // In this case all sites will become A or AC type.
    //        // Then we need to store the last L (the one that was included when constructing the multisite theta).
    //        // This is the L to the right of the last mps in "sites"
    //        auto        pos_back = sites.back();
    //        auto        pos_next = sites.back() + 1;
    //        const auto &mps_back = state.get_mps_site(pos_back);
    //        if(mps_back.isCenter())
    //            lr_edge = {mps_back.get_LC(), mps_back.get_truncation_error_LC(), pos_back};
    //        else if(mps_back.get_label() == "B")
    //            lr_edge = {mps_back.get_L(), mps_back.get_truncation_error(), pos_back};
    //        else if(pos_next < state.get_length()) {
    //            // In this case we know mps_back == A and that the next site can be
    //            // either A or AC, but should not be a B (that would be an error).
    //            const auto &mps_next = state.get_mps_site(pos_next);
    //            if(mps_next.get_label() != "B")
    //                lr_edge = {mps_next.get_L(), mps_next.get_truncation_error(), pos_back};
    //            else
    //                throw except::logic_error("merge_multisite_mps: wrong label [{}] on mps_next | expected A or AC", mps_next.get_label());
    //        }
    //    }
    //    if(num::cmp_less(center_position, sites.back())) {
    //        // In this case all sites will become B type
    //        // Then we need to store the L on the left of all sites in multisite_mps.
    //        auto        pos_frnt = sites.front();
    //        auto        pos_prev = sites.front() - 1;
    //        const auto &mps_frnt = state.get_mps_site(pos_frnt);
    //        if(mps_frnt.get_label() != "B")
    //            ll_edge = {mps_frnt.get_L(), mps_frnt.get_truncation_error(), pos_frnt}; // This one becomes LC on left-move if mps_frnt is AC now
    //        else if(pos_prev > 0) {
    //            const auto &mps_prev = state.get_mps_site(pos_prev);
    //            // In this case we know that mps_frnt is B type, therefore mps_prev must be AC or B
    //            if(mps_prev.isCenter())
    //                ll_edge = {mps_prev.get_LC(), mps_prev.get_truncation_error_LC(), pos_frnt};
    //            else if(mps_prev.get_label() == "B")
    //                ll_edge = {mps_prev.get_L(), mps_prev.get_truncation_error(), pos_frnt};
    //            else
    //                throw except::logic_error("merge_multisite_mps: wrong label [{}] on mps_prev | expected AC or B", mps_prev.get_label());
    //        }
    //    }

    // Hold LC if moving. This should be placed in an L-slot later
    std::optional<stash<Eigen::Tensor<cplx, 1>>> lc_move = std::nullopt;
    if(center_position != current_position and current_position >= 0) {
        auto &mps      = state.get_mps_site(current_position); // Guaranteed to have LC since that is the definition of current_position
        auto  pos_back = static_cast<long>(sites.back());
        auto  pos_frnt = static_cast<long>(sites.front());
        auto  pos_curr = static_cast<size_t>(current_position);

        // Detect right-move
        if(center_position > current_position) { // This AC will become an A (AC moves to the right)
            if(center_position != std::clamp(center_position, pos_frnt, pos_back))
                throw std::logic_error(fmt::format("merge_multisite_mps: right-moving new center position {} must be in sites {}", center_position, sites));

            // Case 1, right-move: LC[3]B[4] -> L[4]A[4]LC[4]V[5], current_position == 3, center_position == 4. Then LC[3] becomes L[4] on A[4]
            // Case 2, swap-move: A[3]LC[3]B[4] -> A[3]A[4]LC[4]V, current_position == 3, center_position == 4. Then LC[3] is thrown away
            // Case 4, deep-move: A[3]A[4]LC[4]B[5]B[6]B[7] -> A[3]A[4]A[5]A[6]LC[6]B[7], current_position == 5, center_position == 6. Then LC[4] is thrown
            // Takeaway: LC is only held when LC is on the left edge, turning a B into an A which needs an L
            // It's important that the last V is a diagonal matrix, otherwise it would truncate the site to the right.
            if(current_position + 1 == pos_frnt) lc_move = stash<Eigen::Tensor<cplx, 1>>{mps.get_LC(), mps.get_truncation_error_LC(), sites.front()};
        }
        // Detect left-move
        if(center_position < current_position) { // This AC position will become a B (AC moves to the left)
            if(center_position < pos_frnt - 1)
                throw std::logic_error(
                    fmt::format("merge_multisite_mps: left-moving new center position {} is out of range [{}]+{}", center_position, pos_frnt - 1, sites));
            if(current_position > pos_back + 1)
                throw std::logic_error(
                    fmt::format("merge_multisite_mps: left-moving current position {} is out of range {}+[{}]", current_position, sites, pos_back + 1));

            // Case 1, left-move: A[3]LC[3]     -> U[2]LC[2]B[3]    , current_position == 3, center_position == 2. Then LC[3] becomes L[3] on B[3]
            // Case 2, swap-move: A[3]A[4]LC[4] -> A[3]LC[3]B[4]    , current_position == 3, center_position == 4. Then LC[4] becomes L[4] on B[4]
            // Case 3, full-move: A[3]A[4]LC[4] -> U[2]LC[2]B[3]B[4], current_position == 3, center_position == 4. Then LC[4] becomes L[4] on B[4]
            // Case 4, deep-move: A[3]A[4]LC[4]B[5]B[6]B[7] -> A[3]LC[4]B[4]B[5]B[6]B[7], current_position == 4 center_position == 3. Then LC[4] is thrown
            // Takeaway: LC is only held when LC is on the right edge, turning an AC into a B which needs an L.
            // It's important that the front U is a diagonal matrix, otherwise it would truncate the site to the left.

            if(current_position == pos_back) lc_move = stash<Eigen::Tensor<cplx, 1>>{mps.get_LC(), mps.get_truncation_error_LC(), pos_curr};
        }
        // Note that one of the positions on the split may contain a new center, so we need to unset
        // the center in our current state so we don't get duplicate centers
        mps.unset_LC();
    }

    if constexpr(settings::debug_merge)
        if(svd_settings) tools::log->trace("merge_multisite_mps: splitting sites {} | {}", sites, svd_settings->to_string());

    // Split the multisite mps into single-site mps objects
    auto mps_list = tools::common::split::split_mps(multisite_mps, spin_dims, sites, center_position, bond_limit, svd_settings);

    // Sanity checks
    if(sites.size() != mps_list.size())
        throw std::runtime_error(
            fmt::format("merge_multisite_mps: number of sites mismatch: sites.size() {} != mps_list.size() {}", sites.size(), mps_list.size()));

    // Fuse the split-up mps components into the current state
    for(auto &mps_src : mps_list) {
        auto  pos     = mps_src.get_position();
        auto &mps_tgt = state.get_mps_site(pos);
        //        tools::log->info("merge_multisite_mps: fusing src {} --> tgt {} | ll {} | lr {} | lc {}",
        //                         mps_src.get_tag(), mps_tgt.get_tag(),
        //                         ll_edge and ll_edge->pos_dst == pos,
        //                         lr_edge and lr_edge->pos_dst == pos,
        //                         lc_move and lc_move->pos_dst == pos
        //                         );
        //        // Now, for split v2, multiply L⁻¹ to generate either A or B from LGL.
        //        if(mps_src.get_label() == "AL") {
        //            tools::log->info("merge_multisite_mps: converting AL to A at pos {} | center {} | sites {}", pos, center_position, sites);
        //            if(lr_edge and lr_edge->pos_dst == pos) {
        //                mps_src.convert_AL_to_A(lr_edge->data);
        //            } else {
        //                Eigen::Tensor<cplx, 1> LR(1);
        //                LR.setConstant(1.0);
        //                mps_src.convert_AL_to_A(LR);
        //            }
        //            // Now it's possible that mps_src is also supposed to become an AC. Then lr_edge has LC.
        //            if(num::cmp_equal(pos, center_position)) mps_src.set_LC(lr_edge->data, lr_edge->error);
        //        }
        //        if(mps_src.get_label() == "LB") {
        //            tools::log->info("merge_multisite_mps: converting LB to B at pos {} | center {} | sites {}", pos, center_position, sites);
        //            if(ll_edge and ll_edge->pos_dst == pos) {
        //                mps_src.convert_LB_to_B(ll_edge->data);
        //            } else {
        //                Eigen::Tensor<cplx, 1> LL(1);
        //                LL.setConstant(1.0);
        //                mps_src.convert_LB_to_B(LL);
        //            }
        //            // Now it's possible that the site left of mps_src is also supposed to become an AC. Then ll_edge has LC.
        //            if(num::cmp_equal(pos - 1, center_position) and ll_edge->pos_dst == pos) {
        //                auto &mps_left = state.get_mps_site(pos - 1);
        //                mps_left.set_LC(ll_edge->data, ll_edge->error);
        //            }
        //        }
        // inject lc_move if there is any waiting
        if(lc_move and pos == lc_move->pos_dst) { mps_src.set_L(lc_move->data, lc_move->error); }

        mps_tgt.fuse_mps(mps_src);
        state.tag_site_normalized(pos, true); // Fused site is normalized

        // Now take stashes for neighboring sites
        if(pos < state.get_length() - 1) state.get_mps_site(pos + 1).take_stash(mps_src); // Take stashed S,V (and possibly LC)
        if(pos > 0) state.get_mps_site(pos - 1).take_stash(mps_src);                      // Take stashed U,S (and possibly LC)
        mps_src.drop_stash(); // Discard whatever is left stashed at the edge (this normalizes the state)
    }

    current_position = state.get_position<long>();
    if(current_position != center_position)
        throw std::logic_error(fmt::format("Center position mismatch {} ! {}\nLabels: {}", current_position, center_position, state.get_labels()));
    state.clear_cache(LogPolicy::QUIET);
    state.clear_measurements(LogPolicy::QUIET);

    if constexpr(settings::debug) {
        auto t_dbg = tid::tic_scope("debug");
        state.assert_validity();
        for(auto &pos : sites) state.get_mps_site(pos).assert_identity();
    }
    return moves;
}

bool tools::finite::mps::normalize_state(StateFinite &state, std::optional<long> bond_limit, std::optional<svd::settings> svd_settings,
                                         NormPolicy norm_policy) {
    // When a state needs to be normalized it's enough to "move" the center position around the whole chain.
    // Each move performs an SVD decomposition which leaves unitaries behind, effectively normalizing the state.
    // NOTE! It IS important to start with the current position.

    if(norm_policy == NormPolicy::IFNEEDED) {
        // We may only go ahead with a normalization if its really needed.
        auto norm = tools::finite::measure::norm(state);
        tools::log->trace("Norm: {:.16f}", norm);
        if(std::abs(norm - 1.0) < settings::precision::max_norm_error) return false;
    }
    // Otherwise we just do the normalization
    if(not bond_limit) bond_limit = state.find_largest_bond();
    // Save the current position, direction and center status
    auto dir   = state.get_direction();
    auto pos   = state.get_position<long>();
    auto cnt   = pos >= 0;
    auto steps = 0;
    if(tools::log->level() == spdlog::level::trace)
        tools::log->trace("Normalizing state | old norm = {:.16f} | pos {} | dir {} | bond_limit {} | bond dims {}", tools::finite::measure::norm(state), pos,
                          dir, bond_limit.value(), tools::finite::measure::bond_dimensions(state));

    // Start with SVD at the current center position
    // NOTE: You have thought that this is unnecessary and removed it, only to find bugs much later.
    //       In particular, the bond dimension will shrink too much when doing projections, if this step is skipped.
    //       This makes sure chiL and chiR differ at most by factor spin_dim when we start the normalization
    if(pos >= 0) {
        auto &mps = state.get_mps_site(pos);
        // Make sure that the bond dimension does not increase faster than spin_dim per site
        long bond_new = std::min(bond_limit.value(), mps.spin_dim() * std::min(mps.get_chiL(), mps.get_chiR()));
        tools::finite::mps::merge_multisite_mps(state, mps.get_M(), {static_cast<size_t>(pos)}, pos, bond_new, svd_settings, LogPolicy::QUIET);
        if constexpr(settings::debug) mps.assert_identity();
    }
    // Now we can move around the chain until we return to the original status
    while(steps++ < 2 or not state.position_is_at(pos, dir, cnt)) move_center_point_single_site(state, bond_limit.value(), svd_settings);
    state.clear_measurements();
    state.clear_cache();
    auto norm = tools::finite::measure::norm(state);
    if(tools::log->level() == spdlog::level::trace)
        tools::log->trace("Normalized  state | new norm = {:.16f} | pos {} | dir {} | bond_limit {} | bond dims {}", norm, pos, dir, bond_limit.value(),
                          tools::finite::measure::bond_dimensions(state));
    if(std::abs(norm - 1) > settings::precision::max_norm_error) {
        for(const auto &mps : state.mps_sites) {
            tools::log->warn("L ({}) | norm {:.16f} \n {}", mps->get_position(), tenx::VectorMap(mps->get_L()).norm(), mps->get_L());
            if(mps->isCenter()) tools::log->warn("LC({}) | norm {:.16f} \n {}", mps->get_position(), tenx::VectorMap(mps->get_LC()).norm(), mps->get_LC());
            mps->assert_identity();
        }
        throw std::runtime_error(fmt::format("Norm too far from unity: {:.16f} | max allowed norm error {}", norm, settings::precision::max_norm_error));
    }
    return true;
}

void tools::finite::mps::randomize_state(StateFinite &state, StateInit init, StateInitType type, std::string_view sector, long bond_limit,
                                         bool use_eigenspinors, std::optional<long> bitfield) {
    switch(init) {
        case StateInit::RANDOM_PRODUCT_STATE: return init::random_product_state(state, type, sector, use_eigenspinors, bitfield);
        case StateInit::RANDOM_ENTANGLED_STATE: return init::random_entangled_state(state, type, sector, bond_limit, use_eigenspinors);
        case StateInit::RANDOMIZE_PREVIOUS_STATE: return init::randomize_given_state(state, type);
        case StateInit::PRODUCT_STATE_ALIGNED: return init::set_product_state_aligned(state, type, sector);
        case StateInit::PRODUCT_STATE_NEEL: return init::set_product_state_neel(state, type, sector);
    }
}

void tools::finite::mps::apply_random_paulis(StateFinite &state, const std::vector<Eigen::Matrix2cd> &paulimatrices) {
    auto [mpos, L, R] = qm::mpo::sum_of_pauli_mpo(paulimatrices, state.get_length(), RandomizerMode::SELECT1);
    tools::finite::ops::apply_mpos(state, mpos, L, R);
}

void tools::finite::mps::apply_random_paulis(StateFinite &state, const std::vector<std::string> &paulistrings) {
    std::vector<Eigen::Matrix2cd> paulimatrices;
    for(const auto &str : paulistrings) paulimatrices.emplace_back(init::get_pauli(str));
    apply_random_paulis(state, paulimatrices);
}

void tools::finite::mps::truncate_all_sites(StateFinite &state, long bond_limit, std::optional<svd::settings> svd_settings) {
    tools::log->trace("Truncating all sites to bond dimension {}", bond_limit);

    auto original_position  = state.get_position();
    auto original_direction = state.get_direction();
    // Start by truncating at the current position.
    while(true) {
        move_center_point(state, bond_limit, svd_settings);
        if(state.get_position() == original_position and state.get_direction() == original_direction) {
            // Check if all bond dimensions less than or equal to below bond_limit
            auto bond_dimensions = tools::finite::measure::bond_dimensions(state);
            if(std::all_of(bond_dimensions.begin(), bond_dimensions.end(), [bond_limit](const long &chi) { return chi <= bond_limit; })) break;
        }
    }
    state.clear_cache();
    state.clear_measurements();
    tools::log->trace("Truncated all sites");
    tools::log->warn("MUST REBUILD EDGES AFTER TRUNCATING ALL SITES");
}

void tools::finite::mps::truncate_active_sites([[maybe_unused]] StateFinite &state, [[maybe_unused]] long bond_limit,
                                               [[maybe_unused]] std::optional<svd::settings> svd_settings) {
    tools::log->warn("Truncate active sites needs an implementation");
    throw std::runtime_error("Truncate active sites needs an implementation");
}

void tools::finite::mps::truncate_next_sites([[maybe_unused]] StateFinite &state, [[maybe_unused]] long bond_limit, [[maybe_unused]] size_t num_sites,
                                             [[maybe_unused]] std::optional<svd::settings> svd_settings) {
    tools::log->warn("Truncate next sites needs an implementation");
    throw std::runtime_error("Truncate next sites needs an implementation");
}

template<typename GateType>
std::vector<size_t> generate_gate_sequence(const StateFinite &state, const std::vector<GateType> &gates, bool reverse) {
    // Generate a list of staggered indices, without assuming that gates are sorted in any way
    // Consider a sequence of short-range gates such as [0,1], [1,2], [2,3], then this function is used to generate a new sequence without overlaps:
    //
    //  * 2-site gates:
    //         * layer 0: [0,1],[2,3],[4,5]... and so on,
    //         * layer 1: [1,2],[3,4],[5,6]..., i.e. even sites first, then odd sites,
    //  * 3-site gates:
    //         * layer 0: [0,1,2], [3,4,5]...
    //         * layer 1: [1,2,3], [4,5,6]...,
    //         * layer 2: [2,3,4], [5,6,7], and so on.
    //
    // To avoid having to move/swap all the way back between layers, one can flip odd layers to generate a zig-zag pattern.
    // Then, when applying the inverse operation, both layers as well as gates within a layer are flipped
    //
    // Now consider a sequence of long-range gates such as [0,1], [0,2], ... [0,L-1], [1,2], [1,3]...., [1,L-1]
    // To generate a performant sequence of gates with respect to swaps, it's important to note which gates commute.
    // In this implementation we make the assumption that gates [i,j] and [i,k] always commute (i != j != k).
    // Thus, we get the following sequence
    //  * 2-site gates:
    //         * layer 0: [0,1],[0,2],[0,3]... then [2,3],[2,4],[2,5]... then [4,5],[4,6][4,7]... and so on
    //         * layer 1: [1,2],[1,3],[1,4]... then [3,4],[3,5],[3,6]... then [5,6],[5,7][5,8]... and so on
    //
    // This way we get as many reverse swap cancellations as possible while maintaining some support
    // for non-commutativity on short-range interactions

    // Performance note:
    // If the state is at position L-1, and the list generated has to start from 0, then L-1 moves have
    // to be done before even starting. Additionally, if unlucky, we have to move/swap L-1 times again to return
    // to the original position.
    //    if(state.get_direction() < 0 and past_middle) state.flip_direction(); // Turn direction away from middle
    //    if(state.get_direction() > 0 and not past_middle) state.flip_direction(); // Turn direction away from middle
    auto                             t_gen     = tid::tic_scope("gen_gate_seq");
    auto                             state_len = state.get_length<long>();
    auto                             state_pos = state.get_position<long>();
    std::vector<std::vector<size_t>> layers;
    std::vector<size_t>              idx = num::range<size_t>(0, gates.size());
    while(not idx.empty()) {
        std::vector<size_t> layer;
        size_t              at = gates.at(idx.front()).pos.front();
        for(size_t i = 0; i < idx.size(); i++) {
            const auto &gate_i = gates.at(idx.at(i));
            if(gate_i.pos.front() >= at) {
                auto back_i = gate_i.pos.back();
                for(size_t j = i; j < idx.size(); j++) { // Look ahead
                    // In this part we accept a gate if pos.front() == at or pos.back() >= back_i
                    const auto &gate_j = gates.at(idx.at(j)); // gate_i == gate_j on first iteration of j, so we always accept
                    if(gate_j.pos.front() == at and gate_j.pos.back() >= back_i) layer.emplace_back(idx.at(j));
                }
                at = gate_i.pos.back() + 1;
            }
        }
        // Append the layer
        layers.emplace_back(layer);

        // Remove the used elements from idx by value
        for(const auto &i : iter::reverse(layer)) idx.erase(std::remove(idx.begin(), idx.end(), i), idx.end());
    }

    // Reverse every other layer to get a zigzag pattern
    // If the state is past the middle, reverse layers 0,2,4... otherwise 1,3,5...
    size_t past_middle = state_pos > state_len / 2 ? 0 : 1;
    for(const auto &[i, l] : iter::enumerate(layers)) {
        if(num::mod<size_t>(i, 2) == past_middle) std::reverse(l.begin(), l.end());
    }
    // Move the layers into a long sequence of positions
    std::vector<size_t> gate_sequence;
    for(auto &l : layers) gate_sequence.insert(gate_sequence.end(), std::make_move_iterator(l.begin()), std::make_move_iterator(l.end()));

    // To apply inverse layers we reverse the whole sequence
    if(reverse) std::reverse(gate_sequence.begin(), gate_sequence.end());
    if(gate_sequence.size() != gates.size())
        throw std::logic_error(fmt::format("Expected gate_sequence.size() {} == gates.size() {}", gate_sequence.size(), gates.size()));
    return gate_sequence;
}

template std::vector<size_t> generate_gate_sequence(const StateFinite &state, const std::vector<qm::Gate> &gates, bool reverse);
template std::vector<size_t> generate_gate_sequence(const StateFinite &state, const std::vector<qm::SwapGate> &gates, bool reverse);

void tools::finite::mps::apply_gate(StateFinite &state, const qm::Gate &gate, Eigen::Tensor<cplx, 3> &temp, bool reverse, long bond_limit, GateMove gm,
                                    std::optional<svd::settings> svd_settings) {
    if(gate.pos.back() >= state.get_length()) throw std::logic_error(fmt::format("The last position of gate is out of bounds: {}", gate.pos));
    if(gate.was_used()) throw std::runtime_error(fmt::format("gate was already used: pos {} ", gate.pos));
    if(gm == GateMove::AUTO) gm = GateMove::OFF; // Most likely it does not pay to enable moving of center site. This is best reserved for swap gates.
    if(gm == GateMove::ON) {
        // Applying gate operators will automatically move around the center position IF either posL-1, posL or posR is a center.
        // To get this process going we need to move the center into one of these positions first.
        // Note 1: posL and posR here refers to the front/back of gate.pos.
        // Note 2: the apply operation will set the current center the left-most site.
        auto current_position = state.get_position<long>();
        long move_to_position = std::clamp<long>(current_position, static_cast<long>(gate.pos.front()) - 1l, static_cast<long>(gate.pos.back()));
        auto steps            = move_center_point_to_pos_dir(state, move_to_position, 1, bond_limit, svd_settings);
        if constexpr(settings::debug_gates)
            if(steps > 0) tools::log->info("apply_gate: moved center point {} -> {} | {} steps", current_position, move_to_position, steps);
    }

    auto center_position = gm == GateMove::ON ? static_cast<long>(gate.pos.front()) : state.get_position<long>();
    auto multisite_mps   = state.get_multisite_mps(gate.pos);
    {
        auto t_apply = tid::tic_token("apply");
        temp.resize(std::array<long, 3>{gate.op.dimension(0), multisite_mps.dimension(1), multisite_mps.dimension(2)});
        if(reverse)
            temp.device(tenx::omp::getDevice()) = gate.adjoint().contract(multisite_mps, tenx::idx({0}, {0}));
        else
            temp.device(tenx::omp::getDevice()) = gate.op.contract(multisite_mps, tenx::idx({0}, {0}));
    }

    gate.mark_as_used();
    if constexpr(settings::debug_gates)
        tools::log->trace("apply_gate: merging gate {} dims {} | center {}", gate.pos, multisite_mps.dimensions(), center_position);
    tools::finite::mps::merge_multisite_mps(state, temp, gate.pos, center_position, bond_limit, svd_settings, LogPolicy::QUIET);
}

void tools::finite::mps::apply_gates(StateFinite &state, const std::vector<Eigen::Tensor<cplx, 2>> &nsite_tensors, size_t gate_size, bool reverse,
                                     long bond_limit, GateMove gm, std::optional<svd::settings> svd_settings) {
    // Pack the two-site operators into a vector of qm::Gates
    std::vector<qm::Gate> gates;
    gates.reserve(nsite_tensors.size());
    for(const auto &[idx, op] : iter::enumerate(nsite_tensors)) {
        auto pos = num::range<size_t>(idx, idx + gate_size, 1);
        auto dim = std::vector<long>(pos.size(), 2);
        gates.emplace_back(qm::Gate(nsite_tensors[idx], pos, dim));
    }
    apply_gates(state, gates, reverse, bond_limit, gm, svd_settings);
}

void tools::finite::mps::apply_gates(StateFinite &state, const std::vector<qm::Gate> &gates, bool reverse, long bond_limit, GateMove gm,
                                     std::optional<svd::settings> svd_settings) {
    auto t_apply_gates = tid::tic_scope("apply_gates");

    if(gates.empty()) return;
    auto svd_count     = svd::solver::count ? svd::solver::count.value() : 0ll;
    auto gate_sequence = generate_gate_sequence(state, gates, reverse);
    if constexpr(settings::debug_gates)
        tools::log->trace("apply_gates: current pos {} dir {} | gate_sequence {}", state.get_position<long>(), state.get_direction(), gate_sequence);

    state.clear_cache(LogPolicy::QUIET);
    Eigen::Tensor<cplx, 3> gate_mps;
    for(const auto &idx : gate_sequence) { apply_gate(state, gates.at(idx), gate_mps, reverse, bond_limit, gm, svd_settings); }

    move_center_point_to_pos_dir(state, 0, 1, bond_limit, svd_settings);
    tools::finite::mps::normalize_state(state, bond_limit, svd_settings, NormPolicy::IFNEEDED);

    svd_count = (svd::solver::count ? svd::solver::count.value() : 0ll) - svd_count;
    tools::log->debug("apply_gates: applied {} gates | svds {} | time {:.4f}", gates.size(), svd_count, t_apply_gates->get_last_interval());
}

void tools::finite::mps::swap_sites(StateFinite &state, size_t posL, size_t posR, std::vector<size_t> &order, GateMove gm,
                                    std::optional<svd::settings> svd_settings) {
    /* The swap operation takes two neighboring sites
     *
     * (1)chiL---[ mpsL ]---chiC(2)    chiC(1)---[ mpsR ]---(2)chiR
     *              |                               |
     *           (0)dL                            (0)dR
     *
     *
     * and joins them in multisite standard form, i.e., with order
     *
     * 1) physical indices di (any number of them, contracted from left to right)
     * 2) left bond index chiL
     * 4) right bond index chiR
     *
     * (1)chiL ---[tensor]--- (2)chiR
     *               |
     *          (0)d*d*d...
     *
     * For the swap operation, we reshape this into a rank-4, splitting up the physical dimensions
     *
     * (2)chiL---[  tensor  ]---(3)chiR
     *            |        |
     *        (0)dL      (1)dR
     *
     * Then the swap operator is applied
     *
     * (2)chiL---[  tensor  ]---(3)chiR
     *            |        |
     *             \      /
     *               \  /
     *                /\      <---- Swap
     *              /    \
     *            /        \
     *            |        |
     *         (0)dR     (1)dL
     *
     *
     * For a schmidt decomposition we need to shuffle this object into standard form for schmidt decomposition
     *
     *   (1)chiL---[  tensor  ]---(3)chiR
     *              |        |
     *          (0)dL      (2)dR
     *
     * Then a schmidt decomposition returns

     * (1)chiL---[  U  ]---chiC(2) chiC(0)---[S]---chiC(1)  chiC(1)---[ V ]---(2)chiR
     *              |                                                   |
     *           (0)dR                                                (0)dL
     *
     *
     * In practice the swap operation can be performed as a reshape + shuffle
     *
     * (1)chiL ---[tensor]--- (2)chiR     reshape      (2)chiL---[  tensor  ]---(3)chiR    shuffle          (1)chiL---[  tensor  ]---(3)chiR
     *               |                      ---->                 |        |               ----->                      |        |
     *          (0)d*d*d...                                    (0)dL      (1)dR                                      (0)dR      (2)dL
     *
     *
     *  Center position:
     *  Before swapping, one of the sites is supposed to be an "AC" site:
     *      Swap: posL is AC, posR is B
     *      Rwap: posL is A, posR is AC
     *
     *  After swapping, the AC-site moves:
     *      Swap: posL becomes A, posR becomes AC
     *      Rwap: posL becomes AC, posR becomes A
     *
     *  This is handled by the merge operation, which takes the "new" center position as argument
     *      Swap: center_position = posR
     *      Rwap: center_position = posL
     *
     */
    auto t_swap = tid::tic_scope("swap");
    if(posR != posL + 1) throw std::logic_error(fmt::format("Expected posR == posL+1. Got posL {} and posR {}", posL, posR));
    if(posR != std::clamp(posR, 0ul, state.get_length<size_t>() - 1ul))
        throw std::logic_error(fmt::format("Expected posR in [0,{}]. Got {}", state.get_length() - 1, posR));
    if(posL != std::clamp(posL, 0ul, state.get_length<size_t>() - 1ul))
        throw std::logic_error(fmt::format("Expected posL in [0,{}]. Got {}", state.get_length() - 1, posL));

    if(gm == GateMove::AUTO) gm = GateMove::ON;
    auto center_position = gm == GateMove::ON ? static_cast<long>(posL) : state.get_position<long>();

    if constexpr(settings::debug_gates)
        tools::log->info("swapping sites ({}, {}) | center {} --> {} | svds {}", posL, posR, state.get_position<long>(), center_position,
                         svd::solver::count.value());

    auto                   dimL        = state.get_mps_site(posL).dimensions();
    auto                   dimR        = state.get_mps_site(posR).dimensions();
    auto                   dL          = dimL[0];
    auto                   dR          = dimR[0];
    auto                   chiL        = dimL[1];
    auto                   chiR        = dimR[2];
    auto                   bond_limit  = std::max(dL * chiL, dR * chiR);
    Eigen::Tensor<cplx, 3> swapped_mps = state.get_multisite_mps({posL, posR})
                                             .reshape(tenx::array4{dL, dR, chiL, chiR})
                                             .shuffle(tenx::array4{1, 0, 2, 3})           // swap
                                             .reshape(tenx::array3{dR * dL, chiL, chiR}); // prepare for merge

    merge_multisite_mps(state, swapped_mps, {posL, posR}, center_position, bond_limit, svd_settings, LogPolicy::QUIET);
    std::swap(order[posL], order[posR]);

    // Sanity check
    if constexpr(settings::debug_gates) {
        state.clear_cache(LogPolicy::QUIET);
        auto norm = tools::finite::measure::norm(state, true);
        tools::log->debug("after swapping            {}: labels {} | order {} | norm {:.16f}", std::vector<size_t>{posL, posR}, state.get_labels(), order,
                          norm);
        //            Eigen::Tensor<cplx, 2> idL, idR;
        //            for(const auto &mps : state.mps_sites) {
        //                idR = mps->get_M_bare().contract(mps->get_M_bare().conjugate(), tenx::idx({0, 2}, {0, 2}));
        //                idL = mps->get_M_bare().contract(mps->get_M_bare().conjugate(), tenx::idx({0, 1}, {0, 1}));
        //                tools::log->info("pos {:<2}: L {:<6} R {:<6}", mps->get_position(), tenx::isIdentity(idL, 1e-10), tenx::isIdentity(idR, 1e-10));
        //            }
        for(const auto &mps : state.mps_sites) mps->assert_identity();
    }
}

void tools::finite::mps::apply_swap_gate(StateFinite &state, qm::SwapGate &gate, Eigen::Tensor<cplx, 3> &temp, bool reverse, long bond_limit,
                                         std::vector<size_t> &order, GateMove gm, std::optional<svd::settings> svd_settings) {
    if(gate.was_used()) return;
    if(gate.pos.back() >= state.get_length()) throw std::logic_error(fmt::format("The last position of gate is out of bounds: {}", gate.pos));
    if constexpr(settings::debug_gates) tools::log->trace("apply_swap_gate: pos {} | order {} | svds {}", gate.pos, order, svd::solver::count.value());

    // Sanity check
    if constexpr(settings::debug_gates) {
        state.clear_cache(LogPolicy::QUIET);
        auto norm = tools::finite::measure::norm(state, true);
        //    for(const auto &mps : state.mps_sites) mps->assert_identity();
        tools::log->info("before applying swap gate {}: labels {} | order {} | norm {:.16f} | svds {}", gate.pos, state.get_labels(), order, norm,
                         svd::solver::count.value());
    }

    if(gm == GateMove::ON) {
        // Applying swap operators will automatically move around the center position IF either posL-1, posL or posR is a center.
        // To get this process going we need to move the center into one of these positions first.
        // Note 1: posL and posR here refers to the first swap in the gate if it exists. Otherwise they are front/back of gate.pos.
        // Note 2: both the swap and apply operations will set the current center the left-most site.
        auto current_position = state.get_position<long>();
        long move_to_position = current_position;
        if(gate.swaps.empty())
            move_to_position = std::clamp<long>(current_position, static_cast<long>(gate.pos.front()) - 1l, static_cast<long>(gate.pos.back()));
        else
            move_to_position = std::clamp<long>(current_position, static_cast<long>(gate.swaps[0].posL) - 1l, static_cast<long>(gate.swaps[0].posR));
        auto steps = move_center_point_to_pos_dir(state, move_to_position, 1, bond_limit, svd_settings);
        if constexpr(settings::debug_gates)
            if(steps > 0) tools::log->trace("apply_swap_gate: moved center point (v1) {} -> {} | {} steps", current_position, move_to_position, steps);
    }

    // with i<j, start by applying all the swap operators S(i,j), which move site i to j-1
    for(const auto &s : gate.swaps) swap_sites(state, s.posL, s.posR, order, gm, svd_settings);

    // First find the actual positions given the current order after having swapped a lot.
    std::vector<size_t> pos = gate.pos;
    if(not order.empty()) {
        for(auto &p : pos) {
            auto p_it = std::find(order.begin(), order.end(), p);
            if(p_it == order.end()) throw std::logic_error(fmt::format("state position of gate pos {} not found in {}", p, order));
            p = static_cast<size_t>(std::distance(order.begin(), p_it));
        }
        // Now pos has the positions on the state which correspond to the gate positions.
        // Example:
        //      order    == [2,3,1,4,5,0]
        //      gate.pos == [1,4]
        //      pos      == [2,3] <--- valid
        // Note that the point of swap-gates is to apply gates on nearest neighbors. Therefore, "pos" should always
        // be a vector with increments by 1 for the gate to be applicable. For instance, with the order above, gate.pos == [1,5] would yield
        // an invalid pos == [2,4].

        // Check that pos is valid
        if(pos.size() >= 2) {
            for(size_t idx = 1; idx < pos.size(); idx++) {
                if(pos[idx] != pos[idx - 1] + 1) throw std::logic_error(fmt::format("Tried to apply invalid swap pos {}", pos));
            }
        }
    }

    auto multisite_mps = state.get_multisite_mps(pos);
    {
        auto t_apply = tid::tic_token("apply");
        temp.resize(std::array<long, 3>{gate.op.dimension(0), multisite_mps.dimension(1), multisite_mps.dimension(2)});
        if(reverse)
            temp.device(tenx::omp::getDevice()) = gate.adjoint().contract(multisite_mps, tenx::idx({0}, {0}));
        else
            temp.device(tenx::omp::getDevice()) = gate.op.contract(multisite_mps, tenx::idx({0}, {0}));
    }

    gate.mark_as_used();
    long center_position = gm == GateMove::ON ? static_cast<long>(pos.front()) : state.get_position<long>(); // Previously state position

    if constexpr(settings::debug_gates)
        tools::log->trace("merging applied gate | pos {} | swapped pos {} | order {} | svds {}", gate.pos, pos, order, svd::solver::count.value());
    tools::finite::mps::merge_multisite_mps(state, temp, pos, center_position, bond_limit, svd_settings, LogPolicy::QUIET);

    // Now swap site j-1 in reverse back to i
    for(const auto &r : gate.rwaps) swap_sites(state, r.posL, r.posR, order, gm, svd_settings);

    // Sanity check
    if constexpr(settings::debug_gates) {
        state.clear_cache(LogPolicy::QUIET);
        auto norm = tools::finite::measure::norm(state, true);
        tools::log->info("after  applying swap gate {}: labels {} | order {} | norm {:.16f} | svds {}", gate.pos, state.get_labels(), order, norm,
                         svd::solver::count.value());
        //        Eigen::Tensor<cplx, 2> idR, idL;
        //        for(const auto &mps : state.mps_sites) {
        //            idR = mps->get_M_bare().contract(mps->get_M_bare().conjugate(), tenx::idx({0, 2}, {0, 2}));
        //            idL = mps->get_M_bare().contract(mps->get_M_bare().conjugate(), tenx::idx({0, 1}, {0, 1}));
        //            tools::log->info("pos {:<2}: L {:<6} R {:<6}", mps->get_position(), tenx::isIdentity(idL, 1e-10), tenx::isIdentity(idR, 1e-10));
        //        }
    }
}

void tools::finite::mps::apply_swap_gates(StateFinite &state, std::vector<qm::SwapGate> &gates, bool reverse, long bond_limit, GateMove gm,
                                          std::optional<svd::settings> svd_settings) {
    auto t_swapgate = tid::tic_scope("apply_swap_gates");
    if(gates.empty()) return;
    state.clear_cache(LogPolicy::QUIET); // So that multisite_mps does not use cache
    // Sanity check
    if constexpr(settings::debug_gates) {
        auto norm = tools::finite::measure::norm(state, true);
        tools::log->info("before applying swap gates: labels {} | norm {:.16f}", state.get_labels(), norm);
        for(const auto &mps : state.mps_sites) mps->assert_identity();
        //        Eigen::Tensor<cplx, 2> idL, idR;
        //        for(const auto &mps : state.mps_sites) {
        //            idR = mps->get_M_bare().contract(mps->get_M_bare().conjugate(), tenx::idx({0, 2}, {0, 2}));
        //            idL = mps->get_M_bare().contract(mps->get_M_bare().conjugate(), tenx::idx({0, 1}, {0, 1}));
        //            tools::log->info("pos {:<2}: L {:<6} R {:<6}", mps->get_position(), tenx::isIdentity(idL, 1e-10), tenx::isIdentity(idR, 1e-10));
        //        }
    }

    if(gm == GateMove::AUTO) {
        // It only pays to move center point if we are actually using swaps
        // If there are any swaps or rwaps to be made, then switch on moving of center sites.
        gm = GateMove::OFF;
        for(const auto &gate : gates) {
            if(not gate.swaps.empty() or not gate.rwaps.empty()) {
                gm = GateMove::ON;
                break;
            }
        }
    }

    auto                   svd_count = svd::solver::count.value();
    auto                   order     = num::range<size_t>(0ul, state.get_length<size_t>(), 1ul);
    Eigen::Tensor<cplx, 3> temp;
    size_t                 skip_count = 0;
    size_t                 swap_count = 0;
    size_t                 rwap_count = 0;

    auto gate_sequence = generate_gate_sequence(state, gates, reverse);
    for(const auto &[i, gate_idx] : iter::enumerate(gate_sequence)) {
        auto &gate = gates.at(gate_idx);
        if(i + 1 < gate_sequence.size()) skip_count += gate.cancel_rwaps(gates[gate_sequence[i + 1]].swaps);
        swap_count += gate.swaps.size();
        rwap_count += gate.rwaps.size();
        apply_swap_gate(state, gate, temp, reverse, bond_limit, order, gm, svd_settings);
    }
    move_center_point_to_pos_dir(state, 0, 1, bond_limit, svd_settings);
    tools::finite::mps::normalize_state(state, bond_limit, svd_settings, NormPolicy::IFNEEDED);

    svd_count = svd::solver::count.value() - svd_count;
    if constexpr(settings::debug_gates or settings::debug)
        tools::log->debug("apply_swap_gates: applied {} gates | swaps {} | rwaps {} | total {} | skips {} | svds {} | time {:.4f}", gates.size(), swap_count,
                          rwap_count, swap_count + rwap_count, skip_count, svd_count, t_swapgate->get_last_interval());
}
