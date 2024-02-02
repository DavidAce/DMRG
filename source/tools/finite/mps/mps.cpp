#include "math/tenx.h"
// -- (textra first)
#include "../mps.h"
#include "config/enums.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "math/cast.h"
#include "math/linalg/tensor.h"
#include "math/num.h"
#include "math/svd.h"
#include "qm/gate.h"
#include "qm/mpo.h"
#include "qm/spin.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"
#include "tools/common/split.h"
#include "tools/finite/measure.h"
#include "tools/finite/ops.h"
#include <fmt/ranges.h>

namespace settings {
    inline constexpr bool debug_merge   = false;
    inline constexpr bool debug_gates   = false;
    inline constexpr bool debug_moves   = false;
    inline constexpr bool verbose_merge = false;
    inline constexpr bool verbose_gates = false;
    inline constexpr bool verbose_moves = false;
}

bool tools::finite::mps::init::bitfield_is_valid(size_t bitfield) { return bitfield != -1ul and init::used_bitfields.count(bitfield) == 0; }

size_t tools::finite::mps::move_center_point_single_site(StateFinite &state, std::optional<svd::config> svd_cfg) {
    auto t_move = tid::tic_scope("move", tid::level::higher);
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
        long pos  = state.get_position<long>();  // If all sites are B's, then this is -1. Otherwise, this is the current "A*LC" site
        long posC = pos + state.get_direction(); // This is the site which becomes the new center position
        if(pos < -1 or pos >= state.get_length<long>()) throw except::range_error("pos out of bounds: {}", pos);
        if(posC < -1 or posC >= state.get_length<long>()) throw except::range_error("posC out of bounds: {}", posC);
        if(state.get_direction() != posC - pos) throw except::logic_error("Expected posC - pos == {}. Got {}", state.get_direction(), posC - pos);

        if constexpr(settings::verbose_moves) {
            if(posC > pos) tools::log->trace("Moving {} -> {}", pos, posC);
            if(posC < pos) tools::log->trace("Moving {} <- {}", posC, pos);
        }

        Eigen::Tensor<cplx, 1> LC(1);
        LC.setConstant(1); // Store the LC bond in a temporary. It will become a regular "L" bond later
        if(pos >= 0) LC = state.get_mps_site(pos).get_LC();

        if(state.get_direction() == 1) {
            auto  posC_ul = safe_cast<size_t>(posC);     // Cast to unsigned
            auto &mpsC    = state.get_mps_site(posC);    // This becomes the new AC (currently B)
            auto  trnc    = mpsC.get_truncation_error(); // Truncation error of the old B/new AC, i.e. bond to the right of posC,
            // Construct a single-site tensor. This is equivalent to state.get_multisite_mps(...) but avoid normalization checks.
            auto onesite_tensor = tools::common::contraction::contract_bnd_mps_temp(LC, mpsC.get_M());
            tools::finite::mps::merge_multisite_mps(state, onesite_tensor, {posC_ul}, posC, svd_cfg, LogPolicy::QUIET);
            mpsC.set_truncation_error_LC(std::max(trnc, mpsC.get_truncation_error_LC()));
        } else if(state.get_direction() == -1) {
            auto  pos_ul         = safe_cast<size_t>(pos);     // Cast to unsigned
            auto &mps            = state.get_mps_site(pos);    // This AC becomes the new B
            auto  trnc           = mps.get_truncation_error(); // Truncation error of old AC/new B, i.e. bond to the left of pos,
            auto  onesite_tensor = mps.get_M(); // No need to contract anything this time. Note that we must take a copy! Not a reference (LC is unset later)
            tools::finite::mps::merge_multisite_mps(state, onesite_tensor, {pos_ul}, posC, svd_cfg, LogPolicy::QUIET);
            if(posC >= 0) {
                auto &mpsC = state.get_mps_site(posC); // This old A is now an AC
                mpsC.set_truncation_error_LC(std::max(trnc, mpsC.get_truncation_error_LC()));
            }
        }
        state.clear_cache(LogPolicy::QUIET);
        state.clear_measurements(LogPolicy::QUIET);
        return 1; // Moved once, so return 1
    }
}

size_t tools::finite::mps::move_center_point(StateFinite &state, std::optional<svd::config> svd_cfg) {
    auto t_move = tid::tic_scope("move");
    if(state.position_is_outward_edge(2)) {
        state.flip_direction(); // Instead of moving out of the chain, just flip the direction and return
        return 0;               // No moves this time, return 0
    } else {
        long pos = state.get_position<long>();
        if(pos < -1 or pos >= state.get_length<long>()) throw except::range_error("pos out of bounds: {}", pos);

        long  posL    = state.get_direction() == 1 ? pos + 1 : pos - 1;
        long  posR    = state.get_direction() == 1 ? pos + 2 : pos;
        auto  posL_ul = safe_cast<size_t>(posL);
        auto  posR_ul = safe_cast<size_t>(posR);
        auto &mps     = state.get_mps_site();
        auto &mpsL    = state.get_mps_site(posL); // Becomes the new center position
        auto &mpsR    = state.get_mps_site(posR); // The site to the right of the new center position
        // Store the special LC bond in a temporary. It needs to be put back afterwards
        // Do the same with its truncation error
        Eigen::Tensor<cplx, 1> LC                  = mps.get_LC();
        double                 truncation_error_LC = mps.get_truncation_error_LC();
        auto                   twosite_tensor      = state.get_multisite_mps({posL_ul, posR_ul});
        tools::finite::mps::merge_multisite_mps(state, twosite_tensor, {static_cast<size_t>(posL), static_cast<size_t>(posR)}, safe_cast<long>(posL), svd_cfg,
                                                LogPolicy::QUIET);
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

size_t tools::finite::mps::move_center_point_to_pos(StateFinite &state, long pos, std::optional<svd::config> svd_cfg) {
    if(pos != std::clamp<long>(pos, -1l, state.get_length<long>() - 1))
        throw except::logic_error("move_center_point_to_pos: Given pos [{}]. Expected range [-1,{}]", pos, state.get_length<long>() - 1);
    if((state.get_direction() < 0 and pos > state.get_position<long>()) or //
       (state.get_direction() > 0 and pos < state.get_position<long>()))   //
        state.flip_direction();                                            // Turn direction towards new position

    size_t moves = 0;
    while(not state.position_is_at(pos)) moves += move_center_point_single_site(state, svd_cfg);
    return moves;
}

size_t tools::finite::mps::move_center_point_to_pos_dir(StateFinite &state, long pos, int dir, std::optional<svd::config> svd_cfg) {
    if(pos != std::clamp<long>(pos, -1l, state.get_length<long>() - 1))
        throw except::logic_error("move_center_point_to_pos_dir: Given pos [{}]. Expected range [-1,{}]", pos, state.get_length<long>() - 1);
    if(std::abs(dir) != 1) throw except::logic_error("move_center_point_to_pos_dir: dir must be 1 or -1");
    if((state.get_direction() < 0 and pos > state.get_position<long>()) or //
       (state.get_direction() > 0 and pos < state.get_position<long>()))   //
        state.flip_direction();                                            // Turn direction towards new position
    size_t moves = 0;
    while(not state.position_is_at(pos)) moves += move_center_point_single_site(state, svd_cfg);
    if(dir != state.get_direction()) state.flip_direction();
    return moves;
}

size_t tools::finite::mps::move_center_point_to_inward_edge(StateFinite &state, std::optional<svd::config> svd_cfg) {
    // Flip if past middle and going in the other direction
    if(state.get_direction() < 0 and state.get_position() > state.get_length() / 2) state.flip_direction();
    if(state.get_direction() > 0 and state.get_position() < state.get_length() / 2) state.flip_direction();
    size_t moves = 0;
    while(not state.position_is_inward_edge()) moves += move_center_point_single_site(state, svd_cfg);
    return moves;
}

size_t tools::finite::mps::move_center_point_to_middle(StateFinite &state, std::optional<svd::config> svd_cfg) {
    size_t moves = 0;
    while(not state.position_is_the_middle_any_direction()) moves += move_center_point_single_site(state, svd_cfg);
    return moves;
}

size_t tools::finite::mps::merge_multisite_mps(StateFinite &state, const Eigen::Tensor<cplx, 3> &multisite_mps, const std::vector<size_t> &sites,
                                               long center_position, std::optional<svd::config> svd_cfg, std::optional<LogPolicy> logPolicy) {
    auto t_merge          = tid::tic_scope("merge", tid::level::higher);
    auto current_position = state.get_position<long>();
    auto moves            = static_cast<size_t>(std::abs(center_position - current_position));
    if constexpr(settings::debug)
        if(logPolicy == LogPolicy::NORMAL)
            tools::log->trace("merge_multisite_mps: sites {} | dimensions {} | center {} -> {} | {}", sites, multisite_mps.dimensions(), current_position,
                              center_position, state.get_labels());
    if constexpr(settings::verbose_merge)
        tools::log->trace("merge_multisite_mps: sites {} | dimensions {} | center {} -> {} | {}", sites, multisite_mps.dimensions(), current_position,
                          center_position, state.get_labels());

    // Some sanity checks
    if(multisite_mps.dimension(1) != state.get_mps_site(sites.front()).get_chiL())
        throw except::logic_error("merge_multisite_mps: mps dim1 {} != chiL {} on left-most site", multisite_mps.dimension(1),
                                  state.get_mps_site(sites.front()).get_chiL(), sites.front());

    if(multisite_mps.dimension(2) != state.get_mps_site(sites.back()).get_chiR())
        throw except::logic_error("merge_multisite_mps: mps dim2 {} != chiR {} on right-most site", multisite_mps.dimension(2),
                                  state.get_mps_site(sites.back()).get_chiR(), sites.back());
    if constexpr(settings::debug_merge or settings::debug) {
        auto t_dbg = tid::tic_scope("debug");
        // Never allow nan's in the multisite_mps
        if(tenx::hasNaN(multisite_mps))
            throw except::runtime_error("merge_multisite_mps: multisite_mps has nan's:\n"
                                        "sites            :{}\n"
                                        "center_position  :{}\n"
                                        "current_position :{}\n"
                                        "multisite_mps    :\n{}",
                                        sites, center_position, current_position, linalg::tensor::to_string(multisite_mps, 3, 6));

        if(state.has_nan())
            throw except::runtime_error("merge_multisite_mps: state has nan's:\n"
                                        "sites            :{}\n"
                                        "center_position  :{}\n"
                                        "current_position :{}\n"
                                        "multisite_mps    :\n{}",
                                        sites, center_position, current_position, linalg::tensor::to_string(multisite_mps, 3, 6));

        // We have to allow non-normalized multisite mps! Otherwise, we won't be able to make them normalized
        auto norm = tenx::VectorCast(multisite_mps).norm();
        if(std::abs(norm - 1) > 1e-8) tools::log->debug("merge_multisite_mps: Multisite mps for positions {} has norm far from unity: {:.16f}", sites, norm);
    }

    // Can't set center on one of sites if the current center is too far away: we would end up with interleaved A's and B sites
    bool center_in_sites = center_position == std::clamp<long>(center_position, safe_cast<long>(sites.front()), safe_cast<long>(sites.back()));
    bool center_in_range = current_position == std::clamp<long>(current_position, safe_cast<long>(sites.front()) - 1, safe_cast<long>(sites.back()));
    if(center_in_sites and not center_in_range)
        throw except::runtime_error("merge_multisite_mps: cannot merge multisite_mps {} with new center at {}: current center {} is too far", sites,
                                    center_position, current_position);

    long              spin_prod = 1;
    std::vector<long> spin_dims;
    spin_dims.reserve(sites.size());
    for(const auto &pos : sites) {
        spin_dims.emplace_back(state.get_mps_site(pos).spin_dim());
        spin_prod *= spin_dims.back();
    }
    if(spin_prod != multisite_mps.dimension(0))
        throw except::runtime_error("merge_multisite_mps: multisite mps dim0 {} != spin_prod {}", multisite_mps.dimension(0), spin_prod);

    // Hold LC if moving. This should be placed in an L-slot later
    std::optional<stash<Eigen::Tensor<cplx, 1>>> lc_move = std::nullopt;
    if(center_position != current_position and current_position >= 0) {
        auto &mps      = state.get_mps_site(current_position); // Guaranteed to have LC since that is the definition of current_position
        auto  pos_back = safe_cast<long>(sites.back());
        auto  pos_frnt = safe_cast<long>(sites.front());
        auto  pos_curr = safe_cast<size_t>(current_position);

        // Detect right-move
        if(center_position > current_position) { // This AC will become an A (AC moves to the right)
            if(center_position != std::clamp(center_position, pos_frnt, pos_back))
                throw except::logic_error("merge_multisite_mps: right-moving new center position {} must be in sites {}", center_position, sites);

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
                throw except::logic_error("merge_multisite_mps: left-moving new center position {} is out of range [{}]+{}", center_position, pos_frnt - 1,
                                          sites);
            if(current_position > pos_back + 1)
                throw except::logic_error("merge_multisite_mps: left-moving current position {} is out of range {}+[{}]", current_position, sites,
                                          pos_back + 1);

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

    if constexpr(settings::verbose_merge)
        if(svd_cfg) tools::log->trace("merge_multisite_mps: splitting sites {} | {}", sites, svd_cfg->to_string());

    // Split the multisite mps into single-site mps objects
    auto mps_list = tools::common::split::split_mps(multisite_mps, spin_dims, sites, center_position, svd_cfg);

    // Sanity checks
    if(sites.size() != mps_list.size())
        throw std::runtime_error(
            fmt::format("merge_multisite_mps: number of sites mismatch: sites.size() {} != mps_list.size() {}", sites.size(), mps_list.size()));

    // Fuse the split-up mps components into the current state
    for(auto &mps_src : mps_list) {
        auto  pos     = mps_src.get_position();
        auto &mps_tgt = state.get_mps_site(pos);

        // inject lc_move if there is any waiting
        if(lc_move and pos == lc_move->pos_dst) { mps_src.set_L(lc_move->data, lc_move->error); }

        // No need to copy the svd truncation errors if no svd config was given. This
        if(not svd_cfg) mps_src.drop_stashed_errors();

        mps_tgt.fuse_mps(mps_src);
        state.tag_site_normalized(pos, true); // Fused site is normalized

        // Now take stashes for neighboring sites
        // Note that if U or V are rectangular and pushed onto the next site, that next site loses its normalization, unless
        // we are pushing across the center matrix.
        // Tagging it as non-normalized lets us determine whether a full normalization is required later.
        if(pos < state.get_length() - 1) {
            auto &mps_nbr  = state.get_mps_site(pos + 1);
            auto  old_chiL = mps_nbr.get_chiL();
            mps_nbr.take_stash(mps_src);                                                                                 // Take stashed S,V (and possibly LC)
            if(mps_nbr.get_label() != "B" and mps_nbr.get_chiL() != old_chiL) state.tag_site_normalized(pos + 1, false); // Normalization may have been ruined
        }
        if(pos > 0) {
            auto &mps_nbr  = state.get_mps_site(pos - 1);
            auto  old_chiR = mps_nbr.get_chiR();
            mps_nbr.take_stash(mps_src);                                                                                 // Take stashed U,S (and possibly LC)
            if(mps_nbr.get_label() == "B" and mps_nbr.get_chiR() != old_chiR) state.tag_site_normalized(pos - 1, false); // Normalization may have been ruined
        }
        mps_src.drop_stash(); // Discard whatever is left stashed at the edge (this normalizes the state)
    }

    current_position = state.get_position<long>();
    if(current_position != center_position)
        throw except::logic_error("Center position mismatch {} ! {}\nLabels: {}", current_position, center_position, state.get_labels());
    state.clear_cache(LogPolicy::QUIET);
    state.clear_measurements(LogPolicy::QUIET);
    if constexpr(settings::debug or settings::debug_merge) {
        auto t_dbg = tid::tic_scope("debug");
        for(const auto &pos : sites) state.get_mps_site(pos).assert_validity();
        for(const auto &pos : sites) state.get_mps_site(pos).assert_normalized();
    }
    return moves;
}

bool tools::finite::mps::normalize_state(StateFinite &state, std::optional<svd::config> svd_cfg, NormPolicy norm_policy) {
    // When a state needs to be normalized it's enough to "move" the center position around the whole chain.
    // Each move performs an SVD decomposition which leaves unitaries behind, effectively normalizing the state.
    // NOTE! It IS important to start with the current position.

    if(norm_policy == NormPolicy::IFNEEDED) {
        // We may only go ahead with a normalization if it's really needed.
        tools::log->trace("normalize_state: checking if needed");
        if(state.is_normalized_on_all_sites()) return false; // Return false, i.e. did "not" perform a normalization.
        // Otherwise, we just do the normalization
    }

    // Save the current position, direction and center status
    auto dir   = state.get_direction();
    auto pos   = state.get_position<long>();
    auto cnt   = pos >= 0;
    auto steps = 0;
    if(tools::log->level() <= spdlog::level::debug)
        tools::log->debug("normalize_state: old local norm = {:.16f} | pos {} | dir {} | bond dims {}", tools::finite::measure::norm(state), pos, dir,
                          tools::finite::measure::bond_dimensions(state));

    // Start with SVD at the current center position
    // NOTE: You have thought that this is unnecessary and removed it, only to find bugs much later.
    //       In particular, the bond dimension will shrink too much when doing projections, if this step is skipped.
    //       This makes sure chiL and chiR differ at most by factor spin_dim when we start the normalization
    if(pos >= 0) {
        auto &mps = state.get_mps_site(pos);
        // Make sure that the bond dimension does not increase faster than spin_dim per site
        tools::finite::mps::merge_multisite_mps(state, mps.get_M(), {static_cast<size_t>(pos)}, pos, svd_cfg, LogPolicy::QUIET);
    }
    // Now we can move around the chain until we return to the original status
    while(steps++ < 2 or not state.position_is_at(pos, dir, cnt)) move_center_point_single_site(state, svd_cfg);
    state.clear_measurements();
    state.clear_cache();
    if(not state.is_normalized_on_all_sites()) {
        for(const auto &mps : state.mps_sites) {
            bool normalized_tag = state.get_normalization_tags()[mps->get_position<size_t>()];
            tools::log->warn("{} | is_normalized {:<7} | L norm {:.16f} | norm tag {}", mps->get_tag(), mps->is_normalized(),
                             tenx::VectorMap(mps->get_L()).norm(), normalized_tag);
            if(mps->isCenter()) tools::log->warn("LC({}) | norm {:.16f}", mps->get_position(), tenx::VectorMap(mps->get_LC()).norm());
        }
        throw except::runtime_error("normalize_state: normalization failed. state norm {:.16f} | max allowed norm error {:.2e} | norm tags {}",
                                    tools::finite::measure::norm(state), settings::precision::max_norm_error, state.get_normalization_tags());
    }

    if(svd_cfg and svd_cfg->rank_max and state.find_largest_bond() > svd_cfg->rank_max.value())
        throw except::logic_error("normalize_state: a bond dimension exceeds bond limit: {} > {}", tools::finite::measure::bond_dimensions(state),
                                  svd_cfg->rank_max.value());
    if(tools::log->level() <= spdlog::level::debug)
        tools::log->debug("normalize_state: new local norm = {:.16f} | pos {} | dir {} | bond dims {}", tools::finite::measure::norm(state), pos, dir,
                          tools::finite::measure::bond_dimensions(state));
    return true;
}

void tools::finite::mps::initialize_state(StateFinite &state, StateInit init, StateInitType type, std::string_view axis, bool use_eigenspinors, long bond_lim,
                                          std::string &pattern) {
    switch(init) {
        case StateInit::RANDOM_PRODUCT_STATE: return init::random_product_state(state, type, axis, use_eigenspinors, pattern);
        case StateInit::RANDOM_ENTANGLED_STATE: return init::random_entangled_state(state, type, axis, bond_lim, use_eigenspinors);
        case StateInit::RANDOMIZE_PREVIOUS_STATE: return init::randomize_given_state(state, type);
        case StateInit::MIDCHAIN_SINGLET_NEEL_STATE: return init::set_midchain_singlet_neel_state(state, type, axis, pattern);
        case StateInit::PRODUCT_STATE_DOMAIN_WALL: return init::set_product_state_domain_wall(state, type, axis, pattern);
        case StateInit::PRODUCT_STATE_ALIGNED: return init::set_product_state_aligned(state, type, axis, pattern);
        case StateInit::PRODUCT_STATE_NEEL: return init::set_product_state_neel(state, type, axis, pattern);
        case StateInit::PRODUCT_STATE_NEEL_SHUFFLED: return init::set_product_state_neel_shuffled(state, type, axis, pattern);
        case StateInit::PRODUCT_STATE_NEEL_DISLOCATED: return init::set_product_state_neel_dislocated(state, type, axis, pattern);
        case StateInit::PRODUCT_STATE_PATTERN: return init::set_product_state_on_axis_using_pattern(state, type, axis, pattern);
    }
}

void tools::finite::mps::apply_random_paulis(StateFinite &state, const std::vector<Eigen::Matrix2cd> &paulimatrices) {
    auto [mpos, L, R] = qm::mpo::sum_of_pauli_mpo(paulimatrices, state.get_length(), RandomizerMode::SELECT1);
    tools::finite::ops::apply_mpos(state, mpos, L, R);
}

void tools::finite::mps::apply_random_paulis(StateFinite &state, const std::vector<std::string> &paulistrings) {
    std::vector<Eigen::Matrix2cd> paulimatrices;
    for(const auto &str : paulistrings) paulimatrices.emplace_back(qm::spin::half::get_pauli(str));
    apply_random_paulis(state, paulimatrices);
}

template<typename GateType>
std::vector<size_t> tools::finite::mps::generate_gate_sequence(const StateFinite &state, const std::vector<GateType> &gates, CircuitOp cop,
                                                               bool range_long_to_short) {
    // Generate a list of staggered indices, without assuming that gates are sorted in any way
    // Consider a sequence of short-range gates such as [0,1], [1,2], [2,3], then this function is used to generate a new sequence without overlaps:
    //
    //  * short range 2-site gates:
    //         * layer 0: [0,1],[2,3],[4,5]... and so on,
    //         * layer 1: [1,2],[3,4],[5,6]..., i.e. even sites first, then odd sites,
    //  * short range 3-site gates:
    //         * layer 0: [0,1,2], [3,4,5]...
    //         * layer 1: [1,2,3], [4,5,6]...,
    //         * layer 2: [2,3,4], [5,6,7], and so on.
    //
    // To avoid having to move/swap all the way back between layers, one can flip odd layers to generate a zigzag pattern.
    // Then, when applying the adjoint operation, both the layer and gates within that layer are flipped
    //
    // Now consider a sequence of long-range gates such as [0,1], [0,2], ... [0,L-1], [1,2], [1,3]...., [1,L-1]
    // To generate a performant sequence of gates with respect to swaps, it's important to note which gates commute.
    // In this implementation we make the assumption that gates [i,j] and [i,k] always commute (i != j != k).
    // Thus, we get the following sequence
    //  * long range 2-site gates:
    //         * layer 0: [0,1],[0,2],[0,3]... then [2,3],[2,4],[2,5]... then [4,5],[4,6][4,7]... and so on
    //         * layer 1: [1,2],[1,3],[1,4]... then [3,4],[3,5],[3,6]... then [5,6],[5,7][5,8]... and so on
    //
    // This way we get as many reverse swap cancellations as possible while maintaining some support
    // for non-commutativity on short-range interactions
    //
    // For performance reasons we may want to apply gates from long to short range:
    // * flipped long range 2-site gates:
    //        * layer 0: ... [0,3],[0,2],[0,1] then ... [2,5],[2,4],[2,3] then ... [4,7],[4,6][4,5] and so on
    //        * layer 1: ... [1,4],[1,3],[1,2] then ... [3,6],[3,5],[3,4] then ... [5,8],[5,7][5,6] and so on
    //
    // Performance note:
    // If the state is at position L-1, and the list generated has to start from 0, then L-1 moves have
    // to be done before even starting. Additionally, if unlucky, we have to move/swap L-1 times again to return
    // to the original position.
    //    if(state.get_direction() < 0 and past_middle) state.flip_direction(); // Turn direction away from middle
    //    if(state.get_direction() > 0 and not past_middle) state.flip_direction(); // Turn direction away from middle
    auto                                          t_gen     = tid::tic_scope("gen_gate_seq", tid::level::highest);
    auto                                          state_len = state.get_length<long>();
    auto                                          state_pos = state.get_position<long>();
    std::vector<std::vector<std::vector<size_t>>> layers; // Layer --> group --> idx
    std::vector<size_t>                           idx = num::range<size_t>(0, gates.size());
    while(not idx.empty()) {
        std::vector<std::vector<size_t>> layer; // Group --> idx
        std::vector<size_t>              idx_used;
        size_t                           at = gates.at(idx.front()).pos.front(); // Left most leg of the gate at idx[0].
        for(size_t i = 0; i < idx.size(); i++) {
            const auto &gate_i = gates.at(idx.at(i)); // Gate at idx[i]
            if(gate_i.pos.front() >= at) {
                auto                back_i = gate_i.pos.back(); // Right-most leg of gate idx[i]
                std::vector<size_t> group_i;                    // List of idx with gates starting at the same position (i.e. same "at")
                for(size_t j = i; j < idx.size(); j++) {        // Look ahead
                    // In this part we accept a gate if gate j is further ahead than i, without overlap.
                    const auto &gate_j  = gates.at(idx.at(j)); // gate_i == gate_j on first iteration of j, so we always accept
                    auto        front_j = gate_j.pos.front();
                    auto        back_j  = gate_j.pos.back();
                    if(front_j == at and back_j >= back_i) {
                        group_i.emplace_back(idx.at(j));
                        idx_used.emplace_back(idx.at(j));
                    }
                }
                // Append the group in reverse order to the current layer
                if(range_long_to_short) std::reverse(group_i.begin(), group_i.end());
                layer.emplace_back(group_i);
                at = back_i + 1; // Move the position forward to a gate that doesn't overlap with the gate at idx[0]
            }
        }
        // Append the layer
        layers.emplace_back(layer);

        // Remove the used elements from idx by value
        for(const auto &i : iter::reverse(idx_used)) idx.erase(std::remove(idx.begin(), idx.end(), i), idx.end());
    }

    // Reverse every other layer to get a zigzag pattern
    // If the state is past the middle, reverse layers 0,2,4... otherwise 1,3,5...
    size_t past_middle = state_pos > state_len / 2 ? 0 : 1;
    for(const auto &[i, layer] : iter::enumerate(layers)) {
        if(num::mod<size_t>(i, 2) == past_middle) std::reverse(layer.begin(), layer.end());
    }

    // To apply inverse we reverse the layers
    // Note that we don't need to reverse groups because we assumed that these commute
    bool reverse = cop != CircuitOp::NONE;
    if(reverse) {
        for(auto &layer : layers) std::reverse(layer.begin(), layer.end());
        std::reverse(layers.begin(), layers.end());
    }

    // Move the layers into a long sequence of positions
    std::vector<size_t> gate_sequence;
    for(auto &layer : layers)
        for(auto &group : layer) gate_sequence.insert(gate_sequence.end(), std::make_move_iterator(group.begin()), std::make_move_iterator(group.end()));

    if(gate_sequence.size() != gates.size()) throw except::logic_error("gate_sequence.size() {} != gates.size() {}", gate_sequence.size(), gates.size());
    return gate_sequence;
}

template std::vector<size_t> tools::finite::mps::generate_gate_sequence(const StateFinite &state, const std::vector<qm::Gate> &gates, CircuitOp cop,
                                                                        bool range_long_to_short);
template std::vector<size_t> tools::finite::mps::generate_gate_sequence(const StateFinite &state, const std::vector<qm::SwapGate> &gates, CircuitOp cop,
                                                                        bool range_long_to_short);

void tools::finite::mps::apply_gate(StateFinite &state, const qm::Gate &gate, Eigen::Tensor<cplx, 3> &temp, GateOp gop, GateMove gm,
                                    std::optional<svd::config> svd_cfg) {
    if(gate.pos.back() >= state.get_length()) throw except::logic_error("The last position of gate is out of bounds: {}", gate.pos);
    if(gm == GateMove::AUTO) gm = GateMove::ON;
    auto old_posC = state.get_position<long>();
    auto old_svds = svd::solver::get_count();

    if constexpr(settings::verbose_gates)
        tools::log->trace("apply_gate: status  pos {} | op {} | gm {} | center {} | svds {} | labels {}", gate.pos, enum2sv(gop), enum2sv(gm), old_posC,
                          old_svds, state.get_labels());
    if(gm == GateMove::ON) {
        // When gm == GateMove::ON, applying swap operators automatically moves the center position onto posL. This can
        // happen IFF the center position is already on one of posL-1, posL or posR.
        // To get this process going we need to move the center into one of these positions first.
        // Note 1: posL and posR here refer to the first swap in the gate if it exists. Otherwise, they are front/back of gate.pos.
        // Note 2: both the swap and apply operations will set the current center on the left-most site of the gate.
        auto tgt_gate = gate.pos;
        auto tgt_posC = std::clamp<long>(old_posC, safe_cast<long>(tgt_gate.front()) - 1l, safe_cast<long>(tgt_gate.back()));
        auto num_step = move_center_point_to_pos_dir(state, tgt_posC, 1, svd_cfg);
        if constexpr(settings::verbose_gates)
            if(num_step > 0)
                tools::log->trace("apply_gate: moved   pos {} | op {} | gm {} | center {} -> {} | tgt {} | steps {} | svds {}", gate.pos, enum2sv(gop),
                                  enum2sv(gm), old_posC, tgt_posC, tgt_gate, num_step, svd::solver::get_count());
        old_posC = tgt_posC; // Update the current position
    }

    // Generate the mps object on which to apply the gate
    auto multisite_mps = state.get_multisite_mps(gate.pos);
    /* Apply the gate
        1 --- [mps] --- 2
                |
                0
                1
                |
               [G]
                |
                0

     */

    {
        auto        t_apply = tid::tic_token("apply", tid::level::highest);
        auto        &threads = tenx::threads::get();
        const auto &gateop  = gate.unaryOp(gop); // Apply any pending modifier like transpose, adjoint or conjugation
        temp.resize(std::array<long, 3>{gateop.dimension(1), multisite_mps.dimension(1), multisite_mps.dimension(2)});
        temp.device(*threads->dev) = gateop.contract(multisite_mps, tenx::idx({1}, {0}));
    }
    if constexpr(settings::verbose_gates)
        tools::log->trace("apply_gate: applied pos {} | op {} | gm {} | svds {} | bond {} | trnc {:.3e}", gate.pos, enum2sv(gop), enum2sv(gm),
                          svd::solver::get_count(), state.get_mps_site(gate.pos.front()).get_chiR(), state.get_truncation_error(gate.pos.front()));

    auto new_posC = gm == GateMove::ON ? safe_cast<long>(gate.pos.front()) : old_posC;
    tools::finite::mps::merge_multisite_mps(state, temp, gate.pos, new_posC, svd_cfg, LogPolicy::QUIET);
    if constexpr(settings::verbose_gates)
        tools::log->trace("apply_gate: merged  pos {} | op {} | gm {} | center {} -> {} | svds {} | bond {} | trnc {:.3e}", gate.pos, enum2sv(gop), enum2sv(gm),
                          old_posC, new_posC, svd::solver::get_count(), state.get_mps_site(gate.pos.front()).get_chiR(),
                          state.get_truncation_error(gate.pos.front()));
}

void tools::finite::mps::apply_gates(StateFinite &state, const std::vector<Eigen::Tensor<cplx, 2>> &nsite_tensors, size_t gate_size, CircuitOp cop,
                                     bool moveback, GateMove gm, std::optional<svd::config> svd_cfg) {
    // Pack the two-site operators into a vector of qm::Gates
    std::vector<qm::Gate> gates;
    gates.reserve(nsite_tensors.size());
    for(const auto &[idx, op] : iter::enumerate(nsite_tensors)) {
        auto pos = num::range<size_t>(idx, idx + gate_size, 1);
        auto dim = std::vector<long>(pos.size(), 2);
        gates.emplace_back(qm::Gate(nsite_tensors[idx], pos, dim));
    }
    apply_gates(state, gates, cop, moveback, gm, svd_cfg);
}

void tools::finite::mps::apply_gates(StateFinite &state, const std::vector<qm::Gate> &gates, CircuitOp cop, bool moveback, GateMove gm,
                                     std::optional<svd::config> svd_cfg) {
    auto t_apply_gates = tid::tic_scope("apply_gates", tid::level::higher);

    if(gates.empty()) return;
    auto svd_count     = svd::solver::get_count();
    auto gate_sequence = generate_gate_sequence(state, gates, cop);
    if constexpr(settings::verbose_gates)
        tools::log->trace("apply_gates: current pos {} dir {} | gate_sequence {}", state.get_position<long>(), state.get_direction(), gate_sequence);

    state.clear_cache(LogPolicy::QUIET);
    Eigen::Tensor<cplx, 3> gate_mps;
    GateOp                 gop = GateOp::NONE;
    switch(cop) {
        case CircuitOp::NONE: gop = GateOp::NONE; break;
        case CircuitOp::ADJ: gop = GateOp::ADJ; break;
        case CircuitOp::TRN: gop = GateOp::TRN; break;
    }
    for(const auto &idx : gate_sequence) apply_gate(state, gates.at(idx), gate_mps, gop, gm, svd_cfg);

    if(moveback) move_center_point_to_pos_dir(state, 0, 1, svd_cfg);
    svd_count = svd::solver::get_count() - svd_count;
    if constexpr(settings::verbose_gates)
        tools::log->debug("apply_gates: applied {} gates | svds {} | time {:.4f}", gates.size(), svd_count, t_apply_gates->get_last_interval());
}

void tools::finite::mps::apply_circuit(StateFinite &state, const std::vector<std::vector<qm::Gate>> &circuit, CircuitOp cop, bool moveback, GateMove gm,
                                       std::optional<svd::config> svd_cfg) {
    //    tools::log->debug("Applying circuit | adjoint: {} | gm {}", adjoint, enum2sv(gm));
    switch(cop) {
        case CircuitOp::ADJ:
        case CircuitOp::TRN: {
            for(const auto &gates : iter::reverse(circuit)) apply_gates(state, gates, cop, moveback, gm, svd_cfg);
            break;
        }
        case CircuitOp::NONE: {
            for(const auto &gates : circuit) apply_gates(state, gates, cop, moveback, gm, svd_cfg);
            break;
        }
    }
}

template<typename T>
std::vector<T> get_site_idx(const std::vector<size_t> &sites, const std::vector<size_t> &pos) {
    std::vector<T> idx;
    if(sites.empty()) {
        for(const auto &p : pos) idx.emplace_back(static_cast<T>(p));
    } else {
        for(const auto &p : pos) {
            auto p_it = std::find(sites.begin(), sites.end(), p);
            if(p_it == sites.end()) throw except::logic_error("state position of gate pos {} not found in {}", p, sites);
            idx.emplace_back(static_cast<T>(std::distance(sites.begin(), p_it)));
        }
    }
    return idx;
}

void tools::finite::mps::swap_sites(StateFinite &state, size_t posL, size_t posR, std::vector<size_t> &sites, GateMove gm) {
    /* The swap operation takes two neighboring sites
     *
     * (1)chiL---[ mpsL ]---chiC(2)    chiC(1)---[ mpsR ]---(2)chiR
     *              |                               |
     *           (0)dL                            (0)dR
     *
     *
     * and joins them in multisite standard form, i.e., in the order:
     *
     *      1) physical indices di (any number of them, contracted from left to right)
     *      2) left bond index chiL
     *      3) right bond index chiR
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
    if(posR != posL + 1) throw except::logic_error("Expected posR == posL+1. Got posL {} and posR {}", posL, posR);
    if(posR != std::clamp(posR, 0ul, state.get_length<size_t>() - 1ul))
        throw except::logic_error("Expected posR in [0,{}]. Got {}", state.get_length() - 1, posR);
    if(posL != std::clamp(posL, 0ul, state.get_length<size_t>() - 1ul))
        throw except::logic_error("Expected posL in [0,{}]. Got {}", state.get_length() - 1, posL);

    if(gm == GateMove::AUTO) gm = GateMove::ON;
    auto           old_svd           = svd::solver::get_count();
    auto           old_pos           = state.get_position<long>();
    auto           dimL              = state.get_mps_site(posL).dimensions();
    auto           dimR              = state.get_mps_site(posR).dimensions();
    auto           dL                = dimL[0];
    auto           dR                = dimR[0];
    auto           chiL              = dimL[1];
    auto           chiR              = dimR[2];
    auto           rsh4              = tenx::array4{dL, dR, chiL, chiR};
    constexpr auto shf4              = tenx::array4{1, 0, 2, 3};
    auto           rsh3              = tenx::array3{dR * dL, chiL, chiR};
    auto           swapped_mps       = Eigen::Tensor<cplx, 3>(rsh3);
    auto           &threads          = tenx::threads::get();
    swapped_mps.device(*threads->dev) = state.get_multisite_mps({posL, posR})
                                           .reshape(rsh4)
                                           .shuffle(shf4)  // swap
                                           .reshape(rsh3); // prepare for merge
    auto new_pos = old_pos;
    if(gm == GateMove::ON)
        new_pos = safe_cast<long>(posL); // The benefit of GateMove::ON is to prefer "AC-B" splits that require a single SVD as often as possible

    // This SVD shouldn't modify the current mps, so no truncation here (no svd config)
    merge_multisite_mps(state, swapped_mps, {posL, posR}, new_pos, std::nullopt, LogPolicy::QUIET);
    std::swap(sites[posL], sites[posR]);

    // Sanity check
    if constexpr(settings::verbose_gates) {
        tools::log->trace("swap_sites     : swapped pos [{}, {}] | idx {} | gm {} | center {} -> {} | svds {} -> {} | labels {} | sites {}", posL, posR,
                          get_site_idx<size_t>(sites, {posL, posR}), enum2sv(gm), old_pos, new_pos, old_svd, svd::solver::get_count(), state.get_labels(),
                          sites);
    }
    if constexpr(settings::debug_gates)
        for(const auto &mps : state.mps_sites) mps->assert_normalized();
}

void tools::finite::mps::apply_swap_gate(StateFinite &state, qm::SwapGate &gate, Eigen::Tensor<cplx, 3> &temp, GateOp gop, std::vector<size_t> &sites,
                                         GateMove gm, std::optional<svd::config> svd_cfg) {
    if(gate.pos.back() >= state.get_length()) throw except::logic_error("The last position of gate is out of bounds: {}", gate.pos);
    auto old_posC = state.get_position<long>();
    auto old_svds = svd::solver::get_count();
    auto pos_idxs = get_site_idx<size_t>(sites, gate.pos); // Site on the state where the gate positions are currently located

    // pos_idx has the indices pointing to where the gate positions are located on the chain right now, after many swaps.
    // Example:
    //      sites    == [0,2,3,1,4,5]
    //      gate.pos == [1,4]
    //      pos_idx  == [3,4] <--- valid if diff(pos_idx) == 1, before a gate is applied
    // Note that the point of swap gates is to apply gates on nearest neighbors. Therefore, pos_idx should always
    // be a vector of elements with increments by 1 for the gate to be applicable.

    if constexpr(settings::verbose_gates)
        tools::log->trace("apply_swap_gate: status  pos {} | idx {} | gm {} | center {} | svds {} | labels {} | sites {}", gate.pos, pos_idxs, enum2sv(gm),
                          old_posC, old_svds, state.get_labels(), sites);

    if(gm == GateMove::ON) {
        // When gm == GateMove::ON, applying swap operators automatically moves the center position onto posL. This can
        // happen IFF the center position is already on one of posL-1, posL or posR.
        // To get this process going we need to move the center into one of these positions first.
        // Note 1: posL and posR here refer to the first swap in the gate if it exists. Otherwise, they are front/back of gate.pos.
        // Note 2: both the swap and apply operations will set the current center on the left-most site of the gate.
        auto tgt_gate = gate.swaps.empty() ? gate.pos : std::vector<size_t>{gate.swaps[0].posL, gate.swaps[0].posR};
        auto tgt_posC = std::clamp<long>(old_posC, safe_cast<long>(tgt_gate.front()) - 1l, safe_cast<long>(tgt_gate.back()));
        auto num_step = move_center_point_to_pos_dir(state, tgt_posC, 1, svd_cfg);
        if constexpr(settings::verbose_gates)
            if(num_step > 0)
                tools::log->trace("apply_swap_gate: moved  pos {} | gm {} | center {} -> {} | tgt {} | steps {} | svds {}", gate.pos, enum2sv(gm), old_posC,
                                  tgt_posC, tgt_gate, num_step, svd::solver::get_count());
        old_posC = tgt_posC; // Update the current position
    }

    // with i<j, start by applying all the swap operators S(i,j), which move site i to j-1
    for(const auto &s : gate.swaps) swap_sites(state, s.posL, s.posR, sites, gm);

    // Refresh the site indices of the gate after having swapped a lot.
    pos_idxs = get_site_idx<size_t>(sites, gate.pos);

    // Check that idx elements increase by one
    if(pos_idxs.size() >= 2) {
        for(size_t i = 1; i < pos_idxs.size(); i++) {
            if(pos_idxs[i] != pos_idxs[i - 1] + 1) throw except::logic_error("Invalid pos idx {} | gate.pos {} | sites {}", pos_idxs, gate.pos, sites);
        }
    }

    // Generate the mps object on which to apply the gate
    auto multisite_mps = state.get_multisite_mps(pos_idxs);

    // Apply the gate
    {
        auto        t_apply = tid::tic_token("apply", tid::level::higher);
        const auto &gateop  = gate.unaryOp(gop); // Apply any pending modifier like transpose, adjoint or conjugation
        auto        &threads = tenx::threads::get();
        temp.resize(std::array<long, 3>{gateop.dimension(1), multisite_mps.dimension(1), multisite_mps.dimension(2)});
        temp.device(*threads->dev) = gateop.contract(multisite_mps, tenx::idx({1}, {0}));
    }
    if constexpr(settings::verbose_gates)
        tools::log->trace("apply_swap_gate: applied pos {} | idx {} | gm {} | sites {} | svds {}", gate.pos, pos_idxs, enum2sv(gm), sites,
                          svd::solver::get_count());

    // It's best to do an AC-B type of SVD split, so we put the center position on the left-most site when GateMove::ON
    long new_posC = gm == GateMove::ON ? safe_cast<long>(pos_idxs.front()) : old_posC;

    tools::finite::mps::merge_multisite_mps(state, temp, pos_idxs, new_posC, svd_cfg, LogPolicy::QUIET);
    if constexpr(settings::verbose_gates)
        tools::log->trace("apply_swap_gate: merged  pos {} | idx {} | gm {} | sites {} | svds {}", gate.pos, pos_idxs, enum2sv(gm), sites,
                          svd::solver::get_count());

    // Now swap site j-1 in reverse back to i
    for(const auto &r : gate.rwaps) swap_sites(state, r.posL, r.posR, sites, gm);

    // Sanity check
    if constexpr(settings::verbose_gates) {
        tools::log->trace("apply_swap_gate: return  pos {} | gm {} | center {} -> {} | svds {} | labels {} | sites {}", gate.pos, enum2sv(gm), old_posC,
                          new_posC, svd::solver::get_count(), state.get_labels(), sites);
    }
}

void tools::finite::mps::apply_swap_gates(StateFinite &state, std::vector<qm::SwapGate> &gates, CircuitOp cop, GateMove gm,
                                          std::optional<svd::config> svd_cfg) {
    auto t_swapgate = tid::tic_scope("apply_swap_gates", tid::level::higher);
    if(gates.empty()) return;
    state.clear_cache(LogPolicy::QUIET); // So that multisite_mps does not use cache
    // Sanity check
    if constexpr(settings::debug_gates) {
        for(const auto &mps : state.mps_sites) mps->assert_normalized();
    }
    if constexpr(settings::verbose_gates) tools::log->trace("before applying swap gates: labels {}", state.get_labels());

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

    auto                    sites = num::range<size_t>(0ul, state.get_length<size_t>(), 1ul);
    Eigen::Tensor<cplx, 3>  temp;
    [[maybe_unused]] auto   svds_count = svd::solver::get_count();
    [[maybe_unused]] size_t skip_count = 0;
    [[maybe_unused]] size_t swap_count = 0;
    [[maybe_unused]] size_t rwap_count = 0;

    auto   gate_sequence = generate_gate_sequence(state, gates, cop);
    GateOp gop           = GateOp::NONE;
    switch(cop) {
        case CircuitOp::NONE: gop = GateOp::NONE; break;
        case CircuitOp::ADJ: gop = GateOp::ADJ; break;
        case CircuitOp::TRN: gop = GateOp::TRN; break;
    }

    for(const auto &[i, gate_idx] : iter::enumerate(gate_sequence)) {
        auto &gate = gates.at(gate_idx);
        if(i + 1 < gate_sequence.size()) skip_count += gate.cancel_rwaps(gates[gate_sequence[i + 1]].swaps);
        swap_count += gate.swaps.size();
        rwap_count += gate.rwaps.size();
        apply_swap_gate(state, gate, temp, gop, sites, gm, svd_cfg);
    }

    move_center_point_to_pos_dir(state, 0, 1, svd_cfg);

    if constexpr(settings::verbose_gates) {
        svds_count = svd::solver::get_count() - svds_count;
        tools::log->debug("apply_swap_gates: applied {} gates | swaps {} | rwaps {} | total {} | skips {} | svds {} | time {:.4f}", gates.size(), swap_count,
                          rwap_count, swap_count + rwap_count, skip_count, svds_count, t_swapgate->get_last_interval());
    }
}
