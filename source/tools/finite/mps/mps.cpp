#include <math/tenx.h>
// -- (textra first)
#include "../mps.h"
#include <config/enums.h>
#include <config/settings.h>
#include <general/iter.h>
#include <math/num.h>
#include <math/svd.h>
#include <qm/gate.h>
#include <qm/mpo.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/common/split.h>
#include <tools/finite/measure.h>
#include <tools/finite/ops.h>

bool tools::finite::mps::init::bitfield_is_valid(std::optional<long> bitfield) {
    return bitfield.has_value() and bitfield.value() > 0 and init::used_bitfields.count(bitfield.value()) == 0;
}

size_t tools::finite::mps::move_center_point_single_site(StateFinite &state, long chi_lim, std::optional<svd::settings> svd_settings) {
    auto t_move = tid::tic_scope("move");
    if(state.position_is_outward_edge()) {
        if(state.get_direction() == -1 and state.get_mps_site(0l).get_chiL() != 1)
            throw std::logic_error(fmt::format("chiL at position 0 must have dimension 1, but it has dimension {}. Mps dims {}",
                                               state.get_mps_site(0l).get_chiL(), state.get_mps_site(0l).dimensions()));
        if(state.get_direction() == 1 and state.get_mps_site().get_chiR() != 1)
            throw std::logic_error(fmt::format("chiR at position {} must have dimension 1, but it has dimension {}. Mps dims {}", state.get_position(),
                                               state.get_mps_site().get_chiR(), state.get_mps_site().dimensions()));
        state.flip_direction(); // Instead of moving out of the chain, just flip the direction and return
        return 0;               // No moves this time, return 0
    } else {
        long pos  = state.get_position<long>();  // If all sites are B's, then this is -1. Otherwise this is the current "A*LC" site
        long posC = pos + state.get_direction(); // This is the site which becomes the new center position
        if(pos < -1 or pos >= state.get_length<long>()) throw std::runtime_error(fmt::format("pos out of bounds: {}", pos));
        if(posC < -1 or posC >= state.get_length<long>()) throw std::runtime_error(fmt::format("posC out of bounds: {}", posC));
        if(state.get_direction() != posC - pos) throw std::logic_error(fmt::format("Expected posC - pos == {}. Got {}", state.get_direction(), posC - pos));

        if constexpr(settings::debug_moves)
            if(posC > pos) tools::log->info("Moving {} -> {}", pos, posC);
        if constexpr(settings::debug_moves)
            if(posC < pos) tools::log->info("Moving {} <- {}", posC, pos);
        Eigen::Tensor<Scalar, 1> LC(1);
        LC.setConstant(1); // Store the LC bond in a temporary. It will become a regular "L" bond later
        if(pos >= 0) LC = state.get_mps_site(pos).get_LC();

        if(state.get_direction() == 1) {
            auto  posC_ul = static_cast<size_t>(posC);                                                       // Cast to unsigned
            auto &mpsC    = state.get_mps_site(posC);                                                        // This becomes the new center position AC
            long  chi_new = std::min(chi_lim, mpsC.spin_dim() * std::min(mpsC.get_chiL(), mpsC.get_chiR())); // Bond dimensions growth limit
            // Construct a single-site tensor. This is equivalent to state.get_multisite_mps(...) but avoid normalization checks.
            Eigen::Tensor<Scalar, 3> onesite_tensor(mpsC.dimensions()); // Allocate for contraction
            onesite_tensor.device(tenx::omp::getDevice()) = tenx::asDiagonal(LC).contract(mpsC.get_M(), tenx::idx({1}, {1})).shuffle(tenx::array3{1, 0, 2});
            tools::finite::mps::merge_multisite_tensor(state, onesite_tensor, {posC_ul}, posC, chi_new, svd_settings, LogPolicy::QUIET);
        } else if(state.get_direction() == -1) {
            auto  pos_ul         = static_cast<size_t>(pos);                                                     // Cast to unsigned
            auto &mps            = state.get_mps_site(pos);                                                      // This becomes the new B
            long  chi_new        = std::min(chi_lim, mps.spin_dim() * std::min(mps.get_chiL(), mps.get_chiR())); // Bond dimensions growth limit
            auto  onesite_tensor = mps.get_M(); // No need to contract anything this time. Note that we must take a copy! Not a reference (MC is unset later)
            tools::finite::mps::merge_multisite_tensor(state, onesite_tensor, {pos_ul}, posC, chi_new, svd_settings, LogPolicy::QUIET);
        }
        state.clear_cache(LogPolicy::QUIET);
        state.clear_measurements(LogPolicy::QUIET);
        return 1; // Moved once, so return 1
    }
}

size_t tools::finite::mps::move_center_point(StateFinite &state, long chi_lim, std::optional<svd::settings> svd_settings) {
    auto t_move = tid::tic_scope("move");
    if(state.position_is_outward_edge(2)) {
        state.flip_direction(); // Instead of moving out of the chain, just flip the direction and return
        return 0;               // No moves this time, return 0
    } else {
        long pos = state.get_position<long>();
        if(pos < -1 or pos >= state.get_length<long>()) throw std::runtime_error(fmt::format("pos out of bounds: {}", pos));
        long  posL = state.get_direction() == 1 ? pos + 1 : pos - 1;
        long  posR = state.get_direction() == 1 ? pos + 2 : pos;
        auto &mps  = state.get_mps_site();
        auto &mpsL = state.get_mps_site(posL); // Becomes the new center position
        auto &mpsR = state.get_mps_site(posR); // The site to the right of the new center position
        long  dL   = mpsL.spin_dim();
        long  dR   = mpsR.spin_dim();
        long  chiL = mpsL.get_chiL();
        long  chiR = mpsR.get_chiR();
        // Store the special LC bond in a temporary. It needs to be put back afterwards
        // Do the same with its truncation error
        Eigen::Tensor<Scalar, 1> LC                  = mps.get_LC();
        double                   truncation_error_LC = mps.get_truncation_error_LC();
        Eigen::Tensor<Scalar, 3> twosite_tensor(tenx::array3{dL * dR, chiL, chiR});
        if(state.get_direction() == 1) {
            // Here both M_bare are B's
            // i.e. mpsL.get_M() = GB * LB
            // and  mpsR.get_M() = GB * GB
            // So we have to attach LC from the left
            twosite_tensor.device(tenx::omp::getDevice()) = tenx::asDiagonal(LC)
                                                                .contract(mpsL.get_M(), tenx::idx({1}, {1}))
                                                                .contract(mpsR.get_M(), tenx::idx({2}, {1}))
                                                                .shuffle(tenx::array4{1, 2, 0, 3})
                                                                .reshape(tenx::array3{dL * dR, chiL, chiR});
        } else {
            // Here both M_bare are A's
            // The right A should be the previous position, so it has an attached
            // LC if we ask for get_M(), i.e. mpsR.get_M() = LA * GA * LC
            // The left A should be a simple A, i.e. mpsL.get_M() = LA * GA

            twosite_tensor.device(tenx::omp::getDevice()) =
                mpsL.get_M().contract(mpsR.get_M(), tenx::idx({2}, {1})).shuffle(tenx::array4{0, 2, 1, 3}).reshape(tenx::array3{dL * dR, chiL, chiR});
        }
        tools::finite::mps::merge_multisite_tensor(state, twosite_tensor, {static_cast<size_t>(posL), static_cast<size_t>(posR)}, static_cast<long>(posL),
                                                   chi_lim, svd_settings, LogPolicy::QUIET);
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

size_t tools::finite::mps::move_center_point_to_pos(StateFinite &state, long pos, long chi_lim, std::optional<svd::settings> svd_settings) {
    if((state.get_direction() < 0 and pos > state.get_position<long>()) or //
       (state.get_direction() > 0 and pos < state.get_position<long>()))   //
        state.flip_direction();                                            // Turn direction towards new position

    size_t moves = 0;
    while(not state.position_is_at(pos)) moves += move_center_point_single_site(state, chi_lim, svd_settings);
    return moves;
}

size_t tools::finite::mps::move_center_point_to_pos_dir(StateFinite &state, long pos, int dir, long chi_lim, std::optional<svd::settings> svd_settings) {
    if((state.get_direction() < 0 and pos > state.get_position<long>()) or //
       (state.get_direction() > 0 and pos < state.get_position<long>()))   //
        state.flip_direction();                                            // Turn direction towards new position
    size_t moves = 0;
    while(not state.position_is_at(pos, dir)) moves += move_center_point_single_site(state, chi_lim, svd_settings);
    return moves;
}

size_t tools::finite::mps::move_center_point_to_edge(StateFinite &state, long chi_lim, std::optional<svd::settings> svd_settings) {
    size_t moves = 0;
    while(not state.position_is_inward_edge()) moves += move_center_point_single_site(state, chi_lim, svd_settings);
    return moves;
}

size_t tools::finite::mps::move_center_point_to_middle(StateFinite &state, long chi_lim, std::optional<svd::settings> svd_settings) {
    size_t moves = 0;
    while(not state.position_is_the_middle_any_direction()) moves += move_center_point_single_site(state, chi_lim, svd_settings);
    return moves;
}

size_t tools::finite::mps::merge_multisite_tensor(StateFinite &state, const Eigen::Tensor<Scalar, 3> &multisite_tensor, const std::vector<size_t> &positions,
                                                  long center_position, long chi_lim, std::optional<svd::settings> svd_settings,
                                                  std::optional<LogPolicy> logPolicy) {
    auto   t_merge          = tid::tic_scope("merge");
    auto   current_position = state.get_position<long>();
    size_t moves            = static_cast<size_t>(std::abs(center_position - current_position));
    if constexpr(settings::debug)
        if(logPolicy == LogPolicy::NORMAL)
            tools::log->trace("Merging multisite tensor for sites {} | chi limit {} | dimensions {} | center {} -> {}", positions, chi_lim,
                              multisite_tensor.dimensions(), current_position, center_position);
    // Some sanity checks
    if(multisite_tensor.dimension(1) != state.get_mps_site(positions.front()).get_chiL())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim1 {} != chiL on left-most site {}",
                                             multisite_tensor.dimension(1), state.get_mps_site(positions.front()).get_chiL(), positions.front()));

    if(multisite_tensor.dimension(2) != state.get_mps_site(positions.back()).get_chiR())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim2 {} != chiR on right-most site {}",
                                             multisite_tensor.dimension(2), state.get_mps_site(positions.back()).get_chiR(), positions.back()));
    if constexpr(settings::debug_merge) {
        // We have to allow non-normalized multisite mps! Otherwise we won't be able to make them normalized
        auto norm = tenx::VectorCast(multisite_tensor).norm();
        if(std::abs(norm - 1) > 1e-8) tools::log->debug("Multisite mps for positions {} has norm far from unity: {:.16f}", positions, norm);
    }

    // Can't merge sites too far away: we would end up with interleaved A's and B sites
    if(current_position < static_cast<long>(positions.front()) - 1 or current_position > static_cast<long>(positions.back()) + 1)
        throw std::runtime_error(fmt::format("Attempting to merge sites {} too far from the current position {}", positions, current_position));

    long              spin_prod = 1;
    std::vector<long> spin_dims;
    spin_dims.reserve(positions.size());
    for(const auto &pos : positions) {
        spin_dims.emplace_back(state.get_mps_site(pos).spin_dim());
        spin_prod *= spin_dims.back();
    }
    if(spin_prod != multisite_tensor.dimension(0))
        throw std::runtime_error(
            fmt::format("Could not merge multisite mps into state: multisite_tensor dim0 {} != spin_prod {}", multisite_tensor.dimension(0), spin_prod));

    // Hold LC if moving. This should be placed in an L-slot later
    std::optional<stash<Eigen::Tensor<Scalar, 1>>> lc_hold = std::nullopt;

    if(current_position >= 0 and current_position != center_position) {
        auto &mps      = state.get_mps_site(current_position); // Guaranteed to have LC since that is the definition of current_position
        auto  pos_back = static_cast<long>(positions.back());
        auto  pos_frnt = static_cast<long>(positions.front());
        auto  pos_curr = static_cast<size_t>(current_position);
        // Detect right-move
        if(current_position < center_position and positions.size() == 1) { // This AC will become an A (AC moves to the right
            if(center_position < pos_frnt or center_position > pos_back)   // Make sure the jump isn't too long
                throw std::logic_error(fmt::format("Cannot right-move from position {} into a center position {} that is not in positions {}", current_position,
                                                   center_position, positions));
            lc_hold = stash<Eigen::Tensor<Scalar, 1>>{mps.get_LC(), mps.get_truncation_error_LC(), positions.front()};
        }
        // Detect left-move
        if(current_position > center_position and positions.size() == 1) {        // This AC position will become a B (AC moves to the left)
            if(center_position < pos_frnt - 1 or current_position > pos_back + 1) // Make sure the jump isn't too long
                throw std::logic_error(fmt::format("Cannot right-move from position {} to a center position {} in a non-neighboring group of positions {}",
                                                   current_position, center_position, positions));
            lc_hold = stash<Eigen::Tensor<Scalar, 1>>{mps.get_LC(), mps.get_truncation_error_LC(), pos_curr};
        }
        // Note that one of the positions on the split may contain a new center, so we need to unset
        // the center in our current state so we don't get duplicate centers
        mps.unset_LC();
    }

    if(svd_settings) tools::log->trace("Splitting multisite_tensor sites {} | {}", positions, svd_settings->to_string());

    // Split the multisite mps into single-site mps objects
    auto mps_list = tools::common::split::split_mps(multisite_tensor, spin_dims, positions, center_position, chi_lim, svd_settings);

    // Sanity checks
    if(positions.size() != mps_list.size())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: number of sites mismatch: positions.size() {} != mps_list.size() {}",
                                             positions.size(), mps_list.size()));

    // In multisite mergers the LC is already where we expect it to be (i.e. on the right-most "A" matrix)
    // (fuse) Copy the split up mps components into the current state
    for(auto &mps_src : mps_list) {
        auto  pos     = mps_src.get_position();
        auto &mps_tgt = state.get_mps_site(pos);

        // inject lc_hold if there is any waiting
        if(lc_hold and pos == lc_hold->pos_dst) { mps_src.set_L(lc_hold->data, lc_hold->error); }

        mps_tgt.fuse_mps(mps_src);
        state.tag_site_normalized(pos, true); // Fused site is normalized

        // Now merge stashes for neighboring sites
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
        for(auto &pos : positions) state.get_mps_site(pos).assert_identity();
    }
    return moves;
}

bool tools::finite::mps::normalize_state(StateFinite &state, std::optional<long> chi_lim, std::optional<svd::settings> svd_settings, NormPolicy norm_policy) {
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
    if(not chi_lim) chi_lim = state.find_largest_chi();
    // Save the current position, direction and center status
    auto dir   = state.get_direction();
    auto pos   = state.get_position<long>();
    auto cnt   = pos >= 0;
    auto steps = 0;
    if(tools::log->level() == spdlog::level::trace)
        tools::log->trace("Normalizing state | Old norm = {:.16f} | pos {} | dir {} | chi_lim {} | bond dims {}", tools::finite::measure::norm(state), pos, dir,
                          chi_lim.value(), tools::finite::measure::bond_dimensions(state));

    // Start with SVD at the current center position
    // NOTE: You have thought that this is unnecessary and removed it, only to find bugs much later.
    //       In particular, the bond dimension will shrink too much when doing projections, if this step is skipped.
    //       This makes sure chiL and chiR differ at most by factor spin_dim when we start the normalization
    if(pos >= 0) {
        auto &mps = state.get_mps_site(pos);
        // Make sure that the bond dimension does not increase faster than spin_dim per site
        long chi_new = std::min(chi_lim.value(), mps.spin_dim() * std::min(mps.get_chiL(), mps.get_chiR()));
        tools::finite::mps::merge_multisite_tensor(state, mps.get_M(), {static_cast<size_t>(pos)}, pos, chi_new, svd_settings, LogPolicy::QUIET);
        if constexpr(settings::debug) mps.assert_identity();
    }
    // Now we can move around the chain until we return to the original status
    while(steps++ < 2 or not state.position_is_at(pos, dir, cnt)) move_center_point_single_site(state, chi_lim.value(), svd_settings);
    state.clear_measurements();
    state.clear_cache();
    auto norm = tools::finite::measure::norm(state);
    if(tools::log->level() == spdlog::level::trace)
        tools::log->trace("Normalized  state | New norm = {:.16f} | pos {} | dir {} | chi_lim {} | bond dims {}", norm, pos, dir, chi_lim.value(),
                          tools::finite::measure::bond_dimensions(state));
    if(std::abs(norm - 1) > settings::precision::max_norm_error) {
        for(const auto &mps : state.mps_sites) {
            tools::log->warn("L ({}) | norm {:.16f} \n {}", mps->get_position(), tenx::VectorMap(mps->get_L()).norm(), mps->get_L());
            if(mps->isCenter()) tools::log->warn("LC({}) | norm {:.16f} \n {}", mps->get_position(), tenx::VectorMap(mps->get_LC()).norm(), mps->get_LC());
            mps->assert_identity();
        }
        throw std::runtime_error(fmt::format("Norm too far from unity: {:.16f} | max allowed norm error {}", norm, settings::precision::max_norm_error));
    }
    state.assert_validity();
    return true;
}

void tools::finite::mps::randomize_state(StateFinite &state, StateInit init, StateInitType type, std::string_view sector, long chi_lim, bool use_eigenspinors,
                                         std::optional<long> bitfield) {
    switch(init) {
        case StateInit::RANDOM_PRODUCT_STATE: return init::random_product_state(state, type, sector, use_eigenspinors, bitfield);
        case StateInit::RANDOM_ENTANGLED_STATE: return init::random_entangled_state(state, type, sector, chi_lim, use_eigenspinors);
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

void tools::finite::mps::truncate_all_sites(StateFinite &state, long chi_lim, std::optional<svd::settings> svd_settings) {
    tools::log->trace("Truncating all sites to bond dimension {}", chi_lim);

    auto original_position  = state.get_position();
    auto original_direction = state.get_direction();
    // Start by truncating at the current position.
    while(true) {
        move_center_point(state, chi_lim, svd_settings);
        if(state.get_position() == original_position and state.get_direction() == original_direction) {
            // Check if all bond dimensions less than or equal to below chi_lim
            auto bond_dimensions = tools::finite::measure::bond_dimensions(state);
            if(std::all_of(bond_dimensions.begin(), bond_dimensions.end(), [chi_lim](const long &chi) { return chi <= chi_lim; })) break;
        }
    }
    state.clear_cache();
    state.clear_measurements();
    tools::log->trace("Truncated all sites");
    tools::log->warn("MUST REBUILD EDGES AFTER TRUNCATING ALL SITES");
}

void tools::finite::mps::truncate_active_sites([[maybe_unused]] StateFinite &state, [[maybe_unused]] long chi_lim,
                                               [[maybe_unused]] std::optional<svd::settings> svd_settings) {
    tools::log->warn("Truncate active sites needs an implementation");
    throw std::runtime_error("Truncate active sites needs an implementation");
}

void tools::finite::mps::truncate_next_sites([[maybe_unused]] StateFinite &state, [[maybe_unused]] long chi_lim, [[maybe_unused]] size_t num_sites,
                                             [[maybe_unused]] std::optional<svd::settings> svd_settings) {
    tools::log->warn("Truncate next sites needs an implementation");
    throw std::runtime_error("Truncate next sites needs an implementation");
}

void tools::finite::mps::apply_gates(StateFinite &state, const std::vector<Eigen::Tensor<Scalar, 2>> &nsite_tensors, size_t gate_size, bool reverse,
                                     long chi_lim, std::optional<svd::settings> svd_settings) {
    // Pack the two-site operators into a vector of qm::Gates
    std::vector<qm::Gate> gates;
    gates.reserve(nsite_tensors.size());
    for(const auto &[idx, op] : iter::enumerate(nsite_tensors)) {
        auto pos = num::range<size_t>(idx, idx + gate_size, 1);
        auto dim = std::vector<long>(pos.size(), 2);
        gates.emplace_back(qm::Gate(nsite_tensors[idx], pos, dim));
    }
    apply_gates(state, gates, reverse, chi_lim, svd_settings);
}

std::vector<size_t> generate_pos_sequence(const StateFinite &state, const std::vector<qm::Gate> &gates, bool reverse) {
    // Generate a list of staggered indices
    // If 2-site gates,
    //      * Apply gates on [0-1],[2-3]... and then [1-2],[3-4]..., i.e. even sites first, then odd sites,
    //      * The corresponing list is [0,2,3,4,6....1,3,5,7,9...]
    // If 3-site gates,
    //      * Apply gates on [0-1-2], [3-4-5]... then on [1-2-3], [4-5-6]..., then on [2-3-4],[5-6-7], and so on.
    //      * The corresponing list is [0,3,6,9...1,4,7,10...2,5,8,11...]
    // So the list contains the index to the "first" or left-most" leg of the unitary.
    // When applying the inverse operation, all the indices are reversed.
    //
    // Performance note:
    // If the state is at position L-1, and the list generated has to start from 0, then L-1 moves have
    // to be done before even starting. Additionally, if unlucky, we have to move L-1 times again to return
    // to the original position.
    //    if(state.get_direction() < 0 and past_middle) state.flip_direction(); // Turn direction away from middle
    //    if(state.get_direction() > 0 and not past_middle) state.flip_direction(); // Turn direction away from middle

    auto                             state_len = state.get_length<long>();
    auto                             state_pos = state.get_position<long>();
    std::vector<std::vector<size_t>> layers;
    std::vector<std::vector<size_t>> pos_list(gates.size());
    for(const auto &[i, g] : iter::enumerate(gates)) pos_list[i] = g.pos;
    while(not pos_list.empty()) {
        // Find a gate where pos[0] == at
        std::vector<size_t> layer;
        std::vector<size_t> used;
        size_t              at = pos_list.front().front();
        for(auto &&[i, pos] : iter::enumerate(pos_list)) {
            if(at == pos.front()) {
                layer.emplace_back(at);
                used.emplace_back(i);
                at = pos.back() + 1;
            }
        }
        // Append the layer
        layers.emplace_back(layer);

        // Remove the used elements from gates_pos
        for(const auto &i : iter::reverse(used)) pos_list.erase(pos_list.begin() + static_cast<long>(i));
    }

    // Reverse every other layer to get a zigzag pattern
    // If the state is past the middle, reverse layers 0,2,4... otherwise 1,3,5...
    size_t past_middle = state_pos > state_len / 2 ? 0 : 1;
    for(auto &&[i, l] : iter::enumerate(layers)) {
        if(num::mod<size_t>(i, 2) == past_middle) std::reverse(l.begin(), l.end());
    }
    // Move the layers into a long sequence of positions
    std::vector<size_t> pos_sequence;
    for(auto &l : layers) pos_sequence.insert(pos_sequence.end(), std::make_move_iterator(l.begin()), std::make_move_iterator(l.end()));

    // To apply inverse layers we reverse the whole sequence
    if(reverse) std::reverse(pos_sequence.begin(), pos_sequence.end());
    return pos_sequence;
}

void tools::finite::mps::apply_gates(StateFinite &state, const std::vector<qm::Gate> &gates, bool reverse, long chi_lim,
                                     std::optional<svd::settings> svd_settings) {
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "  [", "]");
    if constexpr(settings::debug_gates) {
        if(tools::log->level() == spdlog::level::trace and state.get_length() <= 6) {
            auto t_dbg = tid::tic_scope("debug");
            //            tools::log->trace("Before applying gates");
            //            for(const auto &mps : state.mps_sites)
            //                std::cout << "M(" << mps->get_position() << ") dims [" << mps->spin_dim() << "," << mps->get_chiL() << "," << mps->get_chiR() <<
            //                "]:\n"
            //                          << tenx::MatrixMap(mps->get_M_bare(), mps->spin_dim(), mps->get_chiL() * mps->get_chiR()).format(CleanFmt) << std::endl;
        }
    }

    if(gates.empty()) return;
    auto pos_sequence = generate_pos_sequence(state, gates, reverse);
    if constexpr(settings::debug_gates)
        tools::log->trace("current pos {} dir {} | pos_sequence {}", state.get_position<long>(), state.get_direction(), pos_sequence);

    state.clear_cache(LogPolicy::QUIET);
    Eigen::Tensor<Scalar, 3> gate_mps;
    for(const auto &[idx, pos] : iter::enumerate(pos_sequence)) {
        auto &gate = gates[pos];
        if(gate.pos.back() >= state.get_length()) throw std::logic_error(fmt::format("The last position of gate {} is out of bounds: {}", pos, gate.pos));

        move_center_point_to_pos(state, static_cast<long>(gate.pos.front()), chi_lim, svd_settings);
        auto multisite_mps = state.get_multisite_mps(gate.pos);
        gate_mps.resize(std::array<long, 3>{gate.op.dimension(0), multisite_mps.dimension(1), multisite_mps.dimension(2)});
        if(reverse)
            gate_mps.device(tenx::omp::getDevice()) = gate.adjoint().contract(multisite_mps, tenx::idx({0}, {0}));
        else
            gate_mps.device(tenx::omp::getDevice()) = gate.op.contract(multisite_mps, tenx::idx({0}, {0}));
        long min_position = static_cast<long>(gate.pos.front()) - 1;
        long max_position = static_cast<long>(gate.pos.back());
        long tgt_position = static_cast<long>(pos_sequence[std::min<size_t>(idx + 1, pos_sequence.size() - 1)]);
        long new_position = std::clamp<long>(tgt_position, min_position, max_position);
        tools::log->trace("Merging gate sites {} dims {}", gate.pos, multisite_mps.dimensions());
        if constexpr(settings::debug_gates)
            tools::log->trace("pos {} | tgt {} | new {} | from {} - {} | labels {}", gate.pos, tgt_position, new_position, state.get_position<long>(),
                              new_position, state.get_labels());
        tools::finite::mps::merge_multisite_tensor(state, gate_mps, gate.pos, new_position, chi_lim, svd_settings, LogPolicy::NORMAL);
    }

    if constexpr(settings::debug_gates) {
        if(tools::log->level() == spdlog::level::trace and state.get_length() <= 6) {
            auto t_dbg = tid::tic_scope("debug");
            //            tools::log->trace("After applying gates");
            //            for(const auto &mps : state.mps_sites)
            //                std::cout << "M(" << mps->get_position() << ") dims [" << mps->spin_dim() << "," << mps->get_chiL() << "," << mps->get_chiR() <<
            //                "]:\n"
            //                          << tenx::MatrixMap(mps->get_M_bare(), mps->spin_dim(), mps->get_chiL() * mps->get_chiR()).format(CleanFmt) << std::endl;
        }
    }

    move_center_point_to_edge(state, chi_lim, svd_settings);

    auto has_normalized = tools::finite::mps::normalize_state(state, chi_lim, svd_settings, NormPolicy::IFNEEDED);
    if constexpr(settings::debug_gates)
        if(has_normalized and tools::log->level() == spdlog::level::trace and state.get_length() <= 6) {
            auto t_dbg = tid::tic_scope("debug");
            //            tools::log->trace("After normalization");
            //            for(const auto &mps : state.mps_sites)
            //                std::cout << "M(" << mps->get_position() << ") dims [" << mps->spin_dim() << "," << mps->get_chiL() << "," << mps->get_chiR() <<
            //                "]:\n"
            //                          << tenx::MatrixMap(mps->get_M_bare(), mps->spin_dim(), mps->get_chiL() * mps->get_chiR()).format(CleanFmt) << std::endl;
        }

    if(settings::profiling::extra) {
        static long long num_old = 0;
        tools::log->info(
            "Finished step: move {:.4f} | apply {:.4f} | merge {:.4f} | split {:.4f} | merge {:.4f} | svdm {:.4f} | svda {:.4f} | svdb {:.4f} | "
            "svdA {:.4f} | svdB {:.4f} | wrk {:.4f} | adj {:.4f} | jac {:.4f} | svd {:.4f} | num {} (+{})",
            1000 * tid::get("ANY::t_gate_move").get_time(), 1000 * tid::get("ANY::t_gate_apply").get_time(), 1000 * tid::get("ANY::t_gate_merge").get_time(),
            1000 * tid::get("ANY::t_merge_split").get_time(), 1000 * tid::get("ANY::t_merge_merge").get_time(), 1000 * tid::get("ANY::t_split_svdm").get_time(),
            1000 * tid::get("ANY::t_split_svda").get_time(), 1000 * tid::get("ANY::t_split_svdb").get_time(), 1000 * tid::get("ANY::t_splitA_svd").get_time(),
            1000 * tid::get("ANY::t_splitB_svd").get_time(), 1000 * tid::get("ANY::t_svd_wrk").get_time(), 1000 * tid::get("ANY::t_svd_adj").get_time(),
            1000 * tid::get("ANY::t_svd_jac").get_time(), 1000 * tid::get("ANY::t_svd_svd").get_time(), svd::solver::count.value(),
            svd::solver::count.value() - num_old);
        tid::get("ANY::t_gate_move").reset();
        tid::get("ANY::t_gate_apply").reset();
        tid::get("ANY::t_gate_merge").reset();
        tid::get("ANY::t_merge_split").reset();
        tid::get("ANY::t_merge_merge").reset();
        tid::get("ANY::t_split_svdm").reset();
        tid::get("ANY::t_split_svda").reset();
        tid::get("ANY::t_split_svdb").reset();
        tid::get("ANY::t_splitA_svd").reset();
        tid::get("ANY::t_splitB_svd").reset();
        tid::get("ANY::t_svd_wrk").reset();
        tid::get("ANY::t_svd_adj").reset();
        tid::get("ANY::t_svd_jac").reset();
        tid::get("ANY::t_svd_svd").reset();
        num_old = svd::solver::count.value();
    }
}
