//
// Created by david on 2019-01-29.
//

#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>
// -- (textra first)
#include <config/enums.h>
#include <config/nmspc_settings.h>
#include <general/nmspc_iter.h>
#include <iostream>
#include <math/num.h>
#include <math/svd.h>
#include <physics/nmspc_quantum_mechanics.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/common/split.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <utility>



bool tools::finite::mps::internal::bitfield_is_valid(std::optional<long> bitfield) {
    return bitfield.has_value() and bitfield.value() > 0 and internal::used_bitfields.count(bitfield.value()) == 0;
}




    if(state.position_is_outward_edge()){
        if(state.get_direction() == -1 and state.get_mps_site(0l).get_chiL() != 1 )
            throw std::logic_error(fmt::format("chiL at position 0 must have dimension 1, but it has dimension {}. Mps dims {}",
                                                state.get_mps_site(0l).get_chiL(), state.get_mps_site(0l).dimensions()));
        if(state.get_direction() == 1 and state.get_mps_site().get_chiR() != 1 )
            throw std::logic_error(fmt::format("chiR at position {} must have dimension 1, but it has dimension {}. Mps dims {}",
                                               state.get_position(), state.get_mps_site().get_chiR(), state.get_mps_site().dimensions()));
        return state.flip_direction();  // Instead of moving out of the chain, just flip the direction and return
    }else
    {
        long   pos            = state.get_position<long>(); // If all sites are B's, then this is -1. Otherwise this is the current "A*LC" site
        long   posC           = pos + state.get_direction(); // This is the site which becomes the new center position
        if(pos  < -1 or pos >= state.get_length<long>()) throw std::runtime_error(fmt::format("pos out of bounds: {}",pos));
        if(posC < -1 or posC >= state.get_length<long>()) throw std::runtime_error(fmt::format("posC out of bounds: {}",posC));
        if(state.get_direction() != posC - pos) throw std::logic_error(fmt::format("Expected posC - pos == {}. Got {}",state.get_direction(), posC-pos));
//        auto & mps            = state.get_mps_site(pos);
//        if(not mps.isCenter()) posC = pos; // In this case we need to create a new center at position 0 since all sites are B's

        Eigen::Tensor<Scalar, 1> LC(1); LC.setConstant(1); // Store the LC bond in a temporary. It will become a regular "L" bond later
        double                   truncation_error_LC = 0; // Same story with the truncation error
        if(pos >= 0){
            auto & mps          = state.get_mps_site(pos);
            LC                  = mps.get_LC();
            truncation_error_LC = mps.get_truncation_error_LC();
        }
        if(state.get_direction() == 1){
            auto & mpsC = state.get_mps_site(posC); //This becomes the new center position AC
            Eigen::Tensor<Scalar, 3> onesite_tensor(mpsC.dimensions()); // Allocate for contraction of LC * B
            onesite_tensor.device(Textra::omp::getDevice())  = Textra::asDiagonal(LC).contract(mpsC.get_M(), Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
            tools::finite::mps::merge_multisite_tensor(state, onesite_tensor, {static_cast<size_t>(posC)}, posC, chi_lim, svd_threshold, LogPolicy::QUIET);
            mpsC.set_L(LC, truncation_error_LC); // Copy old "LC" into the "L" slot of the new "A" at position "posC"
        }else if (state.get_direction() == -1) {
            auto & mps = state.get_mps_site(pos); //This becomes the new B
            auto &onesite_tensor = mps.get_M(); // No need to contract anything this time.
            tools::finite::mps::merge_multisite_tensor(state, onesite_tensor, {static_cast<size_t>(pos)}, posC, chi_lim, svd_threshold, LogPolicy::QUIET);
            mps.set_L(LC, truncation_error_LC); // Copy old "LC" into the "L" slot of the new "B" at position "pos"}
        }
        state.clear_cache(LogPolicy::QUIET);
        state.clear_measurements(LogPolicy::QUIET);
    }
}


size_t tools::finite::mps::move_center_point(class_state_finite &state, long chi_lim, std::optional<double> svd_threshold) {
    if(state.position_is_outward_edge(2)) {
        state.flip_direction(); // Instead of moving out of the chain, just flip the direction and return
        return 0; // No moves this time, return 0
    }
    else {
        size_t pos  = state.get_position();
        size_t posL = state.get_direction() == 1 ? pos + 1 : pos - 1;
        size_t posR = state.get_direction() == 1 ? pos + 2 : pos;
        auto & mps  = state.get_mps_site();
        auto & mpsL = state.get_mps_site(posL); // Becomes the new center position
        auto & mpsR = state.get_mps_site(posR); // The site to the right of the new center position
        long   dL   = mpsL.spin_dim();
        long   dR   = mpsR.spin_dim();
        long   chiL = mpsL.get_chiL();
        long   chiR = mpsR.get_chiR();
        // Store the special LC bond in a temporary. It needs to be put back afterwards
        // Do the same with its truncation error
        Eigen::Tensor<Scalar, 1> LC                  = mps.get_LC();
        double                   truncation_error_LC = mps.get_truncation_error_LC();
        Eigen::Tensor<Scalar, 3> twosite_tensor(Textra::array3{dL * dR, chiL, chiR});
        if(state.get_direction() == 1) {
            // Here both M_bare are B's
            // i.e. mpsL.get_M() = GB * LB
            // and  mpsR.get_M() = GB * GB
            // So we have to attach LC from the left
            twosite_tensor.device(Textra::omp::getDevice()) = Textra::asDiagonal(LC)
                                                                  .contract(mpsL.get_M(), Textra::idx({1}, {1}))
                                                                  .contract(mpsR.get_M(), Textra::idx({2}, {1}))
                                                                  .shuffle(Textra::array4{1, 2, 0, 3})
                                                                  .reshape(Textra::array3{dL * dR, chiL, chiR});
        } else {
            // Here both M_bare are A's
            // The right A should be the previous position, so it has an attached
            // LC if we ask for get_M(), i.e. mpsR.get_M() = LA * GA * LC
            // The left A should be a simple A, i.e. mpsL.get_M() = LA * GA

            twosite_tensor.device(Textra::omp::getDevice()) =
                mpsL.get_M().contract(mpsR.get_M(), Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{dL * dR, chiL, chiR});
        }
        tools::finite::mps::merge_multisite_tensor(state, twosite_tensor, {posL, posR}, static_cast<long>(posL), chi_lim, svd_threshold, LogPolicy::QUIET);
        state.clear_cache(LogPolicy::QUIET);
        state.clear_measurements(LogPolicy::QUIET);

        // Put LC where it belongs.
        // Recall that mpsL, mpsR are on the new position, not the old one!
        if(state.get_direction() == 1) mpsL.set_L(LC, truncation_error_LC);
        else
            mpsR.set_L(LC, truncation_error_LC);
    }
}

void tools::finite::mps::move_center_point_to_edge(class_state_finite &state, long chi_lim, std::optional<double> svd_threshold) {
    while(not state.position_is_any_edge()) move_center_point_single_site(state, chi_lim, svd_threshold);
}

void tools::finite::mps::move_center_point_to_middle(class_state_finite &state, long chi_lim, std::optional<double> svd_threshold) {
    while(not state.position_is_the_middle()) move_center_point_single_site(state, chi_lim, svd_threshold);
}

void tools::finite::mps::merge_multisite_tensor(class_state_finite &state, const Eigen::Tensor<Scalar, 3> &multisite_mps, const std::vector<size_t> &sites,
                                                long center_position, long chi_lim, std::optional<double> svd_threshold, std::optional<LogPolicy> logPolicy) {
    std::optional<long> current_position = state.get_position<long>();

    if(logPolicy == LogPolicy::NORMAL) tools::log->trace("Merging multisite tensor for sites {} | chi limit {} | dimensions {} | center {} -> {}",
                          sites, chi_lim, multisite_mps.dimensions(),current_position.value(), center_position);
    // Some sanity checks
    if(multisite_mps.dimension(1) != state.get_mps_site(sites.front()).get_chiL())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim1 {} != chiL on left-most site {}", multisite_mps.dimension(1),
                                             state.get_mps_site(sites.front()).get_chiL(), sites.front()));

    if(multisite_mps.dimension(2) != state.get_mps_site(sites.back()).get_chiR())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim2 {} != chiR on right-most site {}", multisite_mps.dimension(2),
                                             state.get_mps_site(sites.back()).get_chiR(), sites.back()));
    if constexpr(settings::debug){
        auto norm = Textra::Tensor_to_Vector(multisite_mps).norm();
        if(std::abs(norm-1) > 1e-8) throw std::runtime_error(fmt::format("Multisite mps norm is too far from unity: {:.16f}",norm));
    }

    long              spin_prod = 1;
    std::vector<long> spin_dims;
    spin_dims.reserve(sites.size());
    for(const auto &site : sites) {
        spin_dims.emplace_back(state.get_mps_site(site).spin_dim());
        spin_prod *= spin_dims.back();
    }
    if(spin_prod != multisite_mps.dimension(0))
        throw std::runtime_error(
            fmt::format("Could not merge multisite mps into state: multisite_mps dim0 {} != spin_prod {}", multisite_mps.dimension(0), spin_prod));

    // Split the multisite mps into single-site mps objects
    tools::common::profile::prof[AlgorithmType::ANY]["t_merge_split"]->tic();
    auto mps_list = tools::common::split::split_mps(multisite_mps, spin_dims, sites, center_position, chi_lim, svd_threshold);
    tools::common::profile::prof[AlgorithmType::ANY]["t_merge_split"]->toc();

    // Sanity checks
    if(sites.size() != mps_list.size())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: number of sites mismatch: positions.size() {} != mps_list.size() {}",
                                             sites.size(), mps_list.size()));

    // Note that one of the positions on the split will be a center, so we need to unset
    // the center in our current state so we don't get duplicate centers
    if(current_position.value() != center_position) {
        state.get_mps_site(current_position.value()).unset_LC();
        current_position = std::nullopt;
    }

    // In multisite mergers the LC is already where we expect it to be (i.e. on the right-most "A" matrix)
    // Copy the split up mps components into the current state
    tools::common::profile::prof[AlgorithmType::ANY]["t_merge_merge"]->tic();
    for(const auto &mps_src : mps_list) {
        auto pos = mps_src.get_position();
        auto &mps_tgt = state.get_mps_site(pos);
        mps_tgt.merge_mps(mps_src);
        state.tag_site_normalized(pos, true); // Merged site is normalized
        if(mps_src.has_stash_V()){
            if(pos == state.get_length()-1) mps_src.unstash(); // Discard whatever is stashed at the edge (this normalizes the state)
            else{
                // Extract the stashed V-matrix and multiply onto the "B" on the right
                if(logPolicy == LogPolicy::NORMAL) tools::log->trace("Merging stash from site {} into {}",pos,pos+1);
                state.get_mps_site(pos+1).merge_stash(mps_src); // Absorb remaining V
            }
        }
        if(mps_src.has_stash_U()){
            if(pos == 0) mps_src.unstash();  // Discard whatever is stashed at the edge (this normalizes the state)
            else{
                if(logPolicy == LogPolicy::NORMAL) tools::log->trace("Merging stash from site {} into {}",pos,pos-1);
                state.get_mps_site(pos-1).merge_stash(mps_src); // Absorb remaining U,S=LC
            }
        }
    }
    tools::common::profile::prof[AlgorithmType::ANY]["t_merge_merge"]->toc();
    if(not current_position) current_position = state.get_position<long>();
    if(current_position.value() != center_position) throw std::logic_error(fmt::format("Center position mismatch {} ! {}\nLabels: {}", current_position.value(),center_position,state.get_labels()));
    state.clear_cache(LogPolicy::QUIET);
    state.clear_measurements(LogPolicy::QUIET);
}

//bool tools::finite::mps::normalize_state(class_state_finite &state, long chi_lim, std::optional<double> svd_threshold, NormPolicy norm_policy) {
//    // When a state needs to be normalized it's enough to "move" the center position around the whole chain.
//    // Each move performs an SVD decomposition which leaves unitaries behind, effectively normalizing the state.
//    // NOTE! It may be important to start with the current position.
//
//    if(norm_policy == NormPolicy::IFNEEDED) {
//        // We may only go ahead with a normalization if its really needed.
//        if(std::abs(tools::finite::measure::norm(state) - 1.0) < settings::precision::max_norm_error) return false;
//    }
//
//    // Otherwise we just do the normalization
//    if(tools::Logger::getLogLevel(tools::log) <= 0) tools::log->trace("Normalizing state | Old norm = {:.16f}", tools::finite::measure::norm(state));
//
//    // Start with normalizing at the current position
//    size_t num_moves = 2 * (state.get_length() - 1);
//    size_t posL      = state.get_position();
//    size_t posR      = posL + 1;
//    auto & mpsL      = state.get_mps_site(posL);
//    auto & mpsR      = state.get_mps_site(posR);
//    long   dL        = mpsL.spin_dim();
//    long   dR        = mpsR.spin_dim();
//    long   chiL      = mpsL.get_chiL();
//    long   chiR      = mpsR.get_chiR();
//
//    Eigen::Tensor<Scalar, 3> twosite_tensor(Textra::array3{dL * dR, chiL, chiR});
//    twosite_tensor.device(Textra::omp::getDevice()) =
//        mpsL.get_M().contract(mpsR.get_M(), Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{dL * dR, chiL, chiR});
//    tools::finite::mps::merge_multisite_tensor(state, twosite_tensor, {posL, posR}, posL, chi_lim, svd_threshold, LogPolicy::QUIET);
//
//    // Now we can move around the chain
//    for(size_t move = 0; move < num_moves; move++) move_center_point(state, chi_lim, svd_threshold);
//    state.clear_measurements();
//    state.clear_cache();
//    if(tools::Logger::getLogLevel(tools::log) <= 0) tools::log->debug("Normalized state | New norm = {:.16f}", tools::finite::measure::norm(state));
//    state.assert_validity();
//    return true;
//}

bool tools::finite::mps::normalize_state(class_state_finite &state, long chi_lim, std::optional<double> svd_threshold, NormPolicy norm_policy) {
    // When a state needs to be normalized it's enough to "move" the center position around the whole chain.
    // Each move performs an SVD decomposition which leaves unitaries behind, effectively normalizing the state.
    // NOTE! It may be important to start with the current position.

    if(norm_policy == NormPolicy::IFNEEDED) {
        // We may only go ahead with a normalization if its really needed.
        auto norm = tools::finite::measure::norm(state);
        tools::log->trace("Norm: {:.16f}",norm);
        if(std::abs(norm - 1.0) < settings::precision::max_norm_error) return false;
    }
    // Otherwise we just do the normalization

    // Save the current position, direction and center status
    auto dir = state.get_direction();
    auto pos = state.get_position<long>();
    auto cnt = state.get_mps_site().isCenter();
    if(tools::Logger::getLogLevel(tools::log) <= 0) tools::log->trace("Normalizing state | Old norm = {:.16f} | pos {} | dir {}", tools::finite::measure::norm(state),pos,dir);

    // Start with SVD at the current position
    auto & mps      = state.get_mps_site(pos);
    tools::finite::mps::merge_multisite_tensor(state, mps.get_M(), {static_cast<size_t>(pos)}, pos, chi_lim, svd_threshold, LogPolicy::QUIET);
    move_center_point_single_site(state, chi_lim, svd_threshold); // Move once to get started
    // Now we can move around the chain until we return to the original status
    while(not state.position_is_at(pos,dir,cnt))
        move_center_point_single_site(state, chi_lim, svd_threshold);
    state.clear_measurements();
    state.clear_cache();
    if(tools::Logger::getLogLevel(tools::log) <= 0) tools::log->debug("Normalized state | New norm = {:.16f}", tools::finite::measure::norm(state));
    state.assert_validity();
    return true;
}



void tools::finite::mps::randomize_state(class_state_finite &state, StateInit init, StateInitType type, const std::string &sector, long chi_lim,
                                         bool use_eigenspinors, std::optional<long> bitfield) {
    switch(init) {
        case StateInit::RANDOM_PRODUCT_STATE: return internal::random_product_state(state, type, sector, use_eigenspinors, bitfield);
        case StateInit::RANDOM_ENTANGLED_STATE: return internal::random_entangled_state(state, type, sector, chi_lim, use_eigenspinors);
        case StateInit::RANDOMIZE_PREVIOUS_STATE: return internal::randomize_given_state(state, type);
        case StateInit::PRODUCT_STATE_ALIGNED: return internal::set_product_state_aligned(state, type, sector);
        case StateInit::PRODUCT_STATE_NEEL: return internal::set_product_state_neel(state, type, sector);
    }
}

void tools::finite::mps::apply_random_paulis(class_state_finite &state, const std::vector<Eigen::Matrix2cd> &paulimatrices) {
    auto [mpos, L, R] = qm::mpo::sum_of_pauli_mpo(paulimatrices, state.get_length(), RandomizerMode::SELECT1);
    tools::finite::ops::apply_mpos(state, mpos, L, R);
}

void tools::finite::mps::apply_random_paulis(class_state_finite &state, const std::vector<std::string> &paulistrings) {
    std::vector<Eigen::Matrix2cd> paulimatrices;
    for(const auto &str : paulistrings) paulimatrices.emplace_back(internal::get_pauli(str));
    apply_random_paulis(state, paulimatrices);
}

void tools::finite::mps::truncate_all_sites(class_state_finite &state, long chi_lim, std::optional<double> svd_threshold) {
    tools::log->trace("Truncating all sites to bond dimension {}", chi_lim);

    auto original_position  = state.get_position();
    auto original_direction = state.get_direction();
    // Start by truncating at the current position.
    while(true) {
        move_center_point(state, chi_lim, svd_threshold);
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

void tools::finite::mps::truncate_active_sites([[maybe_unused]] class_state_finite &state, [[maybe_unused]] long chi_lim,
                                               [[maybe_unused]] std::optional<double> svd_threshold) {
    tools::log->warn("Truncate active sites needs an implementation");
    throw std::runtime_error("Truncate active sites needs an implementation");
}

void tools::finite::mps::truncate_next_sites([[maybe_unused]] class_state_finite &state, [[maybe_unused]] long chi_lim, [[maybe_unused]] size_t num_sites,
                                             [[maybe_unused]] std::optional<double> svd_threshold) {
    tools::log->warn("Truncate next sites needs an implementation");
    throw std::runtime_error("Truncate next sites needs an implementation");
}

void tools::finite::mps::apply_gates(class_state_finite &state, const std::vector<Eigen::Tensor<Scalar, 2>> &nsite_tensors, size_t gate_size, bool reverse,
                                             long chi_lim, std::optional<double> svd_threshold) {
    // Pack the two-site operators into a vector of UnitaryGates
    std::vector<qm::Gate> gates;
    gates.reserve(nsite_tensors.size());
    for(const auto & [idx, op] : iter::enumerate(nsite_tensors)) gates.emplace_back(qm::Gate(nsite_tensors[idx],num::range<size_t>(idx,idx+gate_size,1)));
    apply_gates(state, gates, reverse, chi_lim, svd_threshold);
}


void tools::finite::mps::apply_gates(class_state_finite &state, const std::vector<qm::Gate> &gates, bool reverse, long chi_lim,
                                             std::optional<double> svd_threshold) {

    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "  [", "]");
    if constexpr(settings::debug){
        if(tools::log->level() == spdlog::level::trace and state.get_length() <= 6){
            tools::common::profile::get_default_prof()["t_dbg"]->tic();
            tools::log->trace("Before applying gates");
            for(const auto & mps : state.mps_sites)
                std::cout << "M(" << mps->get_position() << ") dims [" << mps->spin_dim() << "," << mps->get_chiL() << "," << mps->get_chiR() << "]:\n"
                          << Textra::TensorMatrixMap(mps->get_M_bare(), mps->spin_dim(), mps->get_chiL() * mps->get_chiR()).format(CleanFmt) << std::endl;
            tools::common::profile::get_default_prof()["t_dbg"]->toc();
        }
    }



    if(gates.empty()) return;
    auto gate_size = gates.front().pos.size(); // The size of each gate, i.e. how many sites are updated by each gate.
    for (const auto & [idx,gate] : iter::enumerate(gates)) // Check that all gates are of the same size
        if(gate.pos.size() != gate_size)
            throw std::runtime_error(fmt::format("Gate size mismatch: "
                                                 "Gate 0 has {} sites: {} | "
                                                 "gate {} has {} sites: {}",
                                                 gate_size,gates.front().pos, idx,gate.pos.size(),gate.pos ));


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

    bool past_middle = state.get_position<long>() > state.get_length<long>()/2;
    if(state.get_direction() < 0 and past_middle ) state.flip_direction();
    if(state.get_direction() > 0 and not past_middle) state.flip_direction();
    size_t flip = past_middle ? 0 : 1;
    std::vector<size_t> all_idx;
    for(size_t offset = 0; offset < gate_size; offset++){
        if(offset+gate_size > state.get_length()) break;
        auto off_idx = num::range<size_t>(offset, state.get_length()-gate_size+1, gate_size);
        if(num::mod<size_t>(offset,2) == flip) std::reverse(off_idx.begin(), off_idx.end()); // If odd, reverse the sequence
        tools::log->trace("Appending idx {}", off_idx);
        all_idx.insert(all_idx.end(), off_idx.begin(), off_idx.end());
    }
    if(reverse) std::reverse(all_idx.begin(), all_idx.end());

    if constexpr (settings::debug) tools::log->trace("current pos {} dir {} | all_idx {}",state.get_position<long>(),state.get_direction(), all_idx);

    // Save current position and direction so we can leave this function in the same condition
    auto save_pos = state.get_position();
    auto save_dir = state.get_direction();


    state.clear_cache(LogPolicy::QUIET);
    Eigen::Tensor<Scalar, 3> gate_mps;
    for(const auto & idx : all_idx) {
        auto &gate = gates[idx];
        if(gate.pos != num::range<size_t>(gate.pos.front(), gate.pos.front()+gate_size, 1))// Just a sanity check that gate.pos is well defined
            throw std::runtime_error(fmt::format("The positions on gate {} are not well defined for a {}-site gate: {}", idx,gate_size,gate.pos));
        if(gate.pos.back() >= state.get_length())
            throw std::logic_error(fmt::format("The last position of gate {} is out of bounds: {}", idx, gate.pos));

        tools::common::profile::prof[AlgorithmType::ANY]["t_gate_move"]->tic();
        while(state.get_position() != gate.pos.front()) move_center_point_single_site(state, chi_lim, svd_threshold); // Move into position

        tools::common::profile::prof[AlgorithmType::ANY]["t_gate_move"]->toc();
        tools::common::profile::prof[AlgorithmType::ANY]["t_gate_apply"]->tic();
        tools::log->debug("Constructing multisite_mps at sites {}", gate.pos);
        auto multisite_mps = state.get_multisite_mps(gate.pos);
        tools::log->debug("Contracting gate and multisite_mps at sites {}", gate.pos);
        gate_mps.resize({gate.op.dimension(0), multisite_mps.dimension(1), multisite_mps.dimension(2)});
        if(reverse) gate_mps.device(Textra::omp::getDevice()) = gate.adjoint().contract(multisite_mps, Textra::idx({0}, {0}));
        else
            gate_mps.device(Textra::omp::getDevice()) = gate.op.contract(multisite_mps, Textra::idx({0}, {0}));
        tools::common::profile::prof[AlgorithmType::ANY]["t_gate_apply"]->toc();

        // Calculate new center position
//        long center_position = std::clamp(
//            state.get_position<long>() + state.get_direction() * static_cast<long>(gate.pos.size()),
//            0l, static_cast<long>(state.get_length())-1);

        tools::log->debug("Merging at sites {}", gate.pos);
        tools::common::profile::prof[AlgorithmType::ANY]["t_gate_merge"]->tic();
//        tools::finite::mps::merge_multisite_tensor(state, gate_mps, gate.pos, static_cast<long>(gate.pos.front()), chi_lim, svd_threshold, LogPolicy::QUIET);
        tools::finite::mps::merge_multisite_tensor(state, gate_mps, gate.pos, state.get_position<long>(), chi_lim, svd_threshold, LogPolicy::QUIET);
        tools::common::profile::prof[AlgorithmType::ANY]["t_gate_merge"]->toc();
        tools::log->debug("Finished step: move {:.4f} | apply {:.4f} | merge {:.4f} | svdm {:.4f} | svda {:.4f} | svdb {:.4f} | svdA {:.4f} | svdB {:.4f}",
                          1000*tools::common::profile::prof[AlgorithmType::ANY]["t_gate_move"]->get_last_interval(),
                          1000*tools::common::profile::prof[AlgorithmType::ANY]["t_gate_apply"]->get_last_interval(),
                          1000*tools::common::profile::prof[AlgorithmType::ANY]["t_gate_merge"]->get_last_interval(),
                          1000*tools::common::profile::prof[AlgorithmType::ANY]["t_split_svdm"]->get_last_interval(),
                          1000*tools::common::profile::prof[AlgorithmType::ANY]["t_split_svda"]->get_last_interval(),
                          1000*tools::common::profile::prof[AlgorithmType::ANY]["t_split_svdb"]->get_last_interval(),
                          1000*tools::common::profile::prof[AlgorithmType::ANY]["t_splitA_svd"]->get_measured_time(),
                          1000*tools::common::profile::prof[AlgorithmType::ANY]["t_splitB_svd"]->get_measured_time()
                          );
        tools::common::profile::prof[AlgorithmType::ANY]["t_splitA_svd"]->reset();
        tools::common::profile::prof[AlgorithmType::ANY]["t_splitB_svd"]->reset();
    }


    if constexpr(settings::debug){
        if(tools::log->level() == spdlog::level::trace and state.get_length() <= 6 ){
            tools::common::profile::get_default_prof()["t_dbg"]->tic();
            tools::log->trace("After applying gates");
            for(const auto & mps : state.mps_sites)
                std::cout << "M(" << mps->get_position() << ") dims [" << mps->spin_dim() << "," << mps->get_chiL() << "," << mps->get_chiR() << "]:\n"
                          << Textra::TensorMatrixMap(mps->get_M_bare(), mps->spin_dim(), mps->get_chiL() * mps->get_chiR()).format(CleanFmt) << std::endl;
            tools::common::profile::get_default_prof()["t_dbg"]->toc();
        }
    }

    tools::common::profile::prof[AlgorithmType::ANY]["t_gate_return"]->tic();
    move_center_point_to_edge(state,chi_lim,svd_threshold);
    tools::common::profile::prof[AlgorithmType::ANY]["t_gate_return"]->toc();


    auto has_normalized = tools::finite::mps::normalize_state(state, chi_lim, svd_threshold, NormPolicy::IFNEEDED);
    if constexpr(settings::debug)
        if(has_normalized and tools::log->level() == spdlog::level::trace and state.get_length() <= 6){
            tools::common::profile::get_default_prof()["t_dbg"]->tic();
            tools::log->trace("After normalization");
            for(const auto & mps : state.mps_sites)
                std::cout << "M(" << mps->get_position() << ") dims [" << mps->spin_dim() << "," << mps->get_chiL() << "," << mps->get_chiR() << "]:\n"
                          << Textra::TensorMatrixMap(mps->get_M_bare(), mps->spin_dim(), mps->get_chiL() * mps->get_chiR()).format(CleanFmt) << std::endl;
            tools::common::profile::get_default_prof()["t_dbg"]->toc();
        }


}
