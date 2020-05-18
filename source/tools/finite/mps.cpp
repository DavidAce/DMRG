//
// Created by david on 2019-01-29.
//

#include <config/enums.h>
#include <config/nmspc_settings.h>
#include <general/nmspc_quantum_mechanics.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/finite/debug.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <tools/finite/svd.h>

#include "measure.h"
#include <tools/common/prof.h>
#include <tools/common/svd.h>
#include <utility>


void tools::finite::mps::move_center_point(class_state_finite &state, std::optional<long> chi_lim, std::optional<double> svd_threshold) {
    if(state.position_is_any_edge()) {
        // Instead of moving out of the chain, just flip the direction and return
        state.flip_direction();
    } else {
        size_t      pos  = state.get_position();
        size_t      posL = state.get_direction() == 1 ? pos + 1 : pos - 1;
        size_t      posR = state.get_direction() == 1 ? pos + 2 : pos;
        auto &      mps  = state.get_mps_site();
        const auto &mpsL = state.get_mps_site(posL); // Becomes the new center position
        const auto &mpsR = state.get_mps_site(posR); // The site to the right of the new center position
        long        dL   = mpsL.spin_dim();
        long        dR   = mpsR.spin_dim();
        long        chiL = mpsL.get_chiL();
        long        chiR = mpsR.get_chiR();
        // Store the special LC bond in a temporary.
        Eigen::Tensor<Scalar, 1> LC = mps.get_LC();
        Eigen::Tensor<Scalar, 4> twosite_tensor;
        if(state.get_direction() == 1) {
            twosite_tensor = Textra::asDiagonal(LC)
                .contract(mpsL.get_M_bare(), Textra::idx({1}, {1}))
                .contract(mpsR.get_M_bare(), Textra::idx({2}, {1}))
                .shuffle(Textra::array4{1, 2, 0, 3})
                .reshape(Textra::array3{dL * dR, chiL, chiR});
        } else {
            twosite_tensor = mpsL.get_M_bare()
                .contract(mpsR.get_M_bare(), Textra::idx({2}, {1}))
                .contract(Textra::asDiagonal(LC), Textra::idx({3}, {0}))
                .shuffle(Textra::array4{0, 2, 1, 3})
                .reshape(Textra::array3{dL * dR, chiL, chiR});
        }
        tools::finite::mps::merge_multisite_tensor(state, twosite_tensor, {posL, posR}, posL, chi_lim, svd_threshold);
        state.clear_cache();
        state.clear_measurements();
    }
}


void tools::finite::mps::merge_multisite_tensor(class_state_finite &state, const Eigen::Tensor<Scalar, 3> &multisite_mps, const std::list<size_t> &positions,
                                             size_t center_position, std::optional<long> chi_lim, std::optional<double> svd_threshold) {
    // Some sanity checks
    if(multisite_mps.dimension(1) != state.get_mps_site(positions.front()).get_chiL())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim1 {} != chiL on left-most site {}", multisite_mps.dimension(1),
                                             state.get_mps_site(positions.front()).get_chiL(), positions.front()));

    if(multisite_mps.dimension(2) != state.get_mps_site(positions.back()).get_chiR())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim2 {} != chiR on right-most site {}", multisite_mps.dimension(2),
                                             state.get_mps_site(positions.back()).get_chiR(), positions.back()));
    if(not chi_lim)
        chi_lim = state.get_chi_lim();

    std::list<long> spin_dims;
    for(const auto &site : positions) spin_dims.emplace_back(state.get_mps_site(site).spin_dim());

    // Split the multisite mps into single-site mps objects
    auto mps_list = tools::common::svd::split_mps(multisite_mps, spin_dims, positions, center_position, chi_lim.value(), svd_threshold);

    if(positions.size() != mps_list.size())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: number of sites mismatch: positions.size() {} != mps_list.size() {}",
                                             positions.size(), mps_list.size()));

    // Note that one of the positions on the split will be a center, so we need to unset
    // the center in our current state so we don't get duplicate centers
    state.get_mps_site().unset_LC();

    // Copy the split up mps components into the current state
    auto mps_tgt = std::next(state.mps_sites.begin(), static_cast<long>(positions.front()));
    for(const auto &mps_src : mps_list) {
        if(mps_tgt->get_position() != mps_src.get_position())
            throw std::runtime_error(fmt::format("Could not merge multisite mps into state: Position mismatch: mps_tgt pos {} != mps_src pos {}",
                                                 mps_tgt->get_position(), mps_src.get_position()));
        mps_tgt->set_M(mps_src.get_M());
        if(mps_src.get_L().size() > 0) // The edges have empty "L"
            mps_tgt->set_L(mps_src.get_L());
        if(mps_src.isCenter()) mps_tgt->set_LC(mps_src.get_LC());
        mps_tgt++;
    }
    state.clear_cache();
    state.clear_measurements();
}


void tools::finite::mps::normalize_state(class_state_finite &state, std::optional<long> chi_lim, std::optional<double> svd_threshold) {
    // When a state needs to be normalized it's enough to "move" the center position around the whole chain
    // Each move performs an SVD decomposition which leaves unitaries after it, effectively normalizing the state.

    tools::log->trace("Normalizing state");

    size_t num_moves = 2 * (state.get_length() - 2);
    if(state.has_nan()) throw std::runtime_error("State has NAN's before normalization");
    tools::log->info("Norm                 before normalization: {:.16f}", tools::finite::measure::norm(state));
    tools::log->info("Spin components      before normalization: {}", tools::finite::measure::spin_components(state));
    tools::log->info("Bond dimensions      before normalization: {}", tools::finite::measure::bond_dimensions(state));
    tools::log->info("Entanglement entropy before normalization: {}", tools::finite::measure::entanglement_entropies(state));

    for(size_t move = 0; move < num_moves; move++)
        move_center_point(state,chi_lim,svd_threshold);

    if(state.has_nan()) throw std::runtime_error("State has NAN's after normalization");
    tools::log->info("Norm                 after  normalization: {:.16f}", tools::finite::measure::norm(state));
    tools::log->info("Spin components      after  normalization: {}", tools::finite::measure::spin_components(state));
    tools::log->info("Bond dimensions      after  normalization: {}", tools::finite::measure::bond_dimensions(state));
    tools::log->info("Entanglement entropy after  normalization: {}", tools::finite::measure::entanglement_entropies(state));
    std::cerr << "MUST REBUILD ENVIRONMENTS NORMALIZATION" << std::endl;
    //    tools::finite::mps::rebuild_edges(state);
}



void tools::finite::mps::random_product_state(class_state_finite &state, const std::string &axis, const long state_number,
                                              const bool use_pauli_eigenstates)
/*!
 * There are many ways to random_product_state an initial product state state, based on the
 * arguments (axis,state_number,use_pauli_eigenstates) = (string,long,true/false).
 * Let "+-sector" mean one of {"x","+x","-x","y","+y","-y", "z","+z","-z"}.

        a) ("+-sector"  ,+- ,t,f)   Set spinors to a random sequence of eigenvectors (up/down) of either
                                    sx, sy or sz pauli matrices (same pauli for all sites). If the global
                                    sign (+-) is omitted, a random sign is chosen with equal probabilities.
                                    In the x and z cases the full state will turn out to be entirely real,
                                    which improves performance.

        b) ("random"    ,+- ,f,f)   Set each spinor randomly on C2


        c) ("+-sector"  ,+- ,f,f)   Set each spinor randomly on C2 (i.e. case b) and then project the  full state
                                    to the given parity sector. If the global sign (+-) is omitted,  a random
                                    sign is chosen with equal probabilities. As a consequence of this, the
                                    full state will have always have nonzero imaginary part.

        d) ("randomAxis",+- ,f,f)   Randomly select one of {"x","y","z"} and go to case a).
        e) ("none"      ,+- ,f,f)   Does not random_product_state
        f) ("+-sector"  ,>=0,?,t)   Interpret seed_state as bitfield "01100010110..." and interpret these as
                                    up(0)/down(1) of either sx, sy or sz pauli matrices (same pauli for all sites)
 * Note: seed_state is only used if >= 0.
 * Note: we "use" the seed_state only once. Subsequent calls do not keep resetting the seed.
*/
{
    tools::log->debug("Randomizing mps into sector {}", axis);
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_have_been_updated(false);

    if(state_number >= 0)
        internals::set_product_state_in_parity_sector_from_bitset(state, axis, state_number);
    else
        internals::set_product_state_randomly(state, axis, use_pauli_eigenstates);
    std::cerr << "MUST REBUILD ENVIRONMENTS AFTER RANDOM PRODUCT STATE INIT" << std::endl;
    //    tools::finite::mps::rebuild_edges(state);
}

void tools::finite::mps::random_current_state(class_state_finite &state, const std::string &axis1, const std::string &axis2) {
    Eigen::MatrixXcd paulimatrix1;
    Eigen::MatrixXcd paulimatrix2;
    if(axis1 == "x")
        paulimatrix1 = qm::spinOneHalf::sx;
    else if(axis1 == "y")
        paulimatrix1 = qm::spinOneHalf::sy;
    else if(axis1 == "z")
        paulimatrix1 = qm::spinOneHalf::sz;
    else
        paulimatrix1 = qm::spinOneHalf::Id;
    if(axis2 == "x")
        paulimatrix2 = qm::spinOneHalf::sx;
    else if(axis2 == "y")
        paulimatrix2 = qm::spinOneHalf::sy;
    else if(axis2 == "z")
        paulimatrix2 = qm::spinOneHalf::sz;
    else
        paulimatrix2 = qm::spinOneHalf::Id;
    //    auto [mpos,L,R] = qm::mpo::random_pauli_mpos(paulimatrix,state.get_length());
    auto chi_lim      = state.find_largest_chi();
    auto [mpos, L, R] = qm::mpo::random_pauli_mpos_x2(paulimatrix1, paulimatrix2, state.get_length());
    tools::finite::ops::apply_mpos(state, mpos, L, R);
    tools::finite::mps::normalize_state(state, chi_lim);
    tools::finite::debug::check_integrity(state);
    state = tools::finite::ops::get_projection_to_closest_parity_sector(state, "x");
}

void tools::finite::mps::project_to_closest_parity_sector(class_state_finite &state, const std::string &parity_sector) {
    state = tools::finite::ops::get_projection_to_closest_parity_sector(state, parity_sector);
}
