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

#include <tools/common/svd.h>
#include <utility>

void tools::finite::mps::initialize(class_state_finite &state, ModelType model_type, size_t num_sites, size_t position) {
    log->info("Initializing mps with {} sites at position {}", num_sites, position);
    if(num_sites < 2) throw std::logic_error("Tried to initialize MPS with less than 2 sites");
    if(num_sites > 2048) throw std::logic_error("Tried to initialize MPS with more than 2048 sites");
    if(position >= num_sites) throw std::logic_error("Tried to initialize MPS at a position larger than the number of sites");

    size_t spin_dim = 2; // Default is a two-level system
    switch(model_type) {
        case ModelType::ising_tf_rf: spin_dim = settings::model::ising_tf_rf::spin_dim; break;
        case ModelType::ising_sdual: spin_dim = settings::model::ising_sdual::spin_dim; break;
        default: spin_dim = 2;
    }

    state.MPS.clear();

    // Generate a simple MPS with all spins equal
    Eigen::Tensor<Scalar, 3> M(static_cast<long>(spin_dim), 1, 1);
    Eigen::Tensor<Scalar, 1> L(1);
    M(0, 0, 0) = 0;
    M(1, 0, 0) = 1;
    L(0)       = 1;
    for(size_t site = 0; site < num_sites; site++) {
        state.MPS.emplace_back(class_mps_site(M, L, site));
        if(site == position) state.MPS.back().set_LC(L);
    }
    if(state.MPS.size() != num_sites) throw std::logic_error("Initialized MPS with wrong size");
    if(not state.get_mps(position).isCenter()) throw std::logic_error("Initialized center matrix at the wrong position");
    if(state.get_position() != position) throw std::logic_error("Initialized MPS at the wrong position");
    state.site_update_tags = std::vector<bool>(num_sites, false);
}

void tools::finite::mps::random_product_state(class_state_finite &state, const std::string &parity_sector, const long state_number,
                                              const bool use_pauli_eigenstates)
/*!
 * There are many ways to random_product_state an initial product state state, based on the
 * arguments (parity_sector,state_number,use_pauli_eigenstates) = (string,long,true/false).
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
    tools::log->debug("Randomizing mps into sector {}", parity_sector);
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_have_been_updated(false);

    if(state_number >= 0)
        internals::set_product_state_in_parity_sector_from_bitset(state, parity_sector, state_number);
    else
        internals::set_product_state_randomly(state, parity_sector, use_pauli_eigenstates);
    std::cerr << "MUST REBUILD ENVIRONMENTS AFTER RANDOM PRODUCT STATE INIT" << std::endl;
    //    tools::finite::mps::rebuild_edges(state);
}

void tools::finite::mps::random_current_state(class_state_finite &state, const std::string &parity_sector1, const std::string &parity_sector2) {
    Eigen::MatrixXcd paulimatrix1;
    Eigen::MatrixXcd paulimatrix2;
    if(parity_sector1 == "x")
        paulimatrix1 = qm::spinOneHalf::sx;
    else if(parity_sector1 == "y")
        paulimatrix1 = qm::spinOneHalf::sy;
    else if(parity_sector1 == "z")
        paulimatrix1 = qm::spinOneHalf::sz;
    else
        paulimatrix1 = qm::spinOneHalf::Id;
    if(parity_sector2 == "x")
        paulimatrix2 = qm::spinOneHalf::sx;
    else if(parity_sector2 == "y")
        paulimatrix2 = qm::spinOneHalf::sy;
    else if(parity_sector2 == "z")
        paulimatrix2 = qm::spinOneHalf::sz;
    else
        paulimatrix2 = qm::spinOneHalf::Id;
    //    auto [mpos,L,R] = qm::mpo::random_pauli_mpos(paulimatrix,state.get_length());
    auto chi_lim      = state.find_largest_chi();
    auto [mpos, L, R] = qm::mpo::random_pauli_mpos_x2(paulimatrix1, paulimatrix2, state.get_length());
    tools::finite::ops::apply_mpos(state, mpos, L, R);
    tools::finite::mps::normalize(state, chi_lim);
    tools::finite::debug::check_integrity(state);
    state = tools::finite::ops::get_projection_to_closest_parity_sector(state, "x");
}

void tools::finite::mps::project_to_closest_parity_sector(class_state_finite &state, const std::string &parity_sector) {
    state = tools::finite::ops::get_projection_to_closest_parity_sector(state, parity_sector);
}

void tools::finite::mps::move_center_point(class_state_finite &state, std::optional<size_t> chi_lim, std::optional<double> svd_threshold) {
    if(state.position_is_any_edge()) {
        // Instead of moving out of the chain, just flip the direction and return
        state.flip_direction();
    } else {
        auto &MPS = state.MPS;
        if(MPS.empty()) throw std::runtime_error("MPS is empty");

        size_t      pos  = state.get_position();
        size_t      posL = state.get_direction() == 1 ? pos + 1 : pos - 1;
        size_t      posR = state.get_direction() == 1 ? pos + 2 : pos;
        auto &      mps  = state.get_mps();
        const auto &mpsL = state.get_mps(posL); // Becomes the new center position
        const auto &mpsR = state.get_mps(posR); // The site to the right of the new center position
        long        dL   = mpsL.spin_dim();
        long        dR   = mpsR.spin_dim();
        long        chiL = mpsL.get_chiL();
        long        chiR = mpsR.get_chiR();
        // Store the special LC bond in a temporary.
        Eigen::Tensor<Scalar, 1> LC = mps.get_LC();
        Eigen::Tensor<Scalar, 4> mps_2site;
        if(state.get_direction() == 1) {
            mps_2site = Textra::asDiagonal(LC)
                            .contract(mpsL.get_M_bare(), Textra::idx({1}, {1}))
                            .contract(mpsR.get_M_bare(), Textra::idx({2}, {1}))
                            .shuffle(Textra::array4{1, 2, 0, 3})
                            .reshape(Textra::array3{dL * dR, chiL, chiR});
        } else {
            mps_2site = mpsL.get_M_bare()
                            .contract(mpsR.get_M_bare(), Textra::idx({2}, {1}))
                            .contract(Textra::asDiagonal(LC), Textra::idx({3}, {0}))
                            .shuffle(Textra::array4{0, 2, 1, 3})
                            .reshape(Textra::array3{dL * dR, chiL, chiR});
        }
        tools::finite::mps::merge_multisite_mps(state, mps_2site, {posL, posR}, posL, chi_lim, svd_threshold);
        state.clear_cache();
        state.clear_measurements();
    }
}

// void tools::finite::svd::move_center_point(class_state_finite & state, std::optional<size_t> chi_lim, std::optional<double> svd_threshold){
//    if (state.position_is_any_edge()){
//        // Instead of moving out of the chain, just flip the direction and return
//        state.flip_direction();
//    }else{
//        auto & MPS_L  = state.MPS_L;
//        auto & MPS_R  = state.MPS_R;
//        auto & MPO_L  = state.MPO_L;
//        auto & MPO_R  = state.MPO_R;
//        auto & ENV_L  = state.ENV_L;
//        auto & ENV_R  = state.ENV_R;
//        auto & ENV2_L = state.ENV2_L;
//        auto & ENV2_R = state.ENV2_R;
//        if(ENV_L.empty()) throw std::runtime_error("ENVL is empty");
//        if(ENV_R.empty()) throw std::runtime_error("ENVR is empty");
//        if(MPS_L.empty()) throw std::runtime_error("MPSL is empty");
//        if(MPS_R.empty()) throw std::runtime_error("MPSR is empty");
//        if(MPS_L.back().get_position()  != ENV_L.back().get_position())  throw std::runtime_error("MPSL and ENVL have mismatching positions");
//        if(MPS_R.front().get_position() != ENV_R.front().get_position()) throw std::runtime_error("MPSR and ENVR have mismatching positions");
//        if(ENV_L.size() + ENV_R.size() != state.get_length()) throw std::runtime_error("ENVL + ENVR sizes do not add up to chain length");
//        if(MPS_L.size() + MPS_R.size() != state.get_length()) throw std::runtime_error("MPSL + MPSR sizes do not add up to chain length");
//        if(MPO_L.size() + MPO_R.size() != state.get_length()) throw std::runtime_error("MPOL + MPOR sizes do not add up to chain length");
//        assert(ENV_L.back().sites + ENV_R.front().sites == state.get_length() - 2);
//        //Store the special LC bond in a temporary.
//        Eigen::Tensor<Scalar,1> LC = MPS_L.back().get_LC();
//        MPS_L.back().unset_LC();
//
//        if (state.get_direction() == 1){
//            ENV_L .emplace_back(ENV_L .back().enlarge(MPS_L.back(), *MPO_L.back()));
//            ENV2_L.emplace_back(ENV2_L.back().enlarge(MPS_L.back(), *MPO_L.back()));
//            MPS_L.emplace_back(class_mps_site(MPS_R.front().get_M(), LC, MPS_R.front().get_position()));
//            MPO_L.emplace_back(MPO_R.front()->clone());
//            MPS_R.pop_front();
//            MPO_R.pop_front();
//            ENV_R.pop_front();
//            ENV2_R.pop_front();
//            Eigen::Tensor<Scalar,4> theta =
//                Textra::asDiagonal(LC)
//                    .contract(state.MPS_L.back().get_M(), Textra::idx({1},{1}))
//                    .contract(state.MPS_R.front().get_M(), Textra::idx({2},{1}))
//                    .shuffle(Textra::array4{1,0,2,3});
//            tools::finite::svd::truncate_theta(theta,state,chi_lim,svd_threshold);
//        }else{
//            ENV_R .emplace_front(ENV_R .front().enlarge(MPS_R.front(), *MPO_R.front()));
//            ENV2_R.emplace_front(ENV2_R.front().enlarge(MPS_R.front(), *MPO_R.front()));
//            MPS_R.emplace_front(class_mps_site(MPS_L.back().get_M(), LC, MPS_L.back().get_position()));
//            MPO_R.emplace_front(MPO_L.back()->clone());
//            MPS_L.pop_back();
//            MPO_L.pop_back();
//            ENV_L.pop_back();
//            ENV2_L.pop_back();
//            Eigen::Tensor<Scalar,4> theta =
//                state.MPS_L.back().get_M()
//                    .contract(state.MPS_R.front().get_M(), Textra::idx({2},{1}))
//                    .contract(Textra::asDiagonal(LC), Textra::idx({3},{0}));
//            tools::finite::svd::truncate_theta(theta,state,chi_lim,svd_threshold);
//        }
////        size_t pos = state.get_position();
////        tools::log->trace("SVD site {:2} log₁₀ trunc: {:12.8f} χlim: {:4} χ: {:4}", pos, std::log10(state.get_truncation_error(pos)),state.get_chi_lim(),
/// tools::finite::measure::bond_dimension_current(state));
//
//        assert(MPO_L.size() + MPO_R.size() == state.get_length());
//        if(ENV_L.empty()) throw std::runtime_error("ENVL became empty");
//        if(ENV_R.empty()) throw std::runtime_error("ENVR became empty");
//        if(MPS_L.empty()) throw std::runtime_error("MPSL became empty");
//        if(MPS_R.empty()) throw std::runtime_error("MPSR became empty");
//        if(MPS_L.back().get_position()  != ENV_L.back().get_position())  throw std::runtime_error("MPSL and ENVL got mismatching positions");
//        if(MPS_R.front().get_position() != ENV_R.front().get_position()) throw std::runtime_error("MPSR and ENVR got mismatching positions");
//        if(ENV_L.size() + ENV_R.size() != state.get_length()) throw std::runtime_error("ENVL + ENVR sizes do not add up to chain length anymore");
//        if(MPS_L.size() + MPS_R.size() != state.get_length()) throw std::runtime_error("MPSL + MPSR sizes do not add up to chain length anymore");
//        if(MPO_L.size() + MPO_R.size() != state.get_length()) throw std::runtime_error("MPOL + MPOR sizes do not add up to chain length anymore");
//
//        state.clear_cache();
//        state.clear_measurements();
//    }
//}

void tools::finite::mps::normalize_state(class_state_finite &state, std::optional<size_t> chi_lim, std::optional<double> svd_threshold) {
    tools::log->trace("Normalizing state");
    using namespace Textra;
    using Scalar = class_state_finite::Scalar;
    state.clear_cache();
    state.clear_measurements();
    tools::common::profile::t_svd->tic();
    size_t num_moves = 2 * (state.get_length() - 2);
    if(state.has_nan()) throw std::runtime_error("State has NAN's before normalization");
    tools::log->info("Norm                 before normalization: {:.16f}", tools::finite::measure::norm(state));
    tools::log->info("Spin components      before normalization: {}", tools::finite::measure::spin_components(state));
    tools::log->info("Bond dimensions      before normalization: {}", tools::finite::measure::bond_dimensions(state));
    tools::log->info("Entanglement entropy before normalization: {}", tools::finite::measure::entanglement_entropies(state));

    // Start by truncating at the current position.
    tools::finite::svd::truncate_theta(state.get_theta(), state, chi_lim, svd_threshold);
    // Now start moving
    for(size_t moves = 0; moves < num_moves; moves++) {
        //    Check edge
        if(state.position_is_any_edge()) state.flip_direction();
        auto &MPS_L = state.MPS_L;
        auto &MPS_R = state.MPS_R;
        if(MPS_L.empty()) throw std::runtime_error("MPS_L is empty");
        if(MPS_R.empty()) throw std::runtime_error("MPS_R is empty");

        // Store the special LC bond in a temporary.
        Eigen::Tensor<Scalar, 1> LC = MPS_L.back().get_LC();
        MPS_L.back().unset_LC();
        if(state.get_direction() == 1) {
            MPS_L.emplace_back(class_mps_site(MPS_R.front().get_M(), LC, MPS_R.front().get_position()));
            MPS_R.pop_front();
            Eigen::Tensor<Scalar, 4> theta = Textra::asDiagonal(LC)
                                                 .contract(state.MPS_L.back().get_M(), Textra::idx({1}, {1}))
                                                 .contract(state.MPS_R.front().get_M(), Textra::idx({2}, {1}))
                                                 .shuffle(Textra::array4{1, 0, 2, 3});
            tools::finite::svd::truncate_theta(theta, state, chi_lim, svd_threshold);
        } else {
            MPS_R.emplace_front(class_mps_site(MPS_L.back().get_M(), LC, MPS_L.back().get_position()));
            MPS_L.pop_back();
            Eigen::Tensor<Scalar, 4> theta =
                state.MPS_L.back().get_M().contract(state.MPS_R.front().get_M(), Textra::idx({2}, {1})).contract(Textra::asDiagonal(LC), Textra::idx({3}, {0}));
            tools::finite::svd::truncate_theta(theta, state, chi_lim, svd_threshold);
        }

        if(MPS_L.empty()) throw std::runtime_error("MPS_L became empty");
        if(MPS_R.empty()) throw std::runtime_error("MPS_R became empty");
        if(MPS_L.size() + MPS_R.size() != state.get_length()) throw std::runtime_error("MPS_L + MPS_R sizes do not add up to chain length anymore");

        //        state.clear_measurements();
        //        tools::log->trace("Position {}", state.get_position());
        //        tools::log->trace("Bond dimension  after normalization: {}", tools::finite::measure::bond_dimension_current(state));
        //        tools::log->trace("Bond dimensions after normalization: {}", tools::finite::measure::bond_dimensions(state));
    }
    tools::common::profile::t_svd->toc();
    state.clear_measurements();
    state.clear_cache(); // IMPORTANT, because we generate a theta earlier in this function and it needs to be invalidated
    if(state.has_nan()) throw std::runtime_error("State has NAN's after normalization");
    tools::log->info("Norm                 after  normalization: {:.16f}", tools::finite::measure::norm(state));
    tools::log->info("Spin components      after  normalization: {}", tools::finite::measure::spin_components(state));
    tools::log->info("Bond dimensions      after  normalization: {}", tools::finite::measure::bond_dimensions(state));
    tools::log->info("Entanglement entropy after  normalization: {}", tools::finite::measure::entanglement_entropies(state));
    std::cerr << "MUST REBUILD ENVIRONMENTS NORMALIZATION" << std::endl;
    //    tools::finite::mps::rebuild_edges(state);
}

void tools::finite::mps::merge_multisite_mps(class_state_finite &state, const Eigen::Tensor<Scalar, 3> &multisite_mps, const std::list<size_t> &positions,
                                             size_t center_position, std::optional<size_t> chi_lim, std::optional<double> svd_threshold) {
    // Some sanity checks
    if(multisite_mps.dimension(1) != state.get_mps(positions.front()).get_chiL())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim1 {} != chiL on left-most site {}", multisite_mps.dimension(1),
                                             state.get_mps(positions.front()).get_chiL(), positions.front()));

    if(multisite_mps.dimension(2) != state.get_mps(positions.back()).get_chiR())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim2 {} != chiR on right-most site {}", multisite_mps.dimension(2),
                                             state.get_mps(positions.back()).get_chiR(), positions.back()));

    std::list<long> spin_dims;
    for(const auto &site : positions) spin_dims.emplace_back(state.get_mps(site).spin_dim());

    // Split the multisite mps into single-site mps objects
    auto mps_list = tools::common::svd::split_mps(multisite_mps, spin_dims, positions, center_position, state.get_chi_lim(), svd_threshold);

    if(positions.size() != mps_list.size())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: number of sites mismatch: positions.size() {} != mps_list.size() {}",
                                             positions.size(), mps_list.size()));

    // Note that one of the positions on the split will be a center, so we need to unset
    // the center in our current state so we don't get duplicate centers
    state.get_mps().unset_LC();

    // Copy the split up mps components into the current state
    auto mps_tgt = std::next(state.MPS.begin(), static_cast<long>(positions.front()));
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
