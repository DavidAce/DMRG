//
// Created by david on 2019-01-29.
//

#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include <config/enums.h>
#include <config/nmspc_settings.h>
#include <general/nmspc_quantum_mechanics.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/finite/debug.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <tools/finite/svd.h>

#include "measure.h"
#include <math/nmspc_math.h>
#include <tools/common/prof.h>
#include <tools/common/svd.h>
#include <utility>

void tools::finite::mps::move_center_point(class_state_finite &state, long chi_lim, std::optional<double> svd_threshold) {
    if(state.position_is_any_edge()) {
        // Instead of moving out of the chain, just flip the direction and return
        state.flip_direction();
    } else {
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
        Eigen::Tensor<Scalar, 1> LC = mps.get_LC();
        Eigen::Tensor<Scalar, 3> twosite_tensor;
        if(state.get_direction() == 1) {
            // Here both M_bare are B's
            // i.e. mpsL.get_M() = GB * LB
            // and  mpsR.get_M() = GB * GB
            // So we have to attach LC from the left
            twosite_tensor = Textra::asDiagonal(LC)
                                 .contract(mpsL.get_M(), Textra::idx({1}, {1}))
                                 .contract(mpsR.get_M(), Textra::idx({2}, {1}))
                                 .shuffle(Textra::array4{1, 2, 0, 3})
                                 .reshape(Textra::array3{dL * dR, chiL, chiR});
        } else {
            // Here both M_bare are A's
            // The right A should be the previous position, so it has an attached
            // LC if we ask for get_M(), i.e. mpsR.get_M() = LA * GA * LC
            // The left A should be a simple A, i.e. mpsL.get_M() = LA * GA

            twosite_tensor =
                mpsL.get_M().contract(mpsR.get_M(), Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{dL * dR, chiL, chiR});
        }

        tools::finite::mps::merge_multisite_tensor(state, twosite_tensor, {posL, posR}, posL, chi_lim, svd_threshold);
        state.clear_cache();
        state.clear_measurements();

        // Put LC where it belongs.
        // Recall that mpsL, mpsR are on the new position, not the old one!
        if(state.get_direction() == 1) mpsL.set_L(LC);
        else
            mpsR.set_L(LC);
    }
}

void tools::finite::mps::merge_multisite_tensor(class_state_finite &state, const Eigen::Tensor<Scalar, 3> &multisite_mps, const std::list<size_t> &sites,
                                                size_t center_position, long chi_lim, std::optional<double> svd_threshold) {
    // Some sanity checks
    if(multisite_mps.dimension(1) != state.get_mps_site(sites.front()).get_chiL())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim1 {} != chiL on left-most site {}", multisite_mps.dimension(1),
                                             state.get_mps_site(sites.front()).get_chiL(), sites.front()));

    if(multisite_mps.dimension(2) != state.get_mps_site(sites.back()).get_chiR())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim2 {} != chiR on right-most site {}", multisite_mps.dimension(2),
                                             state.get_mps_site(sites.back()).get_chiR(), sites.back()));
    long            spin_prod = 1;
    std::list<long> spin_dims;
    for(const auto &site : sites) {
        spin_dims.emplace_back(state.get_mps_site(site).spin_dim());
        spin_prod *= spin_dims.back();
    }
    if(spin_prod != multisite_mps.dimension(0))
        throw std::runtime_error(
            fmt::format("Could not merge multisite mps into state: multisite_mps dim0 {} != spin_prod {}", multisite_mps.dimension(0), spin_prod));

    // Split the multisite mps into single-site mps objects
    auto mps_list = tools::common::svd::split_mps(multisite_mps, spin_dims, sites, center_position, chi_lim, svd_threshold);

    if(sites.size() != mps_list.size())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: number of sites mismatch: positions.size() {} != mps_list.size() {}",
                                             sites.size(), mps_list.size()));

    // Note that one of the positions on the split will be a center, so we need to unset
    // the center in our current state so we don't get duplicate centers
    state.get_mps_site().unset_LC();

    // Copy the split up mps components into the current state
    auto mps_ptr = std::next(state.mps_sites.begin(), static_cast<long>(sites.front()));
    for(const auto &mps_src : mps_list) {
        auto &mps_tgt = **mps_ptr;
        mps_tgt.merge_mps(mps_src);
        mps_ptr++;
    }
    state.clear_cache();
    state.clear_measurements();
}

bool tools::finite::mps::normalize_state(class_state_finite &state, long chi_lim, std::optional<double> svd_threshold) {
    // When a state needs to be normalized it's enough to "move" the center position around the whole chain.
    // Each move performs an SVD decomposition which leaves unitaries after it, effectively normalizing the state.
    // NOTE! It may be important to start with the current position.

    // We may want to make a quick check on release builds, but more thorough on debug, for performance.
    const double norm = [state]{
      if(settings::debug)
        return tools::finite::measure::norm(state);
      else
        return tools::finite::measure::norm_fast(state);
    }();

    // We may only go ahead with a normalization if its really needed.
    if(std::abs(norm - 1.0) < settings::precision::max_norm_error) return false;

    if(settings::debug) {
        tools::log->trace("Normalizing state");
        tools::log->info("Position             before normalization: {}", state.get_position());
        tools::log->info("Direction            before normalization: {}", state.get_direction());
        tools::log->info("Norm                 before normalization: {:.16f}", tools::finite::measure::norm(state));
        tools::log->info("Spin components      before normalization: {}", tools::finite::measure::spin_components(state));
        tools::log->info("Bond dimensions      before normalization: {}", tools::finite::measure::bond_dimensions(state));
        tools::log->info("Entanglement entropy before normalization: {}", tools::finite::measure::entanglement_entropies(state));
    }
    // Start with normalizing at the current position
    size_t                   num_moves = 2 * (state.get_length() - 1);
    size_t                   posL      = state.get_position();
    size_t                   posR      = posL + 1;
    auto &                   mpsL      = state.get_mps_site(posL);
    auto &                   mpsR      = state.get_mps_site(posR);
    long                     dL        = mpsL.spin_dim();
    long                     dR        = mpsR.spin_dim();
    long                     chiL      = mpsL.get_chiL();
    long                     chiR      = mpsR.get_chiR();
    Eigen::Tensor<Scalar, 3> twosite_tensor =
        mpsL.get_M().contract(mpsR.get_M(), Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{dL * dR, chiL, chiR});
    tools::finite::mps::merge_multisite_tensor(state, twosite_tensor, {posL, posR}, posL, chi_lim, svd_threshold);

    // Now we can move around the chain
    for(size_t move = 0; move < num_moves; move++) move_center_point(state, chi_lim, svd_threshold);
    state.clear_measurements();
    state.clear_cache();
    state.assert_validity();
    if(settings::debug) {
        tools::log->info("Position             after  normalization: {}", state.get_position());
        tools::log->info("Direction            after  normalization: {}", state.get_direction());
        tools::log->info("Norm                 after  normalization: {:.16f}", tools::finite::measure::norm(state));
        tools::log->info("Spin components      after  normalization: {}", tools::finite::measure::spin_components(state));
        tools::log->info("Bond dimensions      after  normalization: {}", tools::finite::measure::bond_dimensions(state));
        tools::log->info("Entanglement entropy after  normalization: {}", tools::finite::measure::entanglement_entropies(state));
    }
    return true;
}

void tools::finite::mps::random_product_state(class_state_finite &state, const std::string &sector, long bitfield, bool use_eigenspinors)
/*!
 * There are many ways to random_product_state an initial product state state, based on the
 * arguments (sector,state_number,use_eigenspinors) = (string,long,true/false).
 * Let sector="+-axis" mean one of {"x","+x","-x","y","+y","-y", "z","+z","-z"}.

        a) ("+-axis"  ,+- ,t,f)     Set spinors to a random sequence of eigenvectors (up/down) of either
                                    sx, sy or sz pauli matrices (same pauli for all sites). If the global
                                    sign (+-) is omitted, a random sign is chosen with equal probabilities.
                                    In the x and z cases the full state will turn out to be entirely real,
                                    which improves performance.

        b) ("random"    ,+- ,f,f)   Set each spinor randomly on C2


        c) ("+-axis"  ,+- ,f,f)   Set each spinor randomly on C2 (i.e. case b) and then project the  full state
                                    to the given parity sector. If the global sign (+-) is omitted,  a random
                                    sign is chosen with equal probabilities. As a consequence of this, the
                                    full state will have always have nonzero imaginary part.

        d) ("randomAxis",+- ,f,f)   Randomly select one of {"x","y","z"} and go to case a).
        e) ("none"      ,+- ,f,f)   Does not random_product_state
        f) ("+-axis"  ,>=0,?,t)   Interpret seed_state as bitfield "01100010110..." and interpret these as
                                    up(0)/down(1) of either sx, sy or sz pauli matrices (same pauli for all sites)
 * Note: seed_state is only used if >= 0.
 * Note: we "use" the seed_state only once. Subsequent calls do not keep resetting the seed.
*/
{
    tools::log->debug("Setting random product state in sector {}", sector);
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_have_been_updated(false);

    if(bitfield >= 0 and bitfield != internals::used_bitfield) {
        internals::set_random_product_state_in_sector_using_bitfield(state, sector, bitfield);
        internals::used_bitfield = bitfield;
    } else
        internals::set_random_product_state(state, sector, use_eigenspinors);
}

void tools::finite::mps::apply_random_paulis(class_state_finite &state, const std::vector<std::string> &paulistrings) {
    std::vector<Eigen::Matrix2cd> paulimatrices;
    for(const auto &str : paulistrings) paulimatrices.emplace_back(internals::get_pauli(str));
    auto [mpos, L, R] = qm::mpo::random_pauli_mpos(paulimatrices, state.get_length());
    tools::finite::ops::apply_mpos(state, mpos, L, R);
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

    tools::log->trace("Truncated all sites");
    std::cerr << "MUST REBUILD EDGES AFTER TRUNCATING ALL SITES" << std::endl;
}

void tools::finite::mps::truncate_active_sites(class_state_finite &state, long chi_lim, std::optional<double> svd_threshold) {
    tools::log->warn("Truncate active sites needs an implementation");
    throw std::runtime_error("Truncate active sites needs an implementation");
}

void tools::finite::mps::truncate_next_sites(class_state_finite &state, long chi_lim, size_t num_sites, std::optional<double> svd_threshold) {
    tools::log->warn("Truncate next sites needs an implementation");
    throw std::runtime_error("Truncate next sites needs an implementation");
}
