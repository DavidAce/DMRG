#include <config/debug.h>
#include <general/iter.h>
#include <io/fmt.h>
#include <math/rnd.h>
#include <math/tenx.h>
#include <qm/mpo.h>
#include <qm/spin.h>
#include <tensors/site/mpo/MpoSite.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>

//
#include <debug/exceptions.h>
#include <math/linalg/tensor.h>
#include <tools/common/contraction.h>

namespace settings {
    inline constexpr bool debug_projection = false;
}

void tools::finite::ops::apply_mpo(StateFinite &state, const Eigen::Tensor<Scalar, 4> &mpo, const Eigen::Tensor<Scalar, 3> &Ledge,
                                   const Eigen::Tensor<Scalar, 3> &Redge) {
    std::vector<Eigen::Tensor<Scalar, 4>> mpos(state.get_length(), mpo);
    apply_mpos(state, mpos, Ledge, Redge);
}

void tools::finite::ops::apply_mpos(StateFinite &state, const std::vector<Eigen::Tensor<Scalar, 4>> &mpos, const Eigen::Tensor<Scalar, 1> &Ledge,
                                    const Eigen::Tensor<Scalar, 1> &Redge) {
    Eigen::Tensor<Scalar, 3> Ledge3, Redge3;
    {
        auto mps_dims = state.mps_sites.front()->dimensions();
        auto mpo_dims = mpos.front().dimensions();

        long mpsDim = mps_dims[1];
        long mpoDim = mpo_dims[0];
        Ledge3.resize(tenx::array3{mpsDim, mpsDim, mpoDim});
        Ledge3.setZero();
        for(long i = 0; i < mpsDim; i++) {
            std::array<long, 1> extent1                     = {mpoDim};
            std::array<long, 3> offset3                     = {i, i, 0};
            std::array<long, 3> extent3                     = {1, 1, mpoDim};
            Ledge3.slice(offset3, extent3).reshape(extent1) = Ledge;
        }
    }
    {
        auto mps_dims = state.mps_sites.back()->dimensions();
        auto mpo_dims = mpos.back().dimensions();

        long mpsDim = mps_dims[2];
        long mpoDim = mpo_dims[1];
        Redge3.resize(tenx::array3{mpsDim, mpsDim, mpoDim});
        Redge3.setZero();
        for(long i = 0; i < mpsDim; i++) {
            std::array<long, 1> extent1                     = {mpoDim};
            std::array<long, 3> offset3                     = {i, i, 0};
            std::array<long, 3> extent3                     = {1, 1, mpoDim};
            Redge3.slice(offset3, extent3).reshape(extent1) = Redge;
        }
    }
    apply_mpos(state, mpos, Ledge3, Redge3);
}

void tools::finite::ops::apply_mpos(StateFinite &state, const std::vector<Eigen::Tensor<Scalar, 4>> &mpos, const Eigen::Tensor<Scalar, 3> &Ledge,
                                    const Eigen::Tensor<Scalar, 3> &Redge) {
    // Apply MPO's on Gamma matrices and
    // increase the size on all Lambdas by chi*mpoDim
    tools::log->trace("Applying MPO's");
    if(mpos.size() != state.get_length()) throw std::runtime_error("Number of mpo's doesn't match the number of sites on the system");

    if constexpr(settings::debug or settings::debug_projection) {
        state.clear_measurements();
        tools::log->debug("Num mpos             before applying mpos: {}", mpos.size());
        tools::log->debug("Norm                 before applying mpos: {:.16f}", tools::finite::measure::norm(state));
        tools::log->debug("Spin components      before applying mpos: {}", tools::finite::measure::spin_components(state));
        tools::log->debug("Bond dimensions      before applying mpos: {}", tools::finite::measure::bond_dimensions(state));
        tools::log->debug("Entanglement entropy before applying mpos: {}", tools::finite::measure::entanglement_entropies(state));
        if(tenx::hasNaN(tools::finite::measure::entanglement_entropies(state)))
            throw except::runtime_error("Entanglement entropy has nans:\n{}", tools::finite::measure::entanglement_entropies(state));
    }
    state.clear_measurements();
    for(const auto &[pos, mpo] : iter::enumerate(mpos)) state.get_mps_site<size_t>(pos).apply_mpo(mpo); // Apply all mpo's

    // Take care of the edges. Apply the left and right MPO-edges on A's and B's
    // so the left- and right-most lambdas become ones.
    // Looking at the pictures, the A/B mps tensors have already been contracted with mpo's
    // in the for loop above, where their horizontal legs have been merged as well.
    // Hence, asking for get_M_bare() will get us the A/B with mpo's connected.

    // Note:
    // If M is an A then M_bare = L0*M,, and the mpo was applied on A.
    // If M is an AC then M_bare = L0*M, (L1 = LC) and the mpo was applied on M_bare, LC dim increased by mpodim
    // If M is a B then M_bare = M * L1, and there is no L0.
    {
        /*
         *    |------ 0   0---[L0]---1   1---[M]---2   0---[L1]---1--|
         *    |                               |              |       |
         *    |                               0              |       |----
         *    |                               2              |       |
         *    |                               |              |       |
         *  [Ledge]--- 2                0---[mpo]---1  0---[ I ]--1--|
         *    |                               |
         *    |                               2
         *    |
         *    |
         *    |------ 1
         */
        auto                    &mps      = state.get_mps_site(0ul);
        auto                     isCenter = mps.isCenter();
        auto                     label    = mps.get_label();
        long                     mpoDimL  = mpos.front().dimension(0);
        auto                     Ldim     = Ledge.dimension(0);
        Eigen::Tensor<Scalar, 3> M_temp =
            Ledge
                .shuffle(tenx::array3{0, 2, 1})                  // Start by shuffling the legs into consecutive order before merge
                .reshape(tenx::array2{Ldim * mpoDimL, Ldim})     // Merge the legs
                .contract(mps.get_M_bare(), tenx::idx({0}, {1})) // Contract with M which already has the mpo on it (not including LC, possibly)
                .shuffle(tenx::array3{1, 0, 2});                 // Shuffle back to convention
        if(isCenter or label != "B") {
            Eigen::Tensor<Scalar, 1> one = Eigen::Tensor<Scalar, 1>(Ldim).constant(1.0);
            mps.set_mps(M_temp, one, 0, label);
        } else {
            // The left edge is a B-site.
            // Then L1 dim has already been increased and L0 doesn't exist. (If it did, it would just be a "1")
            mps.set_M(M_temp);
        }
    }
    {
        // Note:
        // If M is an A then M_bare = (L-1)*M, and the mpo was applied on A (this should never happen)
        // If M is an AC then M_bare = (L-1)*M, (L = LC) and the mpo was applied on M_bare,L dim incresed by mpoDim and LC
        // If M is a B then M_bare = M * L (L-1 belongs on the site to the right), and the mpo was applied on M_bare

        /*
         *  |--0---[L-1]---1  1---[M]---2  0---[L]---1   0 ------|
         *  |                      |                             |
         *  |                      0                             |
         *--|                      2                             |
         *  |                      |                             |
         *  |--0---[ I ]---1  0--[mpo]--1                2 ---[Redge]
         *                          |                             |
         *                          2                             |
         *                                                        |
         *                                                1 ------|
         */
        auto                    &mps      = state.get_mps_site(state.get_length() - 1);
        auto                     label    = mps.get_label();
        bool                     isCenter = mps.isCenter();
        long                     mpoDimR  = mpos.back().dimension(1);
        auto                     Rdim     = Redge.dimension(0);
        Eigen::Tensor<Scalar, 3> M_temp   = Redge.shuffle(tenx::array3{0, 2, 1})
                                              .reshape(tenx::array2{Rdim * mpoDimR, Rdim})
                                              .contract(mps.get_M(), tenx::idx({0}, {2})) // Include LC if it's there
                                              .shuffle(tenx::array3{1, 2, 0});
        Eigen::Tensor<Scalar, 1> one = Eigen::Tensor<Scalar, 1>(Rdim).constant(1.0);
        if(isCenter) {
            mps.set_M(M_temp);
            mps.set_LC(one, 0);
        } else if(label == "A") {
            mps.set_M(M_temp);
        } else if(label == "B") {
            mps.set_M(M_temp);
            mps.set_L(one, 0);
        }
    }

    if constexpr(settings::debug or settings::debug_projection) {
        state.clear_measurements();
        state.clear_cache();
        state.tag_all_sites_normalized(false); // This operation denormalizes all sites
        tools::log->debug("Num mpos             after  applying mpos: {}", mpos.size());
        tools::log->debug("Norm                 after  applying mpos: {:.16f}", tools::finite::measure::norm(state));
        tools::log->debug("Spin components      after  applying mpos: {}", tools::finite::measure::spin_components(state));
        tools::log->debug("Bond dimensions      after  applying mpos: {}", tools::finite::measure::bond_dimensions(state));
        tools::log->debug("Entanglement entropy after  applying mpos: {}", tools::finite::measure::entanglement_entropies(state));
        if(tenx::hasNaN(tools::finite::measure::entanglement_entropies(state)))
            throw except::runtime_error("Entanglement entropy has nans:\n{}", tools::finite::measure::entanglement_entropies(state));
    }

    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
}

void tools::finite::ops::project_to_sector(StateFinite &state, const Eigen::MatrixXcd &paulimatrix, int sign, std::optional<long> bond_limit,
                                           std::optional<svd::settings> svd_settings) {
    // This function applies the projection MPO operator  "0.5 * ( 1 - prod s)", where
    // 1 is understood as a 2^L x 2^L tensor and "prod s" is the outer product of pauli matrices, one for each site.
    // This operation leaves the global norm unchanged (thanks to the 0.5 factor) but locally each MPS loses its
    // norm (norm = 2), and entanglement entropies become doubled.
    // Therefore a proper full normalization is required after this operation, as well as a full
    // rebuild of environments.

    if(std::abs(sign) != 1) throw std::runtime_error(fmt::format("Expected 'sign' +1 or -1. Got [{}]", sign));
    tools::log->debug("Projecting state to sector with sign {}", sign);
    auto t_prj = tid::tic_scope("projection");
    tools::finite::mps::normalize_state(state, bond_limit, svd_settings, NormPolicy::IFNEEDED);

    auto spin_components = tools::finite::measure::spin_components(state);
    tools::log->debug("Spin components before projection : X = {:.16f}  Y = {:.16f}  Z = {:.16f}", spin_components[0], spin_components[1], spin_components[2]);
    state.clear_measurements();
    state.clear_cache();
    // Do the projection
    const auto [mpos, L, R] = qm::mpo::parity_projector_mpos(paulimatrix, state.get_length(), sign);
    apply_mpos(state, mpos, L, R);
    tools::finite::mps::normalize_state(state, bond_limit, svd_settings, NormPolicy::ALWAYS); // Has to be normalized ALWAYS, projection ruins normalization!
    spin_components = tools::finite::measure::spin_components(state);
    tools::log->debug("Spin components after  projection : X = {:.16f}  Y = {:.16f}  Z = {:.16f}", spin_components[0], spin_components[1], spin_components[2]);
    if constexpr(settings::debug) state.assert_validity();
}

std::optional<double> tools::finite::ops::get_spin_component_in_sector(StateFinite &state, std::string_view sector) {
    auto                                      t_align         = tid::tic_scope("align");
    if(mps::init::is_valid_axis(sector)) {
        return tools::finite::measure::spin_component(state, mps::init::get_pauli(sector));
    } else
        return std::nullopt;
}

void tools::finite::ops::project_to_nearest_sector(StateFinite &state, std::string_view sector, std::optional<long> bond_limit,
                                                   std::optional<svd::settings> svd_settings) {
    /*
     * When projecting, there is one bad thing that may happen: that the norm of the state vanishes.
     *
     * This can happen in a couple of different scenarios:
     *      - The global state has spin component X = -1.0 and we project to X = +1.0 (or vice versa)
     *
     * Therefore the projection is only done if the state has a chance of surviving
     * Otherwise emit a warning and return
     *
     */

    tools::log->debug("Projecting state to axis nearest sector {}", sector);
    auto t_prj                    = tid::tic_scope("proj");
    auto spin_component_in_sector = get_spin_component_in_sector(state, sector);
    if(spin_component_in_sector.has_value()) {
        auto sector_sign    = mps::init::get_sign(sector);
        auto paulimatrix    = mps::init::get_pauli(sector);
        auto spin_alignment = sector_sign * spin_component_in_sector.value();
        // Now we have to check that the intended projection is safe
        tools::log->debug("Spin component in sector {}: {:.16f}", sector, spin_component_in_sector.value());
        if(spin_alignment > 0)
            // In this case the state has an aligned component along the requested axis --> safe
            project_to_sector(state, paulimatrix, sector_sign, bond_limit, svd_settings);
        else if(spin_alignment < 0) {
            constexpr auto spin_alignment_threshold = 1e-3;
            // In this case the state has an anti-aligned component along the requested axis --> safe if spin_component < 1 - spin_component_threshold
            // Remember that  spin_alignment == -1 means orthogonal!
            if(std::abs(spin_alignment) < 1.0 - spin_alignment_threshold)
                project_to_sector(state, paulimatrix, sector_sign, bond_limit, svd_settings);
            else
                return tools::log->warn("Skipping projection to [{0}]: State spin is orthogonal to the requested projection axis: <{0}|Ψ> = {1:.16f}", sector,
                                        spin_alignment);
        } else if(spin_alignment == 0) {
            // No sector sign was specified, so we select the one along which there is a component
            if(spin_component_in_sector.value() >= 0)
                sector_sign = 1;
            else
                sector_sign = -1;
            project_to_sector(state, paulimatrix, sector_sign, bond_limit, svd_settings);
        }
    } else if(sector == "randomAxis") {
        std::vector<std::string> possibilities = {"x", "y", "z"};
        std::string              chosen_axis   = possibilities[rnd::uniform_integer_box<size_t>(0, possibilities.size() - 1)];
        project_to_nearest_sector(state, chosen_axis, bond_limit, svd_settings);
    } else if(sector == "random") {
        auto             coeffs    = Eigen::Vector3d::Random().normalized();
        Eigen::Matrix2cd random_c2 = coeffs(0) * qm::spin::half::sx + coeffs(1) * qm::spin::half::sy + coeffs(2) * qm::spin::half::sz;
        return project_to_sector(state, random_c2, 1, bond_limit, svd_settings);
    } else if(sector == "none") {
        return;
    } else
        throw std::runtime_error(fmt::format("Could not parse sector string [{}]", sector));
}

StateFinite tools::finite::ops::get_projection_to_sector(const StateFinite &state, const Eigen::MatrixXcd &paulimatrix, int sign,
                                                         std::optional<long> bond_limit, std::optional<svd::settings> svd_settings) {
    auto state_projected = state;
    project_to_sector(state_projected, paulimatrix, sign, bond_limit, svd_settings);
    return state_projected;
}

StateFinite tools::finite::ops::get_projection_to_nearest_sector(const StateFinite &state, std::string_view sector, std::optional<long> bond_limit,
                                                                 std::optional<svd::settings> svd_settings) {
    auto state_projected = state;
    project_to_nearest_sector(state_projected, sector, bond_limit, svd_settings);
    return state_projected;
}

double tools::finite::ops::overlap(const StateFinite &state1, const StateFinite &state2) {
    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");
    size_t pos     = 0;
    auto   overlap = tools::common::contraction::contract_mps_mps_partial(state1.get_mps_site(pos).get_M(), state2.get_mps_site(pos).get_M(), {0, 1});
    for(pos = 1; pos < state1.get_length(); pos++) {
        Eigen::Tensor<Scalar, 2> temp = overlap.contract(state1.get_mps_site(pos).get_M(), tenx::idx({0}, {1}))
                                            .contract(state2.get_mps_site(pos).get_M().conjugate(), tenx::idx({0, 1}, {1, 0}));
        overlap = temp;
    }

    double norm_chain = std::real(tenx::MatrixMap(overlap).trace());
    return norm_chain;
}
