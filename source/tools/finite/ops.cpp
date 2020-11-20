//
// Created by david on 2019-01-30.
//

#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <math/rnd.h>
#include <physics/nmspc_quantum_mechanics.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/debug.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>

using Scalar = std::complex<double>;
using namespace Textra;

void tools::finite::ops::apply_mpo(class_state_finite &state, const Eigen::Tensor<Scalar, 4> &mpo, const Eigen::Tensor<Scalar, 3> &Ledge,
                                   const Eigen::Tensor<Scalar, 3> &Redge) {
    std::list<Eigen::Tensor<Scalar, 4>> mpos(state.get_length(), mpo);
    apply_mpos(state, mpos, Ledge, Redge);
}

void tools::finite::ops::apply_mpos(class_state_finite &state, const std::list<Eigen::Tensor<Scalar, 4>> &mpos, const Eigen::Tensor<Scalar, 3> &Ledge,
                                    const Eigen::Tensor<Scalar, 3> &Redge) {
    // Apply MPO's on Gamma matrices and
    // increase the size on all Lambdas by chi*mpoDim
    tools::log->trace("Applying MPO's");
    if(mpos.size() != state.get_length()) throw std::runtime_error("Number of mpo's doesn't match the number of sites on the system");

    if constexpr(settings::debug){
        state.clear_measurements();
        tools::log->debug("Norm                 before applying mpos: {:.16f}", tools::finite::measure::norm(state));
        tools::log->debug("Spin components      before applying mpos: {}", tools::finite::measure::spin_components(state));
        tools::log->debug("Bond dimensions      before applying mpos: {}", tools::finite::measure::bond_dimensions(state));
        tools::log->debug("Entanglement entropy before applying mpos: {}", tools::finite::measure::entanglement_entropies(state));
    }
    state.clear_measurements();


    auto mpo = mpos.begin();
    for(size_t pos = 0; pos < state.get_length(); pos++) {
        state.get_mps_site(pos).apply_mpo(*mpo);
        mpo++;
    }

    // Take care of the edges. Apply the left and right MPO-edges on A's and B's
    // so the left- and right-most lambdas become ones.
    // Looking at the pictures, the A/B mps tensors have already been contracted with mpo's
    // in the for loop above, where their horizontal legs have been merged as well.
    // Hence, asking for get_M_bare() will get us the A/B with mpo's connected.
    {
        /*
         *    |------ 0     1---[A]---2
         *    |                  |
         *    |                  0
         *    |                  2
         *    |                  |
         *  [Ledge]--- 2   0---[mpo]---1
         *    |                  |
         *    |                  2
         *    |
         *    |------ 1
         */

        long                     mpoDimL = mpos.front().dimension(0);
        auto                     Ldim    = Ledge.dimension(0);
        Eigen::Tensor<Scalar, 3> M_temp =
            Ledge
                .shuffle(Textra::array3{0, 2, 1})                                       // Start by shuffling the legs into consecutive order before merge
                .reshape(Textra::array2{Ldim * mpoDimL, Ldim})                          // Merge the legs
                .contract(state.mps_sites.front()->get_M_bare(), Textra::idx({0}, {1})) // Contract with A which already has the mpo on it
                .shuffle(Textra::array3{1, 0, 2});                                      // Shuffle back to convention
        Eigen::Tensor<Scalar, 1> one = Eigen::Tensor<Scalar, 1>(Ldim).constant(1.0);
        state.mps_sites.front()->set_mps(M_temp, one, 0);
    }
    {
        /*
         *     1---[B]---2   0 ------|
         *          |                |
         *          0                |
         *          2                |
         *          |                |
         *    0---[mpo]---1  2 ---[Redge]
         *          |                |
         *          2                |
         *                           |
         *                   1 ------|
         */

        long                     mpoDimR = mpos.back().dimension(1);
        auto                     Rdim    = Redge.dimension(0);
        Eigen::Tensor<Scalar, 3> M_temp  = Redge.shuffle(Textra::array3{0, 2, 1})
                                              .reshape(Textra::array2{Rdim * mpoDimR, Rdim})
                                              .contract(state.mps_sites.back()->get_M_bare(), Textra::idx({0}, {2}))
                                              .shuffle(Textra::array3{1, 2, 0});
        Eigen::Tensor<Scalar, 1> one = Eigen::Tensor<Scalar, 1>(Rdim).constant(1.0);
        state.mps_sites.back()->set_mps(M_temp, one, 0);
    }


    if constexpr(settings::debug){
        state.clear_measurements();
        state.clear_cache();
        state.tag_all_sites_normalized(false); // This operation denormalizes all sites
        tools::log->debug("Norm                 after  applying mpos: {:.16f}", tools::finite::measure::norm(state));
        tools::log->debug("Spin components      after  applying mpos: {}", tools::finite::measure::spin_components(state));
        tools::log->debug("Bond dimensions      after  applying mpos: {}", tools::finite::measure::bond_dimensions(state));
        tools::log->debug("Entanglement entropy after  applying mpos: {}", tools::finite::measure::entanglement_entropies(state));
    }

    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
    state.assert_validity();
}



void tools::finite::ops::project_to_sector(class_state_finite &state, const Eigen::MatrixXcd &paulimatrix, int sign) {
    // This function applies the projection MPO operator  "0.5 * ( 1 - prod s)", where
    // 1 is understood as a 2^L x 2^L tensor and "prod s" is the outer product of pauli matrices, one for each site.
    // This operation leaves the global norm unchanged (thanks to the 0.5 factor) but locally each MPS loses its
    // norm (norm = 2), and entanglement entropies become doubled.
    // Therefore a proper full normalization is required after this operation, as well as a full
    // rebuild of environments.

    if(std::abs(sign) != 1) throw std::runtime_error(fmt::format("Expected 'sign' +1 or -1. Got [{}]", sign));
    tools::common::profile::get_default_prof()["t_prj"]->tic();
    tools::log->debug("Projecting state into sector with sign {}", sign);
    auto spin_components = tools::finite::measure::spin_components(state);
    tools::log->debug("Spin components before projection : X = {:.16f}  Y = {:.16f}  Z = {:.16f}", spin_components[0], spin_components[1], spin_components[2]);
    state.clear_measurements();
    state.clear_cache();
    // Do the projection
    const auto [mpos, L, R] = qm::mpo::parity_projector_mpos(paulimatrix, state.get_length(), sign);
    apply_mpos(state, mpos, L, R);
    spin_components = tools::finite::measure::spin_components(state);
    tools::log->debug("Spin components after  projection : X = {:.16f}  Y = {:.16f}  Z = {:.16f}", spin_components[0], spin_components[1], spin_components[2]);
    tools::common::profile::get_default_prof()["t_prj"]->toc();
}

void tools::finite::ops::project_to_nearest_sector(class_state_finite &state, const std::string &sector) {
    /*
     * When projecting, there is one bad thing that may happen: that the norm of the state vanishes.
     *
     * This can happen in a couple of different scenarios:
     *      - The global state has spin component X = -1.0 and we project to +X (or vice versa)
     *
     * Therefore the projection is only done if the state has a chance of surviving
     * Otherwise emit a warning and return
     *
     */

    tools::log->debug("Projecting state to axis nearest sector {}", sector);
    std::vector<std::string> valid_sectors   = {"x", "+x", "-x", "y", "+y", "-y", "z", "+z", "-z"};
    bool                     sector_is_valid = std::find(valid_sectors.begin(), valid_sectors.end(), sector) != valid_sectors.end();
    auto spin_component_threshold = 1e-3;
    if(sector_is_valid) {
        auto sector_sign = mps::internal::get_sign(sector);
        auto paulimatrix = mps::internal::get_pauli(sector);
        auto spin_component_along_requested_axis = tools::finite::measure::spin_component(state, paulimatrix);
        // Now we have to check that the projection intended projection is safe
        auto alignment = sector_sign * spin_component_along_requested_axis;
        if(alignment > 0)
            // In this case the state has an aligned component along the requested axis --> safe
            project_to_sector(state, paulimatrix, sector_sign);
        else if (alignment < 0){
            // In this case the state has an anti-aligned component along the requested axis --> safe if spin_component < 1 - spin_component_threshold
            if(spin_component_along_requested_axis < 1.0 - spin_component_threshold)
                project_to_sector(state, paulimatrix, sector_sign);
            else return tools::log->warn("Skipping projection to [{}]: State spin component is opposite to the requested projection axis: {:.16f}", sector,spin_component_along_requested_axis);
        }
        else if(alignment == 0){
            // No sector sign was specified, so we select the one along which there is a component
            if(spin_component_along_requested_axis >= spin_component_threshold) sector_sign = 1;
            if(spin_component_along_requested_axis <= spin_component_threshold) sector_sign = -1;
            project_to_sector(state, paulimatrix, sector_sign);
        }
    } else if(sector == "randomAxis") {
        std::vector<std::string> possibilities = {"x", "y", "z"};
        std::string              chosen_axis   = possibilities[rnd::uniform_integer_box<size_t>(0, possibilities.size()-1)];
        project_to_nearest_sector(state, chosen_axis);
    } else if(sector == "random") {
        auto             coeffs    = Eigen::Vector3d::Random().normalized();
        Eigen::Matrix2cd random_c2 = coeffs(0) * qm::spinOneHalf::sx + coeffs(1) * qm::spinOneHalf::sy + coeffs(2) * qm::spinOneHalf::sz;
        return project_to_sector(state, random_c2, 1);
    } else if(sector == "none") {
        return;
    } else
        throw std::runtime_error(fmt::format("Could not parse sector string [{}]", sector));
}

class_state_finite tools::finite::ops::get_projection_to_sector(const class_state_finite &state, const Eigen::MatrixXcd &paulimatrix, int sign) {
    auto state_projected = state;
    project_to_sector(state_projected, paulimatrix, sign);
    return state_projected;
}

class_state_finite tools::finite::ops::get_projection_to_nearest_sector(const class_state_finite &state, const std::string &sector) {
    auto state_projected = state;
    project_to_nearest_sector(state_projected,sector);
    return state_projected;
}

double tools::finite::ops::overlap(const class_state_finite &state1, const class_state_finite &state2) {
    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");
    size_t                   pos     = 0;
    Eigen::Tensor<Scalar, 2> overlap = state1.get_mps_site(pos).get_M().contract(state2.get_mps_site(pos).get_M().conjugate(), Textra::idx({0, 1}, {0, 1}));
    for(pos = 1; pos < state1.get_length(); pos++) {
        Eigen::Tensor<Scalar, 2> temp = overlap.contract(state1.get_mps_site(pos).get_M(), Textra::idx({0}, {1}))
                                            .contract(state2.get_mps_site(pos).get_M().conjugate(), Textra::idx({0, 1}, {1, 0}));
        overlap = temp;
    }

    double norm_chain = std::real(Textra::TensorMatrixMap(overlap).trace());
    return norm_chain;
}

double tools::finite::ops::expectation_value(const class_state_finite &state1, const class_state_finite &state2,
                                             const std::list<Eigen::Tensor<std::complex<double>, 4>> &mpos, const Eigen::Tensor<std::complex<double>, 3> &Ledge,
                                             const Eigen::Tensor<std::complex<double>, 3> &Redge) {
    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");
    auto                     mpo_it = mpos.begin();
    Eigen::Tensor<Scalar, 3> L      = Ledge;
    for(size_t pos = 0; pos < state1.get_length(); pos++) {
        Eigen::Tensor<Scalar, 3> temp = L.contract(state1.get_mps_site(pos).get_M(), idx({0}, {1}))
                                            .contract(*mpo_it++, idx({1, 2}, {0, 2}))
                                            .contract(state2.get_mps_site(pos).get_M().conjugate(), idx({0, 3}, {1, 0}))
                                            .shuffle(array3{0, 2, 1});

        L = temp;
    }
    assert(L.dimensions() == Redge.dimensions());
    Eigen::Tensor<Scalar, 0> E_all_sites  = L.contract(Redge, idx({0, 1, 2}, {0, 1, 2}));
    double                   energy_chain = std::real(E_all_sites(0));
    return energy_chain;
}

double tools::finite::ops::exp_sq_value(const class_state_finite &state1, const class_state_finite &state2,
                                        const std::list<Eigen::Tensor<std::complex<double>, 4>> &mpos, const Eigen::Tensor<std::complex<double>, 4> &Ledge,
                                        const Eigen::Tensor<std::complex<double>, 4> &Redge) {
    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");
    auto                     mpo_it = mpos.begin();
    Eigen::Tensor<Scalar, 4> L      = Ledge;
    for(size_t pos = 0; pos < state1.get_length(); pos++) {
        Eigen::Tensor<Scalar, 4> temp = L.contract(state1.get_mps_site(pos).get_M(), idx({0}, {1}))
                                            .contract(*mpo_it, idx({1, 3}, {0, 2}))
                                            .contract(*mpo_it++, idx({1, 4}, {0, 2}))
                                            .contract(state2.get_mps_site(pos).get_M().conjugate(), idx({0, 4}, {1, 0}))
                                            .shuffle(array4{0, 3, 1, 2});

        L = temp;
    }
    assert(L.dimensions() == Redge.dimensions());
    Eigen::Tensor<Scalar, 0> H2_all_sites = L.contract(Redge, idx({0, 1, 2, 3}, {0, 1, 2, 3}));
    return std::real(H2_all_sites(0));
}
