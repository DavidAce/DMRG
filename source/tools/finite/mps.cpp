//
// Created by david on 2019-01-29.
//

#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include <tools/finite/measure.h>
#include <config/enums.h>
#include <config/nmspc_settings.h>
#include <general/nmspc_iter.h>
#include <iostream>
#include <utility>
#include <math/num.h>
#include <physics/nmspc_quantum_mechanics.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/common/split.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>

bool tools::finite::mps::internal::bitfield_is_valid(std::optional<long> bitfield) {
    return bitfield.has_value() and bitfield.value() > 0 and internal::used_bitfields.count(bitfield.value()) == 0;
}

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
        // Do the same with its truncation error
        Eigen::Tensor<Scalar, 1> LC                  = mps.get_LC();
        double                   truncation_error_LC = mps.get_truncation_error_LC();
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
        tools::finite::mps::merge_multisite_tensor(state, twosite_tensor, {posL, posR}, posL, chi_lim, svd_threshold, LogPolicy::QUIET);
        state.clear_cache(LogPolicy::QUIET);
        state.clear_measurements(LogPolicy::QUIET);

        // Put LC where it belongs.
        // Recall that mpsL, mpsR are on the new position, not the old one!
        if(state.get_direction() == 1) mpsL.set_L(LC, truncation_error_LC);
        else
            mpsR.set_L(LC, truncation_error_LC);
    }
}

void tools::finite::mps::move_center_point_to_edge(class_state_finite &state, long chi_lim, std::optional<double> svd_threshold){
    while(not state.position_is_any_edge())
        move_center_point(state,chi_lim,svd_threshold);
    state.flip_direction();
}

void tools::finite::mps::move_center_point_to_middle(class_state_finite &state, long chi_lim, std::optional<double> svd_threshold){
    while(not state.position_is_the_middle())
        move_center_point(state,chi_lim,svd_threshold);
}

void tools::finite::mps::merge_multisite_tensor(class_state_finite &state, const Eigen::Tensor<Scalar, 3> &multisite_mps, const std::vector<size_t> &sites,
                                                size_t center_position, long chi_lim, std::optional<double> svd_threshold, std::optional<LogPolicy> logPolicy) {
    if (not logPolicy or logPolicy == LogPolicy::NORMAL) tools::log->trace("Merging multisite tensor for sites {} | chi limit {}", sites, chi_lim);
    // Some sanity checks
    if(multisite_mps.dimension(1) != state.get_mps_site(sites.front()).get_chiL())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim1 {} != chiL on left-most site {}", multisite_mps.dimension(1),
                                             state.get_mps_site(sites.front()).get_chiL(), sites.front()));

    if(multisite_mps.dimension(2) != state.get_mps_site(sites.back()).get_chiR())
        throw std::runtime_error(fmt::format("Could not merge multisite mps into state: mps dim2 {} != chiR on right-most site {}", multisite_mps.dimension(2),
                                             state.get_mps_site(sites.back()).get_chiR(), sites.back()));
    long              spin_prod = 1;
    std::vector<long> spin_dims;
    for(const auto &site : sites) {
        spin_dims.emplace_back(state.get_mps_site(site).spin_dim());
        spin_prod *= spin_dims.back();
    }
    if(spin_prod != multisite_mps.dimension(0))
        throw std::runtime_error(
            fmt::format("Could not merge multisite mps into state: multisite_mps dim0 {} != spin_prod {}", multisite_mps.dimension(0), spin_prod));

    // Split the multisite mps into single-site mps objects
    auto mps_list = tools::common::split::split_mps(multisite_mps, spin_dims, sites, center_position, chi_lim, svd_threshold);

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
        state.tag_site_normalized(mps_tgt.get_position(), true); // Merged site is normalized
        mps_ptr++;
    }
    state.clear_cache(LogPolicy::QUIET);
    state.clear_measurements(LogPolicy::QUIET);
}

bool tools::finite::mps::normalize_state(class_state_finite &state, long chi_lim, std::optional<double> svd_threshold, NormPolicy norm_policy) {
    // When a state needs to be normalized it's enough to "move" the center position around the whole chain.
    // Each move performs an SVD decomposition which leaves unitaries after it, effectively normalizing the state.
    // NOTE! It may be important to start with the current position.

    if(norm_policy == NormPolicy::IFNEEDED) {
        // We may only go ahead with a normalization if its really needed.
        if(std::abs(tools::finite::measure::norm(state) - 1.0) < settings::precision::max_norm_error) return false;
    }

    // Otherwise we just do the normalization
    if(tools::Logger::getLogLevel(tools::log) <= 0 )
        tools::log->trace("Normalizing state | Old norm = {:.16f}",tools::finite::measure::norm(state));

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
    tools::finite::mps::merge_multisite_tensor(state, twosite_tensor, {posL, posR}, posL, chi_lim, svd_threshold,LogPolicy::QUIET);

    // Now we can move around the chain
    for(size_t move = 0; move < num_moves; move++) move_center_point(state, chi_lim, svd_threshold);
    state.clear_measurements();
    state.clear_cache();
    if(tools::Logger::getLogLevel(tools::log) <= 0 )
        tools::log->trace("Normalizing state | New norm = {:.16f}",tools::finite::measure::norm(state));
    state.assert_validity();
    return true;
}

void tools::finite::mps::randomize_state(class_state_finite &state, StateInit init, StateInitType type, const std::string &sector, long chi_lim, bool use_eigenspinors,
                                         std::optional<long> bitfield) {
    switch(init) {
        case StateInit::RANDOM_PRODUCT_STATE: return internal::random_product_state(state, type, sector, use_eigenspinors, bitfield);
        case StateInit::RANDOM_ENTANGLED_STATE: return internal::random_entangled_state(state,type, sector, chi_lim, use_eigenspinors);
        case StateInit::RANDOMIZE_PREVIOUS_STATE: return internal::randomize_given_state(state,type);
        case StateInit::PRODUCT_STATE: return internal::set_product_state(state, type, sector);
    }
}

void tools::finite::mps::apply_random_paulis(class_state_finite &state, const std::vector<Eigen::Matrix2cd> & paulimatrices) {
    auto [mpos, L, R] = qm::mpo::sum_of_pauli_mpo(paulimatrices, state.get_length(),true);
    tools::finite::ops::apply_mpos(state, mpos, L, R);
}

void tools::finite::mps::apply_random_paulis(class_state_finite &state, const std::vector<std::string> &paulistrings) {
    std::vector<Eigen::Matrix2cd> paulimatrices;
    for(const auto &str : paulistrings) paulimatrices.emplace_back(internal::get_pauli(str));
    apply_random_paulis(state,paulimatrices);
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


void is_unitary(const Eigen::Tensor<std::complex<double>, 2> & U_tensor){
    Eigen::Tensor<std::complex<double>, 2> UU = U_tensor.contract(U_tensor.conjugate(), Textra::idx({0},{0}));
    std::cout << UU << std::endl;
}


void tools::finite::mps::apply_twosite_operators(class_state_finite &state, const std::vector<Eigen::Tensor<Scalar, 2>> &twosite_operators, long chi_lim, std::optional<double> svd_threshold) {
    if(twosite_operators.size()+1 != state.get_length())
        throw std::runtime_error(fmt::format("Size mismatch: Given {} two-site operators for a system of {} sites. Expected {}", twosite_operators.size(), state.get_length(),state.get_length()-1));
    Eigen::IOFormat  CleanFmt(4, 0, ", ", "\n", "  [", "]");
    tools::log->info("Unitary gates");

    Eigen::Tensor<Scalar,3>   zminus_spinor  = Textra::MatrixToTensor(internal::get_spinor("z", -1).normalized(), 2, 1, 1);
    state.get_mps_site(2).set_M(zminus_spinor);

    for(size_t idx = 0; idx < twosite_operators.size(); idx++){
        std::cout << "U(" << idx << ")\n" << Textra::TensorMatrixMap(twosite_operators[idx]).format(CleanFmt) << std::endl;
    }


    tools::log->info("Before applying unitaries");
    for(auto && mps : state.mps_sites)
        std::cout << "M(" << mps->get_position() << ") dims ["<< mps->spin_dim() << "," << mps->get_chiL() << "," << mps->get_chiR() << "]:\n" << Textra::TensorMatrixMap(mps->get_M_bare(),mps->spin_dim(),mps->get_chiL()*mps->get_chiR()).format(CleanFmt) << std::endl;

    while (state.get_position() != 0 or state.get_direction() != 1) {
        move_center_point(state, chi_lim, svd_threshold);
        tools::log->debug("At pos {} dir {}",state.get_position(), state.get_direction());
    }
    tools::log->debug("Starting even at pos {} dir {}",state.get_position(), state.get_direction());
    for(size_t idx = 0; idx < twosite_operators.size(); idx+=2){
        // Start applying from the end point
        if(idx != state.get_position()) throw std::logic_error(fmt::format("idx and position mismatch {} != {}",idx, state.get_position()));
        size_t                   posL      = state.get_position();
        size_t                   posR      = posL + 1;
        auto &                   mpsL      = state.get_mps_site(posL);
        auto &                   mpsR      = state.get_mps_site(posR);
        long                     dL        = mpsL.spin_dim();
        long                     dR        = mpsR.spin_dim();
        long                     chiL      = mpsL.get_chiL();
        long                     chiR      = mpsR.get_chiR();
        tools::log->debug("Applying even operator at index {} | sites [{} {}]",idx,posL,posR);
        Eigen::Tensor<Scalar, 3> twosite_tensor =
            mpsL.get_M().contract(mpsR.get_M(), Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{dL * dR, chiL, chiR});
        Eigen::Tensor<Scalar,3> twosite_tensor_op = twosite_operators[idx].contract(twosite_tensor, Textra::idx({0},{0}));
        tools::finite::mps::merge_multisite_tensor(state, twosite_tensor_op, {posL, posR}, posL, chi_lim, svd_threshold);
        move_center_point(state, chi_lim, svd_threshold);
        move_center_point(state, chi_lim, svd_threshold);
    }
    tools::log->debug("Finished even at pos {} dir {}",state.get_position(), state.get_direction());
    while (state.get_position() != 1 or state.get_direction() != 1) {
        move_center_point(state, chi_lim, svd_threshold);
        tools::log->debug("At pos {} dir {}",state.get_position(), state.get_direction());
    }
    tools::log->debug("Starting odd at pos {} dir {}",state.get_position(), state.get_direction());
    for(size_t idx = 1; idx < twosite_operators.size(); idx+=2){
        // Start applying from the end point
        if(idx != state.get_position()) throw std::logic_error(fmt::format("idx and position mismatch {} != {}",idx, state.get_position()));
        size_t                   posL      = state.get_position();
        size_t                   posR      = posL + 1;
        auto &                   mpsL      = state.get_mps_site(posL);
        auto &                   mpsR      = state.get_mps_site(posR);
        long                     dL        = mpsL.spin_dim();
        long                     dR        = mpsR.spin_dim();
        long                     chiL      = mpsL.get_chiL();
        long                     chiR      = mpsR.get_chiR();
        tools::log->debug("Applying odd operator at index {} | sites [{} {}]",idx,posL,posR);

        Eigen::Tensor<Scalar, 3> twosite_tensor =
            mpsL.get_M().contract(mpsR.get_M(), Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{dL * dR, chiL, chiR});
        Eigen::Tensor<Scalar,3> twosite_tensor_op = twosite_operators[idx].contract(twosite_tensor, Textra::idx({0},{0}));
        tools::finite::mps::merge_multisite_tensor(state, twosite_tensor_op, {posL, posR}, posL, chi_lim, svd_threshold);
        move_center_point(state, chi_lim, svd_threshold);
        move_center_point(state, chi_lim, svd_threshold);
    }
    while (state.get_position() != 0 or state.get_direction() != 1)
        move_center_point(state, chi_lim, svd_threshold);
    state.clear_measurements();
    state.clear_cache();
    state.assert_validity();


    tools::log->info("After applying unitaries");

    for(auto && mps : state.mps_sites)
        std::cout << "M(" << mps->get_position() << ") dims ["<< mps->spin_dim() << "," << mps->get_chiL() << "," << mps->get_chiR() << "]:\n" << Textra::TensorMatrixMap(mps->get_M_bare(),mps->spin_dim(),mps->get_chiL()*mps->get_chiR()).format(CleanFmt) << std::endl;

    exit(0);

}
