//
// Created by david on 2020-05-13.
//
#include "env.h"
#include <config/debug.h>
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <math/linalg/tensor.h>
#include <math/num.h>
#include <math/svd.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/edges/class_env_ene.h>
#include <tensors/edges/class_env_var.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/mps.h>

std::vector<size_t> tools::finite::env::expand_subspace(class_state_finite &state, const class_model_finite &model, class_edges_finite &edges, std::optional<double> alpha,
                                                        long chi_lim, std::optional<svd::settings> svd_settings) {
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw std::runtime_error(
            fmt::format("All lengths not equal: state {} | model {} | edges {}", state.get_length(), model.get_length(), edges.get_length()));
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw std::runtime_error(
            fmt::format("All active sites are not equal: state {} | model {} | edges {}", state.active_sites, model.active_sites, edges.active_sites));
    //    if(alpha < 1e-12) return;
    if(state.active_sites.empty()) throw std::runtime_error("No active sites for subspace expansion");
    std::vector<size_t> pos_expanded;

    if(not alpha){
        // When alpha == std::nullopt this turns into a position query
        // Just return the positions that would get modified if alpha had been defined
        auto posL = state.active_sites.front();
        auto posR = state.active_sites.back();
        // Construct PL = alpha * eneL(i-1) * mps(i-1) * mpo(i-1)
        if(posL > 0 and posL < state.get_length<size_t>()) {
            auto &mpsL = state.get_mps_site(posL - 1); // The site to the left
            if(mpsL.get_label() == "A" or mpsL.get_label() == "AC"){
                pos_expanded.emplace_back(posL - 1);
                pos_expanded.emplace_back(posL);
            }
        }
        // Construct PL = alpha * eneL(i+1) * mps(i+1) * mpo(i+1)
        if(posR < state.get_length<size_t>() - 1) {
            auto &mpsR = state.get_mps_site(posR + 1); // The the site to the right
            if(mpsR.get_label() == "B"){
                pos_expanded.emplace_back(posR);
                pos_expanded.emplace_back(posR + 1);
            }
        }
        return pos_expanded;
    }
    using Scalar = class_state_finite::Scalar;

    // Set up the SVD
    svd::solver svd(svd_settings);

    state.clear_cache();

    // Follows the subspace expansion technique explained in https://link.aps.org/doi/10.1103/PhysRevB.91.155115
    {
        auto posL = state.active_sites.front();
        // Construct PL = alpha * eneL(i-1) * mps(i-1) * mpo(i-1)
        if(posL > 0 and posL < state.get_length<size_t>() and state.get_direction() == 1) {
            auto &mpsL = state.get_mps_site(posL - 1); // The site to the left
            auto &mpsR = state.get_mps_site(posL);     // The site to the right
            if(mpsL.get_label() == "A" or mpsL.get_label() == "AC") {
                pos_expanded.emplace_back(posL - 1);
                pos_expanded.emplace_back(posL);
                auto  chiL_old = mpsL.get_chiR();
                auto  dimL_old = mpsL.dimensions();
                auto  dimR_old = mpsR.dimensions();
                auto &mpoL     = model.get_mpo(posL - 1);
                auto &eneL     = edges.get_eneL(posL - 1);
//                auto &varL     = edges.get_varL(posL - 1);

                mpsL.set_M(state.get_multisite_mps({mpsL.get_position()}));
//                Eigen::Tensor<Scalar, 3> PL1      = eneL.get_expansion_term(mpsL, mpoL, alpha.value());
//                Eigen::Tensor<Scalar, 3> PL2      = varL.get_expansion_term(mpsL, mpoL, alpha.value());
//                Eigen::Tensor<Scalar, 3> PL       = PL1.concatenate(PL2, 2);
//                Eigen::Tensor<Scalar, 3> PL = varL.get_expansion_term(mpsL, mpoL, alpha.value());
                Eigen::Tensor<Scalar, 3> PL = eneL.get_expansion_term(mpsL, mpoL, alpha.value());
//                Eigen::Tensor<Scalar, 3> P0 = Textra::TensorConstant<Scalar>(0.0, mpsR.spin_dim(), PL.dimension(2), mpsR.get_chiR());
                Eigen::Tensor<Scalar, 3> P0 = Textra::TensorRandom<double>(mpsR.spin_dim(), PL.dimension(2), mpsR.get_chiR()).cast<Scalar>();
                P0 = P0 * P0.constant(1e-10);
                Eigen::Tensor<Scalar, 3> ML_PL = mpsL.get_M_bare().concatenate(PL, 2);
                Eigen::Tensor<Scalar, 3> MR_P0 = mpsR.get_M_bare().concatenate(P0, 1);

                chi_lim = std::min(chi_lim, mpsL.spin_dim() * std::min(mpsL.get_chiL(), mpsR.get_chiR())); // Bond dimension can't grow faster than x spin_dim
                auto [U, S, V] = svd.schmidt_into_left_normalized(ML_PL, mpsL.spin_dim(), chi_lim);
                mpsL.set_M(U);
                mpsR.set_M(MR_P0);
                if(mpsL.isCenter()) {
                    // Here we expect mpsL to be "AC" and mpsR to be a "B"
                    mpsL.set_LC(S, svd.truncation_error);
                    mpsL.stash_V(V, mpsR.get_position());
                } else {
                    // Here we expect mpsL to be "A" and mpsR to be an "A" or "AC"
                    mpsL.stash_S(S, svd.truncation_error, mpsR.get_position());
                    mpsL.stash_V(V, mpsR.get_position());
                }
                mpsR.merge_stash(mpsL);
                if constexpr(settings::debug) mpsL.assert_identity();

                // Make mpsR normalized
                {
                    Eigen::Tensor<Scalar, 3> M_tmp; // A one-site "theta" i.e. a we want this to have norm 1, but currently it does not, so we rescale mpsR
                    if(mpsR.isCenter() or posL == state.get_length<size_t>() - 1)
                        M_tmp = mpsR.get_M();
                    else if(mpsR.get_label() == "A" and posL < state.get_length<size_t>() - 1) {
                        auto &LR = state.get_mps_site(posL + 1).get_L();
                        M_tmp    = mpsR.get_M().contract(Textra::asDiagonal(LR), Textra::idx({2}, {0}));
                    } else
                        throw std::logic_error(fmt::format("Can't make mpsR normalized: Unknown case: mpsL = {}({}) | mpsR {}({})", mpsL.get_label(),
                                                           mpsL.get_position(), mpsR.get_label(), mpsR.get_position()));

                    double norm_old = Textra::norm(M_tmp.contract(M_tmp.conjugate(), Textra::idx({0, 1, 2}, {0, 1, 2})));
                    M_tmp = mpsR.get_M_bare() * mpsR.get_M_bare().constant(std::pow(norm_old, -0.5)); // Do the rescale (use the tmp object to save memory)
                    mpsR.set_M(M_tmp);
                    if constexpr(settings::debug) {
                        auto   mpsR_final = state.get_multisite_mps({mpsR.get_position()});
                        double norm_new   = Textra::norm(mpsR_final.contract(mpsR_final.conjugate(), Textra::idx({0, 1, 2}, {0, 1, 2})));
                        tools::log->debug("Normalized expanded mps {}({}): {:.16f} -> {:.16f}",mpsR.get_label(), mpsR.get_position(), norm_old, norm_new);
                    }
                }

                auto mpsL_new = state.get_mps_site(posL - 1); // The site to the left
                auto mpsR_new = state.get_mps_site(posL);     // The site to the right
                auto chiL_new = mpsL_new.get_chiR();
                auto dimL_new = mpsL_new.dimensions();
                auto dimR_new = mpsR_new.dimensions();
                tools::log->debug("Subspace expansion pos L {} {} | alpha {:.2e} | χ {} -> {} -> {} | χlim {} ", posL - 1, posL, alpha.value(), chiL_old,
                                 ML_PL.dimension(2), chiL_new, chi_lim);
                if(dimL_old[1] != dimL_new[1])
                    throw std::runtime_error(fmt::format("mpsL changed chiL during left-moving expansion: {} -> {}", dimL_old, dimL_new));
                if(dimR_old[2] != dimR_new[2])
                    throw std::runtime_error(fmt::format("mpsR changed chiR during left-moving expansion: {} -> {}", dimR_old, dimR_new));
            }
        }
    }
    {
        auto posR = state.active_sites.back();
        // Construct PL = alpha * eneL(i+1) * mps(i+1) * mpo(i+1)
        if(posR < state.get_length<size_t>() - 1 and state.get_direction() == -1) {
            auto &mpsR = state.get_mps_site(posR + 1); // The the site to the right
            auto &mpsL = state.get_mps_site(posR);     // The site to the left
            if(mpsR.get_label() == "B") {
                pos_expanded.emplace_back(posR);
                pos_expanded.emplace_back(posR + 1);
                auto  chiR_old = mpsR.get_chiL();
                auto  dimR_old = mpsR.dimensions();
                auto  dimL_old = mpsL.dimensions();
                auto &mpoR     = model.get_mpo(posR + 1);
                auto &eneR     = edges.get_eneR(posR + 1);
//                auto &varR     = edges.get_varR(posR + 1);

                mpsR.set_M(state.get_multisite_mps({mpsR.get_position()}));
//                Eigen::Tensor<Scalar, 3> PR1      = eneR.get_expansion_term(mpsR, mpoR, alpha.value());
//                Eigen::Tensor<Scalar, 3> PR2      = varR.get_expansion_term(mpsR, mpoR, alpha.value());
//                Eigen::Tensor<Scalar, 3> PR       = PR1.concatenate(PR2, 1);
//                Eigen::Tensor<Scalar, 3> PR = varR.get_expansion_term(mpsR, mpoR, alpha.value());
                Eigen::Tensor<Scalar, 3> PR = eneR.get_expansion_term(mpsR, mpoR, alpha.value());
//                Eigen::Tensor<Scalar, 3> P0 = Textra::TensorConstant<Scalar>(0.0, mpsL.spin_dim(), mpsL.get_chiL(), PR.dimension(1));
                Eigen::Tensor<Scalar, 3> P0 = Textra::TensorRandom<double>(mpsL.spin_dim(), mpsL.get_chiL(), PR.dimension(1)).cast<Scalar>();
                P0 = P0 * P0.constant(1e-10);
                Eigen::Tensor<Scalar, 3> MR_PR = mpsR.get_M_bare().concatenate(PR, 1);
                Eigen::Tensor<Scalar, 3> ML_P0 = mpsL.get_M_bare().concatenate(P0, 2); // Usually an AC

                chi_lim = std::min(chi_lim, mpsR.spin_dim() * std::min(mpsL.get_chiL(), mpsR.get_chiR())); // Bond dimension can't grow faster than x spin_dim
                auto [U, S, V] = svd.schmidt_into_right_normalized(MR_PR, mpsR.spin_dim(), chi_lim);
                mpsL.set_M(ML_P0);
                mpsR.set_M(V);
                if(mpsL.isCenter()) {
                    mpsR.stash_C(S, svd.truncation_error, mpsL.get_position());
                    mpsR.stash_U(U, mpsL.get_position());
                } else {
                    // Here we expect mpsL to be a "B" as well
                    mpsR.stash_S(S, svd.truncation_error, mpsL.get_position());
                    mpsR.stash_U(U, mpsL.get_position());
                }
                mpsL.merge_stash(mpsR);
                if constexpr(settings::debug) mpsR.assert_identity();

                // Make mpsL normalized
                {
                    Eigen::Tensor<Scalar, 3> M_tmp; // A one-site "theta" i.e. we want this object to have norm 1, but currently it does not, so we rescale mpsL
                    if(mpsL.isCenter() or posR == 0)
                        M_tmp = mpsL.get_M();
                    else if(mpsL.get_label() == "B" and posR > 0) {
                        auto &mpsLL = state.get_mps_site(posR - 1);
                        if(mpsLL.isCenter())
                            M_tmp = Textra::asDiagonal(mpsLL.get_LC()).contract(mpsL.get_M(), Textra::idx({1}, {1})).shuffle(Textra::array3{1, 0, 2});
                        else
                            M_tmp = Textra::asDiagonal(mpsLL.get_L()).contract(mpsL.get_M(), Textra::idx({1}, {1})).shuffle(Textra::array3{1, 0, 2});
                    } else
                        throw std::logic_error(fmt::format("Can't make mpsL normalized: Unknown case: mpsL = {}({}) | mpsR {}({})", mpsL.get_label(),
                                                           mpsL.get_position(), mpsR.get_label(), mpsR.get_position()));

                    double norm_old = Textra::norm(M_tmp.contract(M_tmp.conjugate(), Textra::idx({0, 1, 2}, {0, 1, 2})));
                    M_tmp = mpsL.get_M_bare() * mpsL.get_M_bare().constant(std::pow(norm_old, -0.5)); // Do the rescale (use the tmp object to save memory)
                    mpsL.set_M(M_tmp);
                    if constexpr(settings::debug) {
                        auto   mpsL_final = state.get_multisite_mps({mpsL.get_position()});
                        double norm_new   = Textra::norm(mpsL_final.contract(mpsL_final.conjugate(), Textra::idx({0, 1, 2}, {0, 1, 2})));
                        tools::log->debug("Normalized expanded mps {}({}): {:.16f} -> {:.16f}",mpsL.get_label(), mpsL.get_position(), norm_old, norm_new);
                    }
                }

                auto mpsR_new = state.get_mps_site(posR + 1);
                auto mpsL_new = state.get_mps_site(posR);
                auto chiR_new = mpsR_new.get_chiL();
                auto dimR_new = mpsR_new.dimensions();
                auto dimL_new = mpsL_new.dimensions();

                tools::log->debug("Subspace expansion pos R {} {} | alpha {:.2e} | χ {} -> {} -> {} | χlim {} ", posR, posR + 1, alpha.value(), chiR_old,
                                 MR_PR.dimension(1), chiR_new, chi_lim);
                if(dimL_old[1] != dimL_new[1])
                    throw std::runtime_error(fmt::format("mpsL changed chiL during right-moving expansion: {} -> {}", dimL_old, dimL_new));
                if(dimR_old[2] != dimR_new[2])
                    throw std::runtime_error(fmt::format("mpsR changed chiR during right-moving expansion: {} -> {}", dimR_old, dimR_new));
            }
        }
    }
    env::rebuild_edges(state, model, edges);
    return pos_expanded;
}

void tools::finite::env::assert_edges_ene(const class_state_finite &state, const class_model_finite &model, const class_edges_finite &edges) {
    size_t min_pos = 0;
    size_t max_pos = state.get_length() - 1;

    // If there are no active sites then we can build up until current position
    // Otherwise it's enough to check up until the bracket defined by active_sites
    long   current_position = state.get_position<long>();
    size_t posL_active      = static_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    size_t posR_active      = static_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    if(not edges.active_sites.empty()) {
        posL_active = edges.active_sites.front();
        posR_active = edges.active_sites.back();
    }

    tools::log->trace("Asserting edges eneL from [{} to {}]", min_pos, posL_active);
    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &ene = edges.get_eneL(pos);
        if(pos == 0 and not ene.has_block()) throw std::runtime_error(fmt::format("ene L at pos {} does not have a block", pos));
        if(pos >= std::min(posL_active, state.get_length() - 1)) continue;
        auto &mps      = state.get_mps_site(pos);
        auto &mpo      = model.get_mpo(pos);
        auto &ene_next = edges.get_eneL(pos + 1);
        ene_next.assert_unique_id(ene, mps, mpo);
    }
    tools::log->trace("Asserting edges eneR from [{} to {}]", posR_active, max_pos);
    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); pos--) {
        auto &ene = edges.get_eneR(pos);
        if(pos == state.get_length() - 1 and not ene.has_block()) throw std::runtime_error(fmt::format("ene R at pos {} does not have a block", pos));
        if(pos <= std::max(posR_active, 0ul)) continue;
        auto &mps      = state.get_mps_site(pos);
        auto &mpo      = model.get_mpo(pos);
        auto &ene_prev = edges.get_eneR(pos - 1);
        ene_prev.assert_unique_id(ene, mps, mpo);
    }
}

void tools::finite::env::assert_edges_var(const class_state_finite &state, const class_model_finite &model, const class_edges_finite &edges) {
    size_t min_pos = 0;
    size_t max_pos = state.get_length() - 1;

    // If there are no active sites then we can build up until current position
    // Otherwise it's enough to check up until the bracket defined by active_sites
    long   current_position = state.get_position<long>();
    size_t posL_active      = static_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    size_t posR_active      = static_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    if(not edges.active_sites.empty()) {
        posL_active = edges.active_sites.front();
        posR_active = edges.active_sites.back();
    }
    tools::log->trace("Asserting edges varL from [{} to {}]", min_pos, posL_active);
    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &var = edges.get_varL(pos);
        if(pos == 0 and not var.has_block()) throw std::runtime_error(fmt::format("var L at pos {} does not have a block", pos));
        if(pos >= std::min(posL_active, state.get_length() - 1)) continue;
        auto &mps      = state.get_mps_site(pos);
        auto &mpo      = model.get_mpo(pos);
        auto &var_next = edges.get_varL(pos + 1);
        var_next.assert_unique_id(var, mps, mpo);
    }
    tools::log->trace("Asserting edges varR from [{} to {}]", posR_active, max_pos);
    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); pos--) {
        auto &var = edges.get_varR(pos);
        if(pos == state.get_length() - 1 and not var.has_block()) throw std::runtime_error(fmt::format("var R at pos {} does not have a block", pos));
        if(pos <= std::max(posR_active, 0ul)) continue;
        auto &mps      = state.get_mps_site(pos);
        auto &mpo      = model.get_mpo(pos);
        auto &var_prev = edges.get_varR(pos - 1);
        var_prev.assert_unique_id(var, mps, mpo);
    }
}

void tools::finite::env::assert_edges(const class_state_finite &state, const class_model_finite &model, const class_edges_finite &edges) {
    assert_edges_ene(state, model, edges);
    assert_edges_var(state, model, edges);
}

void tools::finite::env::rebuild_edges_ene(const class_state_finite &state, const class_model_finite &model, class_edges_finite &edges) {
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw std::runtime_error(
            fmt::format("All lengths not equal: state {} | model {} | edges {}", state.get_length(), model.get_length(), edges.get_length()));
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw std::runtime_error(
            fmt::format("All active sites are not equal: state {} | model {} | edges {}", state.active_sites, model.active_sites, edges.active_sites));
    auto t_env = tools::common::profile::get_default_prof()["t_env"]->tic_token();

    size_t min_pos = 0;
    size_t max_pos = state.get_length() - 1;

    // If there are no active sites then we can build up until current position
    long   current_position = state.get_position<long>();
    size_t posL_active      = static_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    size_t posR_active      = static_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    if(not edges.active_sites.empty()) {
        posL_active = edges.active_sites.front();
        posR_active = edges.active_sites.back();
    }

    tools::log->trace("Inspecting edges eneL from [{} to {}]", min_pos, posL_active);
    std::vector<size_t> ene_pos_log;
    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &ene_curr = edges.get_eneL(pos);
        if(pos == 0 and not ene_curr.has_block()) {
            ene_pos_log.emplace_back(pos);
            ene_curr.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
            if(not ene_curr.has_block()) throw std::runtime_error("No edge detected after setting edge");
        }
        if(pos >= std::min(posL_active, state.get_length() - 1)) continue;
        auto &ene_next = edges.get_eneL(pos + 1);
        auto  id       = ene_next.get_unique_id();
        ene_next.refresh(ene_curr, state.get_mps_site(pos), model.get_mpo(pos));
        if(id != ene_next.get_unique_id()) ene_pos_log.emplace_back(ene_next.get_position());
    }
    if(not ene_pos_log.empty()) tools::log->trace("Rebuilt L ene edges: {}", ene_pos_log);
    ene_pos_log.clear();
    tools::log->trace("Inspecting edges eneR from [{} to {}]", posR_active, max_pos);
    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); pos--) {
        auto &ene_curr = edges.get_eneR(pos);
        if(pos == state.get_length() - 1 and not ene_curr.has_block()) {
            ene_pos_log.emplace_back(pos);
            ene_curr.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
            if(not ene_curr.has_block()) throw std::runtime_error("No edge detected after setting edge");
        }
        if(pos <= std::max(posR_active, 0ul)) continue;
        auto &ene_prev = edges.get_eneR(pos - 1);
        auto  id       = ene_prev.get_unique_id();
        ene_prev.refresh(ene_curr, state.get_mps_site(pos), model.get_mpo(pos));
        if(id != ene_prev.get_unique_id()) ene_pos_log.emplace_back(ene_prev.get_position());
    }
    std::reverse(ene_pos_log.begin(), ene_pos_log.end());
    if(not ene_pos_log.empty()) tools::log->trace("Rebuilt R ene edges: {}", ene_pos_log);
    if(not edges.get_eneL(posL_active).has_block()) throw std::logic_error(fmt::format("Left active ene edge has undefined block"));
    if(not edges.get_eneR(posR_active).has_block()) throw std::logic_error(fmt::format("Right active ene edge has undefined block"));
}

void tools::finite::env::rebuild_edges_var(const class_state_finite &state, const class_model_finite &model, class_edges_finite &edges) {
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw std::runtime_error(
            fmt::format("All lengths not equal: state {} | model {} | edges {}", state.get_length(), model.get_length(), edges.get_length()));
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw std::runtime_error(
            fmt::format("All active sites are not equal: state {} | model {} | edges {}", state.active_sites, model.active_sites, edges.active_sites));
    auto t_env = tools::common::profile::get_default_prof()["t_env"]->tic_token();

    size_t min_pos = 0;
    size_t max_pos = state.get_length() - 1;

    // If there are no active sites then we can build up until current position
    // Otherwise it's enough to build up until the bracket defined by active_sites
    long   current_position = state.get_position<long>();
    size_t posL_active      = static_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    size_t posR_active      = static_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    if(not edges.active_sites.empty()) {
        posL_active = edges.active_sites.front();
        posR_active = edges.active_sites.back();
    }
    tools::log->trace("Inspecting edges varL from [{} to {}]", min_pos, posL_active);
    std::vector<size_t> var_pos_log;
    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &var_curr = edges.get_varL(pos);
        if(pos == 0 and not var_curr.has_block()) {
            var_pos_log.emplace_back(pos);
            var_curr.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
            if(not var_curr.has_block()) throw std::runtime_error("No edge detected after setting edge");
        }
        if(pos >= std::min(posL_active, state.get_length() - 1)) continue;
        auto &var_next = edges.get_varL(pos + 1);
        auto  id       = var_next.get_unique_id();
        var_next.refresh(var_curr, state.get_mps_site(pos), model.get_mpo(pos));
        if(id != var_next.get_unique_id()) var_pos_log.emplace_back(var_next.get_position());
    }

    if(not var_pos_log.empty()) tools::log->trace("Rebuilt L var edges: {}", var_pos_log);
    var_pos_log.clear();
    tools::log->trace("Inspecting edges varR from [{} to {}]", posR_active, max_pos);
    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); pos--) {
        auto &var_curr = edges.get_varR(pos);
        if(pos == state.get_length() - 1 and not var_curr.has_block()) {
            var_pos_log.emplace_back(pos);
            var_curr.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
            if(not var_curr.has_block()) throw std::runtime_error("No edge detected after setting edge");
        }
        if(pos <= std::max(posR_active, 0ul)) continue;
        auto &var_prev = edges.get_varR(pos - 1);
        auto  id       = var_prev.get_unique_id();
        var_prev.refresh(var_curr, state.get_mps_site(pos), model.get_mpo(pos));
        if(id != var_prev.get_unique_id()) var_pos_log.emplace_back(var_prev.get_position());
    }
    std::reverse(var_pos_log.begin(), var_pos_log.end());
    if(not var_pos_log.empty()) tools::log->trace("Rebuilt R var edges: {}", var_pos_log);
    if(not edges.get_varL(posL_active).has_block()) throw std::logic_error(fmt::format("Left active var edge has undefined block"));
    if(not edges.get_varR(posR_active).has_block()) throw std::logic_error(fmt::format("Right active var edge has undefined block"));
}

void tools::finite::env::rebuild_edges(const class_state_finite &state, const class_model_finite &model, class_edges_finite &edges) {
    rebuild_edges_ene(state, model, edges);
    rebuild_edges_var(state, model, edges);
}
