#include "../env.h"
#include "config/debug.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "math/num.h"
#include "math/svd.h"
#include "math/tenx.h"
#include "tensors/edges/EdgesFinite.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/site/env/EnvEne.h"
#include "tensors/site/env/EnvVar.h"
#include "tensors/site/mpo/MpoSite.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"

namespace settings {
    static constexpr bool debug_edges     = false;
    static constexpr bool debug_expansion = false;
}

std::vector<size_t> tools::finite::env::expand_environment(StateFinite &state, const ModelFinite &model, EdgesFinite &edges, std::optional<double> alpha,
                                                           EnvExpandMode envExpandMode, std::optional<svd::config> svd_cfg) {
    // Follows the environment expansion technique explained in https://link.aps.org/doi/10.1103/PhysRevB.91.155115

    // Going left-to-right with active site == i,  we would update:
    //      A(i-1) --> [A(i-1), P(i-1)],
    //      AC(i)  --> [AC(i) , 0 ]^T  <--- This one loses normalization, but it doesn't matter because it will get optimized
    // where
    //      P(i-1) = alpha * envL(i) * mps(i) * mpo(i) (dimensions d, chi(i-1), chi(i)*m(i))

    // Going right-to-left with active site == i,  we would update:
    //      AC(i)   --> [AC(i), 0 ]     <--- This one loses normalization, but it doesn't matter because it will get optimized
    //      B(i+1)  --> [B(i+1), P(i+1)]^T,
    // where
    //      P(i) = alpha * envR(i) * mps(i) * mpo(i) (dimensions d, chi(i-1)*m(i), chi(i))
    //
    // Note that in this convention we expand the "trailing" bond during a DMRG sweep. Doing it the other way would
    // force us to insert a non-normalized mps into the environment ahead.
    // The rule is to optimize the site which loses normalization during the expansion.

    using Scalar = cplx;
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw except::runtime_error("All lengths not equal: state {} | model {} | edges {}", state.get_length(), model.get_length(), edges.get_length());
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw except::runtime_error("All active sites are not equal: state {} | model {} | edges {}", state.active_sites, model.active_sites,
                                    edges.active_sites);
    if(state.active_sites.empty()) throw except::logic_error("No active sites for environment expansion");
    if(state.active_sites.size() != 1) {
        tools::log->warn("expand_environment_var: expected active_sites.size() == 1. Got: {}", state.active_sites);
        return {};
    }
    if(state.active_sites.front() != state.get_position()) {
        tools::log->warn("expand_environment_var: expected active_sites.front() == current position. Got: {} != {}", state.active_sites, state.get_position());
        return {};
    }
    auto                pos = state.active_sites.front();
    std::vector<size_t> pos_expanded;
    if(state.get_direction() > 0 and pos == std::clamp<size_t>(pos, 0, state.get_length<size_t>() - 2)) pos_expanded = {pos, pos + 1};
    if(state.get_direction() < 0 and pos == std::clamp<size_t>(pos, 1, state.get_length<size_t>() - 1)) pos_expanded = {pos - 1, pos};
    if(not alpha or pos_expanded.empty()) {
        // When alpha == std::nullopt this turns into a position query
        // Just return the positions that would get modified if alpha had been defined.
        return pos_expanded;
    }

    // Set up the SVD
    svd::solver svd(svd_cfg);

    // Define the left and right mps that will get modified
    state.clear_cache();
    auto  posL     = pos_expanded.front();
    auto  posR     = pos_expanded.back();
    auto &mpsL     = state.get_mps_site(posL);
    auto &mpsR     = state.get_mps_site(posR);
    auto  dimL_old = mpsL.dimensions();
    auto  dimR_old = mpsR.dimensions();

    // Set up the svd
    auto bond_lim = mpsR.spin_dim() * std::min(mpsL.get_chiL(), mpsR.get_chiR()); // Bond dimension can't grow faster than x spin_dim
    if(not svd_cfg) svd_cfg = svd::config();
    if(not svd_cfg->rank_max) svd_cfg->rank_max = bond_lim;
    svd_cfg->rank_max = std::min(bond_lim, svd_cfg->rank_max.value());

    if(state.get_direction() > 0) {
        // The expanded bond sits between mpsL and mpsR. When direction is left-to-right:
        //      * mpsL is an "AC" and is the active site
        //      * mpsR is an "B"  and belongs in envR later in the optimization step
        auto &mpoR = model.get_mpo(posR);
        mpsR.set_M(state.get_multisite_mps({posR})); // mpsR absorbs the C bond so that the SVD makes sense later
        Eigen::Tensor<Scalar, 3> PR;
        switch(envExpandMode) {
            case EnvExpandMode::ENE: PR = edges.get_env_eneR(posR).get_expansion_term(mpsR, mpoR, alpha.value()); break;
            case EnvExpandMode::VAR: PR = edges.get_env_varR(posR).get_expansion_term(mpsR, mpoR, alpha.value()); break;
        }
        Eigen::Tensor<Scalar, 3> P0    = tenx::TensorConstant<Scalar>(0.0, mpsL.spin_dim(), mpsL.get_chiL(), PR.dimension(1));
        Eigen::Tensor<Scalar, 3> ML_P0 = mpsL.get_M_bare().concatenate(P0, 2); // mpsL is going to be optimized, zero-padded with P0
        Eigen::Tensor<Scalar, 3> MR_PR = mpsR.get_M_bare().concatenate(PR, 1); // mpsR is going into the environment, enriched with PR.

        auto [U, S, V] = svd.schmidt_into_right_normalized(MR_PR, mpsR.spin_dim(), svd_cfg.value());
        mpsR.set_M(V);
        mpsR.stash_U(U, posL);
        mpsR.stash_C(S, -1.0, posL); // Set a negative truncation error to ignore it.
        mpsL.set_M(ML_P0);
        mpsL.take_stash(mpsR); // normalization of mpsL is lost here.

        {
            // Make mpsL normalized. This is strictly optional, but we do it so that normalization checks can succeed
            double                   norm_old = tools::common::contraction::contract_mps_norm(mpsL.get_M());
            Eigen::Tensor<Scalar, 3> M_tmp    = mpsL.get_M_bare() * mpsL.get_M_bare().constant(std::pow(norm_old, -0.5)); // Rescale
            mpsL.set_M(M_tmp);
            if constexpr(settings::debug or settings::debug_expansion) {
                auto   mpsL_final = state.get_multisite_mps({mpsL.get_position()});
                double norm_new   = tools::common::contraction::contract_mps_norm(mpsL_final);
                tools::log->debug("Normalized expanded mps {}({}): {:.16f} -> {:.16f}", mpsL.get_label(), mpsL.get_position(), norm_old, norm_new);
            }
        }

        tools::log->debug("Environment expansion pos {} | alpha {:.2e} | χ {} -> {} -> {} | χlim {}", pos_expanded, alpha.value(), dimL_old[2],
                          MR_PR.dimension(1), mpsL.get_chiR(), bond_lim);
    }
    if(state.get_direction() < 0) {
        // The expanded bond sits between mpsL and mpsR. When direction is left-to-right:
        //      * mpsL is an "A" and belongs in envL later in the optimization step
        //      * mpsR is an "AC" and is the active site
        auto &mpoL = model.get_mpo(posL);
        mpsL.set_M(state.get_multisite_mps({posL})); // mpsL absorbs the C bond so that the SVD makes sense later
        Eigen::Tensor<Scalar, 3> PL;
        switch(envExpandMode) {
            case EnvExpandMode::ENE: PL = edges.get_env_eneL(posL).get_expansion_term(mpsL, mpoL, alpha.value()); break;
            case EnvExpandMode::VAR: PL = edges.get_env_varL(posL).get_expansion_term(mpsL, mpoL, alpha.value()); break;
        }
        Eigen::Tensor<Scalar, 3> P0    = tenx::TensorConstant<Scalar>(0.0, mpsR.spin_dim(), PL.dimension(2), mpsR.get_chiR());
        Eigen::Tensor<Scalar, 3> ML_PL = mpsL.get_M_bare().concatenate(PL, 2);
        Eigen::Tensor<Scalar, 3> MR_P0 = mpsR.get_M_bare().concatenate(P0, 1);

        auto [U, S, V] = svd.schmidt_into_left_normalized(ML_PL, mpsL.spin_dim(), svd_cfg.value());
        // Recall: the following lines hold because we demand that mpsL is an "A" and mpsR an "AC"
        mpsL.set_M(U);
        mpsL.stash_S(S, -1.0, posR); // Set a negative truncation error to ignore it.
        mpsL.stash_V(V, posR);
        mpsR.set_M(MR_P0);
        mpsR.take_stash(mpsL); // normalization of mpsR is lost here

        {
            // Make mpsR normalized. This is strictly optional, but we do it so that normalization checks can succeed
            double                   norm_old = tools::common::contraction::contract_mps_norm(mpsR.get_M());
            Eigen::Tensor<Scalar, 3> M_tmp    = mpsR.get_M_bare() * mpsR.get_M_bare().constant(std::pow(norm_old, -0.5)); // Rescale
            mpsR.set_M(M_tmp);
            if constexpr(settings::debug or settings::debug_expansion) {
                auto   mpsR_final = state.get_multisite_mps({mpsR.get_position()});
                double norm_new   = tools::common::contraction::contract_mps_norm(mpsR_final);
                tools::log->debug("Normalized expanded mps {}({}): {:.16f} -> {:.16f}", mpsR.get_label(), mpsR.get_position(), norm_old, norm_new);
            }
        }
        tools::log->debug("Environment expansion pos {} | alpha {:.2e} | χ {} -> {} -> {} | χlim {}", pos_expanded, alpha.value(), dimR_old[1],
                          ML_PL.dimension(2), mpsR.get_chiL(), bond_lim);
    }

    if(dimL_old[1] != mpsL.get_chiL()) throw except::runtime_error("mpsL changed chiL during environment expansion: {} -> {}", dimL_old, mpsL.dimensions());
    if(dimR_old[2] != mpsR.get_chiR()) throw except::runtime_error("mpsR changed chiR during environment expansion: {} -> {}", dimR_old, mpsR.dimensions());
    if constexpr(settings::debug or settings::debug_expansion) mpsL.assert_normalized();
    if constexpr(settings::debug or settings::debug_expansion) mpsR.assert_normalized();
    env::rebuild_edges(state, model, edges);
    return pos_expanded;
}

void tools::finite::env::assert_edges_ene(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges) {
    if(state.get_algorithm() == AlgorithmType::fLBIT) throw except::logic_error("assert_edges_var: fLBIT algorithm should never assert energy edges!");
    size_t min_pos = 0;
    size_t max_pos = state.get_length() - 1;

    // If there are no active sites we shouldn't be asserting edges.
    // For instance, the active sites are cleared after a move of center site.
    // We could always keep all edges refreshed but that would be wasteful, since the next iteration
    // may activate other sites and not end up needing those edges.
    // Instead, we force the hand of the algorithm, to only allow edge assertions with active sites defined.
    // Ideally, then, this should be done directly after activating new sites in a new iteration.
    if(edges.active_sites.empty())
        throw except::runtime_error("assert_edges_ene: no active sites.\n"
                                    "Hint:\n"
                                    " One could in principle keep edges refreshed always, but\n"
                                    " that would imply rebuilding many edges that end up not\n"
                                    " being used. Make sure to only run this assertion after\n"
                                    " activating sites.");

    long   current_position = state.get_position<long>();
    size_t posL_active      = edges.active_sites.front();
    size_t posR_active      = edges.active_sites.back();

    //    size_t posL_active      = static_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    //    size_t posR_active      = static_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    //    if(not edges.active_sites.empty()) {
    //        posL_active = edges.active_sites.front();
    //        posR_active = edges.active_sites.back();
    //    }

    if constexpr(settings::debug or settings::debug_edges)
        tools::log->trace("assert_edges_ene: pos {} | dir {} | "
                          "asserting edges eneL from [{} to {}]",
                          current_position, state.get_direction(), min_pos, posL_active);

    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &ene = edges.get_env_eneL(pos);
        if(pos == 0 and not ene.has_block()) throw except::runtime_error("ene L at pos {} does not have a block", pos);
        if(pos >= std::min(posL_active, state.get_length() - 1)) continue;
        auto &mps      = state.get_mps_site(pos);
        auto &mpo      = model.get_mpo(pos);
        auto &ene_next = edges.get_env_eneL(pos + 1);
        ene_next.assert_unique_id(ene, mps, mpo);
    }
    if constexpr(settings::debug or settings::debug_edges)
        tools::log->trace("assert_edges_ene: pos {} | dir {} | "
                          "asserting edges eneR from [{} to {}]",
                          current_position, state.get_direction(), posR_active, max_pos);

    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); --pos) {
        auto &ene = edges.get_env_eneR(pos);
        if(pos == state.get_length() - 1 and not ene.has_block()) throw except::runtime_error("ene R at pos {} does not have a block", pos);
        if(pos <= std::max(posR_active, 0ul)) continue;
        auto &mps      = state.get_mps_site(pos);
        auto &mpo      = model.get_mpo(pos);
        auto &ene_prev = edges.get_env_eneR(pos - 1);
        ene_prev.assert_unique_id(ene, mps, mpo);
    }
}

void tools::finite::env::assert_edges_var(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges) {
    if(state.get_algorithm() == AlgorithmType::fLBIT) throw except::logic_error("assert_edges_var: fLBIT algorithm should never assert variance edges!");
    size_t min_pos = 0;
    size_t max_pos = state.get_length() - 1;

    // If there are no active sites we shouldn't be asserting edges.
    // For instance, the active sites are cleared after a move of center site.
    // We could always keep all edges refreshed but that would be wasteful, since the next iteration
    // may activate other sites and not end up needing those edges.
    // Instead, we force the hand of the algorithm, to only allow edge assertions with active sites defined.
    // Ideally, then, this should be done directly after activating new sites in a new iteration.
    if(edges.active_sites.empty())
        throw except::runtime_error("assert_edges_var: no active sites.\n"
                                    "Hint:\n"
                                    " One could in principle keep edges refreshed always, but\n"
                                    " that would imply rebuilding many edges that end up not\n"
                                    " being used. Make sure to only run this assertion after\n"
                                    " activating sites.");

    long   current_position = state.get_position<long>();
    size_t posL_active      = edges.active_sites.front();
    size_t posR_active      = edges.active_sites.back();

    //    long   current_position = state.get_position<long>();
    //    size_t posL_active      = static_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    //    size_t posR_active      = static_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    //    if(not edges.active_sites.empty()) {
    //        posL_active = edges.active_sites.front();
    //        posR_active = edges.active_sites.back();
    //    }
    if constexpr(settings::debug or settings::debug_edges)
        tools::log->trace("assert_edges_var: pos {} | dir {} | "
                          "asserting edges varL from [{} to {}]",
                          current_position, state.get_direction(), min_pos, posL_active);
    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &var = edges.get_env_varL(pos);
        if(pos == 0 and not var.has_block()) throw except::runtime_error("var L at pos {} does not have a block", pos);
        if(pos >= std::min(posL_active, state.get_length() - 1)) continue;
        auto &mps      = state.get_mps_site(pos);
        auto &mpo      = model.get_mpo(pos);
        auto &var_next = edges.get_env_varL(pos + 1);
        var_next.assert_unique_id(var, mps, mpo);
    }
    if constexpr(settings::debug or settings::debug_edges)
        tools::log->trace("assert_edges_var: pos {} | dir {} | "
                          "asserting edges varR from [{} to {}]",
                          current_position, state.get_direction(), posR_active, max_pos);
    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); --pos) {
        auto &var = edges.get_env_varR(pos);
        if(pos == state.get_length() - 1 and not var.has_block()) throw except::runtime_error("var R at pos {} does not have a block", pos);
        if(pos <= std::max(posR_active, 0ul)) continue;
        auto &mps      = state.get_mps_site(pos);
        auto &mpo      = model.get_mpo(pos);
        auto &var_prev = edges.get_env_varR(pos - 1);
        var_prev.assert_unique_id(var, mps, mpo);
    }
}

void tools::finite::env::assert_edges(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges) {
    if(state.get_algorithm() == AlgorithmType::fLBIT) return;
    assert_edges_ene(state, model, edges);
    assert_edges_var(state, model, edges);
}

void tools::finite::env::rebuild_edges_ene(const StateFinite &state, const ModelFinite &model, EdgesFinite &edges) {
    if(state.get_algorithm() == AlgorithmType::fLBIT) throw except::logic_error("rebuild_edges_ene: fLBIT algorithm should never rebuild energy edges!");
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw except::runtime_error("All lengths not equal: state {} | model {} | edges {}", state.get_length(), model.get_length(), edges.get_length());
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw except::runtime_error("All active sites are not equal: state {} | model {} | edges {}", state.active_sites, model.active_sites,
                                    edges.active_sites);
    auto   t_reb   = tid::tic_scope("rebuild_edges_ene", tid::higher);
    size_t min_pos = 0;
    size_t max_pos = state.get_length() - 1;

    // If there are no active sites then we can build up until current position
    /*
     * LOG:
     * - 2021-10-14:
     *      Just had a terribly annoying bug:
     *      Moving the center position clears active_sites, which caused problems when turning back from the right edge.
     *          1)  active_sites [A(L-1), AC(L)] are updated, left edge exist for A(L-1), right edge exists for AC(L)
     *          2)  move dir -1, clear active sites
     *          3)  assert_edges checks up to AC(L-1), but this site has a stale right edge.
     *      Therefore, one would have to rebuild edges between steps 2) and 3) to solve this issue
     *
     *      One solution would be to always rebuild edges up to the current position from both sides, but that would be
     *      wasteful. Instead, we could just accept that some edges are stale after moving the center-point,
     *      as long as we rebuild those when sites get activated again.
     *
     */

    // If there are no active sites we shouldn't be rebuilding edges.
    // For instance, the active sites are cleared after a move of center site.
    // We could always keep all edges refreshed but that would be wasteful, since the next iteration
    // may activate other sites and not end up needing those edges.
    // Instead, we force the hand of the algorithm, to only allow edge rebuilds with active sites defined.
    // Ideally, then, this should be done directly after activating new sites in a new iteration.
    if(edges.active_sites.empty())
        throw except::runtime_error("rebuild_edges_ene: no active sites.\n"
                                    "Hint:\n"
                                    " One could in principle keep edges refreshed always, but\n"
                                    " that would imply rebuilding many edges that end up not\n"
                                    " being used. Make sure to only run this rebuild after\n"
                                    " activating sites.");

    long   current_position = state.get_position<long>();
    size_t posL_active      = edges.active_sites.front();
    size_t posR_active      = edges.active_sites.back();

    if constexpr(settings::debug or settings::debug_edges)
        tools::log->trace("rebuild_edges_ene: pos {} | dir {} | "
                          "inspecting edges eneL from [{} to {}]",
                          current_position, state.get_direction(), min_pos, posL_active);
    std::vector<size_t> env_pos_log;
    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &env_here = edges.get_env_eneL(pos);
        auto  id_here  = env_here.get_unique_id();
        if(pos == 0) env_here.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
        if(not env_here.has_block()) throw except::runtime_error("rebuild_edges_ene: No eneL block detected at pos {}", pos);
        if(pos >= std::min(posL_active, state.get_length() - 1)) continue;
        auto &env_rght = edges.get_env_eneL(pos + 1);
        auto  id_rght  = env_rght.get_unique_id();
        env_rght.refresh(env_here, state.get_mps_site(pos), model.get_mpo(pos));
        if(id_here != env_here.get_unique_id()) env_pos_log.emplace_back(env_here.get_position());
        if(id_rght != env_rght.get_unique_id()) env_pos_log.emplace_back(env_rght.get_position());
    }
    if(not env_pos_log.empty()) tools::log->trace("rebuild_edges_ene: rebuilt eneL edges: {}", env_pos_log);

    env_pos_log.clear();
    if constexpr(settings::debug or settings::debug_edges)
        tools::log->trace("rebuild_edges_ene: pos {} | dir {} | "
                          "inspecting edges eneR from [{} to {}]",
                          current_position, state.get_direction(), posR_active, max_pos);
    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); --pos) {
        auto &env_here = edges.get_env_eneR(pos);
        auto  id_here  = env_here.get_unique_id();
        if(pos == state.get_length() - 1) env_here.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
        if(not env_here.has_block()) throw except::runtime_error("rebuild_edges_ene: No eneR block detected at pos {}", pos);
        if(pos <= std::max(posR_active, 0ul)) continue;
        auto &env_left = edges.get_env_eneR(pos - 1);
        auto  id_left  = env_left.get_unique_id();
        env_left.refresh(env_here, state.get_mps_site(pos), model.get_mpo(pos));
        if(id_here != env_here.get_unique_id()) env_pos_log.emplace_back(env_here.get_position());
        if(id_left != env_left.get_unique_id()) env_pos_log.emplace_back(env_left.get_position());
    }
    std::reverse(env_pos_log.begin(), env_pos_log.end());
    if(not env_pos_log.empty()) tools::log->trace("rebuild_edges_ene: rebuilt eneR edges: {}", env_pos_log);
    if(not edges.get_env_eneL(posL_active).has_block()) throw except::logic_error("rebuild_edges_ene: active env eneL has undefined block");
    if(not edges.get_env_eneR(posR_active).has_block()) throw except::logic_error("rebuild_edges_ene: active env eneR has undefined block");
}

void tools::finite::env::rebuild_edges_var(const StateFinite &state, const ModelFinite &model, EdgesFinite &edges) {
    if(state.get_algorithm() == AlgorithmType::fLBIT) throw except::logic_error("rebuild_edges_var: fLBIT algorithm should never rebuild variance edges!");
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw except::runtime_error("rebuild_edges_var: All lengths not equal: state {} | model {} | edges {}", state.get_length(), model.get_length(),
                                    edges.get_length());
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw except::runtime_error("rebuild_edges_var: All active sites are not equal: state {} | model {} | edges {}", state.active_sites, model.active_sites,
                                    edges.active_sites);
    auto t_reb = tid::tic_scope("rebuild_edges_var", tid::level::higher);

    size_t min_pos = 0;
    size_t max_pos = state.get_length() - 1;

    // If there are no active sites we shouldn't be rebuilding edges.
    // For instance, the active sites are cleared after a move of center site.
    // We could always keep all edges refreshed but that would be wasteful, since the next iteration
    // may activate other sites and not end up needing those edges.
    // Instead, we force the hand of the algorithm, to only allow edge rebuilds with active sites defined.
    // Ideally, then, this should be done directly after activating new sites in a new iteration.
    if(edges.active_sites.empty())
        throw except::runtime_error("rebuild_edges_var: no active sites.\n"
                                    "Hint:\n"
                                    " One could in principle keep edges refreshed always, but\n"
                                    " that would imply rebuilding many edges that end up not\n"
                                    " being used. Make sure to only run this assertion after\n"
                                    " activating sites.");

    long   current_position = state.get_position<long>();
    size_t posL_active      = edges.active_sites.front();
    size_t posR_active      = edges.active_sites.back();

    if constexpr(settings::debug or settings::debug_edges) {
        tools::log->trace("rebuild_edges_var: pos {} | dir {} | "
                          "inspecting edges varL from [{} to {}]",
                          current_position, state.get_direction(), min_pos, posL_active);
    }

    std::vector<size_t> env_pos_log;
    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &env_here = edges.get_env_varL(pos);
        auto  id_here  = env_here.get_unique_id();
        if(pos == 0) env_here.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
        if(not env_here.has_block()) throw except::runtime_error("rebuild_edges_var: No varL block detected at pos {}", pos);
        if(pos >= std::min(posL_active, state.get_length() - 1)) continue;
        auto &env_rght = edges.get_env_varL(pos + 1);
        auto  id_rght  = env_rght.get_unique_id();
        env_rght.refresh(env_here, state.get_mps_site(pos), model.get_mpo(pos));
        if(id_here != env_here.get_unique_id()) env_pos_log.emplace_back(env_here.get_position());
        if(id_rght != env_rght.get_unique_id()) env_pos_log.emplace_back(env_rght.get_position());
    }

    if(not env_pos_log.empty()) tools::log->trace("rebuild_edges_var: rebuilt varL edges: {}", env_pos_log);
    env_pos_log.clear();
    if constexpr(settings::debug or settings::debug_edges) {
        tools::log->trace("rebuild_edges_var: pos {} | dir {} | "
                          "inspecting edges varR from [{} to {}]",
                          current_position, state.get_direction(), posR_active, max_pos);
    }

    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); --pos) {
        auto &env_here = edges.get_env_varR(pos);
        auto  id_here  = env_here.get_unique_id();
        if(pos == state.get_length() - 1) env_here.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
        if(not env_here.has_block()) throw except::runtime_error("rebuild_edges_var: No varR block detected at pos {}", pos);
        if(pos <= std::max(posR_active, 0ul)) continue;
        auto &env_left = edges.get_env_varR(pos - 1);
        auto  id_left  = env_left.get_unique_id();
        env_left.refresh(env_here, state.get_mps_site(pos), model.get_mpo(pos));
        if(id_here != env_here.get_unique_id()) env_pos_log.emplace_back(env_here.get_position());
        if(id_left != env_left.get_unique_id()) env_pos_log.emplace_back(env_left.get_position());
    }
    std::reverse(env_pos_log.begin(), env_pos_log.end());
    if(not env_pos_log.empty()) tools::log->trace("rebuild_edges_var: rebuilt varR edges: {}", env_pos_log);
    if(not edges.get_env_varL(posL_active).has_block()) throw except::logic_error("rebuild_edges_var: active env varL has undefined block");
    if(not edges.get_env_varR(posR_active).has_block()) throw except::logic_error("rebuild_edges_var: active env varR has undefined block");
}

void tools::finite::env::rebuild_edges(const StateFinite &state, const ModelFinite &model, EdgesFinite &edges) {
    if(state.get_algorithm() == AlgorithmType::fLBIT) return;
    rebuild_edges_ene(state, model, edges);
    rebuild_edges_var(state, model, edges);
}
