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

/*!
    Follows the environment (aka subspace) expansion technique as explained in https://link.aps.org/doi/10.1103/PhysRevB.91.155115
    Note that in this convention we expand/enrich the trailing bond during a DMRG sweep.
    Going left-to-right with active site == i, we update:

        - AC(i)   --> SVD( [AC(i)  , PL(i)] ) = U*S*V, to update AC(i) = U, LC = S and V is contracted onto B(i+1) below
        - B(i+1)  --> V * [B(i+1) , 0   ]^T  : This site loses right-normalization when we contract V from the SVD above.

    where PL(i) = alpha * ENVL(i) * AC(i) * MPO(i) (dimensions d, AC.dimension(2), AC.dimension(0)*MPO(i).dimension(1))
    Then, immediately after calling this function, we should move the current center site i --> i+1

    Similarly, going right-to-left with active site == i, we update:

        - AC(i)   --> SVD( [AC(i) , PR(i) ]^T ) = U*S*V, to update AC(i)=S*V, and U is contracted onto A(i-1) below
        - A(i-1)  -->      [A(i-1), 0    ] * U  : This one loses left-normalization when we contract U from the SVD above.

    where PR(i) = alpha * ENVR(i) * AC(i) * MPO(i) (dimensions d, AC.dimension(1), AC.dimension(2) * MPO(i).dimension(0))
    Then, immediately after calling this function, we should move the current center site i --> i-1
*/
std::vector<size_t> tools::finite::env::expand_environment_backward(StateFinite &state, const ModelFinite &model, EdgesFinite &edges,
                                                                    std::optional<double> alpha, EnvExpandMode envExpandMode,
                                                                    std::optional<svd::config> svd_cfg) {
    auto                pos = state.get_position<size_t>();
    std::vector<size_t> pos_expanded;
    if(state.get_direction() > 0 and pos == std::clamp<size_t>(pos, 0, state.get_length<size_t>() - 2)) pos_expanded = {pos, pos + 1};
    if(state.get_direction() < 0 and pos == std::clamp<size_t>(pos, 1, state.get_length<size_t>() - 1)) pos_expanded = {pos - 1, pos};
    if(not alpha or pos_expanded.empty()) {
        // When alpha == std::nullopt this turns into a position query
        // Just return the positions that would get modified if alpha had been defined.
        return pos_expanded;
    }

    // Define the left and right mps that will get modified
    state.clear_cache();
    auto  posL     = pos_expanded.front();
    auto  posR     = pos_expanded.back();
    auto &mpsL     = state.get_mps_site(posL);
    auto &mpsR     = state.get_mps_site(posR);
    auto  dimL_old = mpsL.dimensions();
    auto  dimR_old = mpsR.dimensions();

    // Set up the SVD
    auto bond_lim = std::min(mpsL.spin_dim() * mpsL.get_chiL(), mpsR.spin_dim() * mpsR.get_chiR()); // Bond dimension can't grow faster than x spin_dim
    svd_cfg       = svd_cfg.value_or(svd::config());
    svd_cfg->truncation_limit = svd_cfg->truncation_limit.value_or(settings::precision::svd_truncation_lim);
    svd_cfg->rank_max         = std::min(bond_lim, svd_cfg->rank_max.value_or(bond_lim));
    svd::solver svd(svd_cfg);

    if(state.get_direction() > 0) {
        // The expanded bond sits between mpsL and mpsR. When direction is left-to-right:
        //      * mpsL is "AC(i)" and is the active site that will get moved into envL later
        //      * mpsR is "B(i+1)" and will become the active site after moving center position
        auto                  &mpoL = model.get_mpo(posL);
        Eigen::Tensor<cplx, 3> PL;
        switch(envExpandMode) {
            case EnvExpandMode::ENE: PL = edges.get_env_eneL(posL).get_expansion_term(mpsL, mpoL, alpha.value()); break;
            case EnvExpandMode::VAR: PL = edges.get_env_varL(posL).get_expansion_term(mpsL, mpoL, alpha.value()); break;
        }
        Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsR.spin_dim(), PL.dimension(2), mpsR.get_chiR());
        Eigen::Tensor<cplx, 3> ML_PL = mpsL.get_M().concatenate(PL, 2); // mpsL is going to be optimized, enriched with PL
        Eigen::Tensor<cplx, 3> MR_P0 = mpsR.get_M().concatenate(P0, 1); // mpsR is going into the environment, padded with zeros

        auto [U, S, V] = svd.schmidt_into_left_normalized(ML_PL, mpsL.spin_dim(), svd_cfg.value());
        mpsL.set_M(U);
        mpsL.set_LC(S, -1.0); // Set a negative truncation error to ignore it.
        mpsL.stash_V(V, posR);
        mpsR.set_M(MR_P0);
        mpsR.take_stash(mpsL); // right normalization of mpsR is lost here -- but we will move right soon, thereby left-normalizing it into an AC
        tools::log->debug("Environment expansion backward pos {} | {} | alpha {:.2e} | svd_ε {:.2e} | χlim {} | χ {} -> {} -> {}", pos_expanded,
                          enum2sv(envExpandMode), alpha.value(), svd_cfg->truncation_limit.value_or(std::numeric_limits<double>::quiet_NaN()), bond_lim,
                          dimL_old[2], ML_PL.dimension(2), mpsL.get_chiR());
    }
    if(state.get_direction() < 0) {
        // The expanded bond sits between mpsL and mpsR. When direction is left-to-right:
        //      * mpsL is "A(i-1)" and will become the active site after moving center position
        //      * mpsR is "AC(i)" and is the active site that will get moved into envR later
        auto                  &mpoR = model.get_mpo(posR);
        Eigen::Tensor<cplx, 3> PR;
        switch(envExpandMode) {
            case EnvExpandMode::ENE: PR = edges.get_env_eneR(posR).get_expansion_term(mpsR, mpoR, alpha.value()); break;
            case EnvExpandMode::VAR: PR = edges.get_env_varR(posR).get_expansion_term(mpsR, mpoR, alpha.value()); break;
        }
        Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsL.spin_dim(), mpsL.get_chiL(), PR.dimension(1));
        Eigen::Tensor<cplx, 3> ML_P0 = mpsL.get_M_bare().concatenate(P0, 2); // [ A   0 ]
        Eigen::Tensor<cplx, 3> MR_PR = mpsR.get_M_bare().concatenate(PR, 1); // [ AC  P ]^T including LC
        auto [U, S, V]               = svd.schmidt_into_right_normalized(MR_PR, mpsR.spin_dim(), svd_cfg.value());
        // We have right-normalized the a bare AC = L * A (not including * LC)
        // Therefore the resulting the new truncated AC should include the singular values as well, that is: AC = L * A --> S * V
        auto SV = Eigen::Tensor<cplx, 3>(tenx::asDiagonal(S).contract(V, tenx::idx({1}, {1})).shuffle(tenx::array3{1, 0, 2}));
        mpsR.set_M(SV); // This AC*LC becomes right-normalized! Note that L*A*LC is still a valid normalized "multisite_mps"
        mpsR.set_L(S);  // We keep track of L in an A type tensor, even though it is included in A already
        mpsR.stash_U(U, posL);
        mpsL.set_M(ML_P0);
        mpsL.take_stash(mpsR);
        tools::log->debug("Environment expansion backward pos {} | {} | alpha {:.2e} | svd_ε {:.2e} | χlim {} | χ {} -> {} -> {}", pos_expanded,
                          enum2sv(envExpandMode), alpha.value(), svd_cfg->truncation_limit.value_or(std::numeric_limits<double>::quiet_NaN()), bond_lim,
                          dimL_old[2], MR_PR.dimension(1), mpsL.get_chiR());
    }

    if(dimL_old[1] != mpsL.get_chiL()) throw except::runtime_error("mpsL changed chiL during environment expansion: {} -> {}", dimL_old, mpsL.dimensions());
    if(dimR_old[2] != mpsR.get_chiR()) throw except::runtime_error("mpsR changed chiR during environment expansion: {} -> {}", dimR_old, mpsR.dimensions());
    // if constexpr(settings::debug or settings::debug_expansion) mpsL.assert_normalized();
    // if constexpr(settings::debug or settings::debug_expansion) mpsR.assert_normalized();
    env::rebuild_edges(state, model, edges);
    return pos_expanded;
}

/*!
    This is a modified (aka subspace) expansion technique, compared to the one explained in https://link.aps.org/doi/10.1103/PhysRevB.91.155115
    In this convention we expand/enrich the forward bond during a DMRG sweep, i.e. the bond one step ahead in the sweep direction.

    Going left-to-right with active site == i, we update:

        - AC(i), LC(i) -->      [AC(i)     ,  0 ]  * U, S(i)      : AC(i) is bare, without LC(i) (which is used below). Loses left-normalization due to U.
        - LC*B(i+1)    --> SVD( [LC*B(i+1) , PR(i+1) ]^T ) = U*S*V: update B(i+1) = V, LC(i) = S and U is contracted onto AC(i) above

    where PR(i+1) = alpha * ENVR(i+1) * B(i+1) * MPO(i+1) (dimensions d(i-1), B(i+1).dimension(1), B(i+1).dimension(2) * MPO(i+1).dimension(0))
    Then, immediately after calling this function, we should optimize the current site AC(i) with a DMRG step.

    Similarly, going right-to-left with active site == i, we update:

         - A(i-1)*L(i)  --> SVD( [A(i-1)*L(i), PL(i-1)] ) = U*S*V : update A(i-1) = U, L(i) = S and S*V is contracted onto AC(i) below.
         - L(i), AC(i)  -->  V*[AC(i) ,  0      ]^T               : AC(i) is bare, without LC(i). Loses left-normalization due to V.

    where PL(i-1) = alpha * ENVL(i-1) * A(i-1) * MPO(i-1) (dimensions d(i-1), A(i-1).dimension(2), A(i-1).dimension(0)*MPO(i-1).dimension(1))
    Then, immediately after calling this function, we should optimize the current site AC(i) with a DMRG step.


    Thus, one step of the DMRG algorithm proceeds as

    1. expand between sites i, i+1
    2. optimize $\Psi^{i} = A_C^{i} \Lambda_C^{i}$.
    3. split $\Psi^{i} \stackrel{\text{SVD}}{\rightarrow}A_C^{i}  \Lambda_C^{i} V^\dagger B^{i+1}$ normally and update the MPS on sites $i,i+1$.
    4. move: $A_C^i \Lambda_C^i \rightarrow A^i$ and $\Lambda_C^i B^{i+1} \rightarrow \Lambda^i A_C^{i+1} \Lambda_C^{i+1}$, and repeat from step 1.
    5. update $\alpha$: decrease by $0.01/L$ if **if the variance decreased** in the last step, else increase by $10/L$, bounded by $\alpha \in [\alpha_{\min}, \alpha_{\max}]$, where:

      - $\alpha_{\min} = \text{Var}(H) \cdot 10^{-3}$,
      - $\alpha_{\max}= \min{(10^{-4}, \text{Var}(H))}$.

    The main differences compared with the original:
    - The SVD in the expansion can be done with high precision and bond dimension: the real truncation happens after the optimization in step 3.
    - Therefore, the optimization step sees a much richer environment, which speeds up convergence.
    - Backwards expansion pushes a non-optimal/degraded site-tensor into the trailing environment. Forward expansion optimizes it first.

*/

std::vector<size_t> tools::finite::env::expand_environment_forward(StateFinite &state, const ModelFinite &model, EdgesFinite &edges,
                                                                   std::optional<double> alpha, EnvExpandMode envExpandMode,
                                                                   std::optional<svd::config> svd_cfg) {
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw except::runtime_error("expand_environment_forward: All lengths not equal: state {} | model {} | edges {}", state.get_length(), model.get_length(),
                                    edges.get_length());
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw except::runtime_error("expand_environment_forward: All active sites are not equal: state {} | model {} | edges {}", state.active_sites,
                                    model.active_sites, edges.active_sites);
    if(state.active_sites.empty()) throw except::logic_error("No active sites for environment expansion");
    if(state.active_sites.size() != 1) {
        tools::log->warn("expand_environment_forward: expected active_sites.size() == 1. Got: {}", state.active_sites);
        return {};
    }
    if(state.active_sites.front() != state.get_position()) {
        tools::log->warn("expand_environment_forward: expected active_sites.front() == current position. Got: {} != {}", state.active_sites,
                         state.get_position());
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

    // Define the left and right mps that will get modified
    state.clear_cache();
    auto  posL     = pos_expanded.front();
    auto  posR     = pos_expanded.back();
    auto &mpsL     = state.get_mps_site(posL);
    auto &mpsR     = state.get_mps_site(posR);
    auto  dimL_old = mpsL.dimensions();
    auto  dimR_old = mpsR.dimensions();

    // Set up the SVD
    // Bond dimension can't grow faster than x spin_dim, but we can generate a highly enriched environment here for optimization,
    // and let the proper truncation happen after optimization instead.
    auto bond_lim             = std::min(mpsL.spin_dim() * mpsL.get_chiL(), mpsR.spin_dim() * mpsR.get_chiR());
    svd_cfg                   = svd_cfg.value_or(svd::config());
    svd_cfg->truncation_limit = svd_cfg->truncation_limit.value_or(settings::precision::svd_truncation_lim);
    svd_cfg->rank_max         = std::min(bond_lim, svd_cfg->rank_max.value_or(bond_lim));
    svd::solver svd(svd_cfg);

    if(state.get_direction() > 0) {
        // The expanded bond sits between mpsL and mpsR. When direction is left-to-right:
        //      * mpsL is "AC(i)"  and is the active site
        //      * mpsR is "B(i+1)" and belongs in envR later in the optimization step
        auto &mpoR = model.get_mpo(posR);
        mpsR.set_M(state.get_multisite_mps({posR})); // mpsR absorbs the LC(i) bond from the left so that the SVD makes sense later
        Eigen::Tensor<cplx, 3> PR;
        switch(envExpandMode) {
            case EnvExpandMode::ENE: PR = edges.get_env_eneR(posR).get_expansion_term(mpsR, mpoR, alpha.value()); break;
            case EnvExpandMode::VAR: PR = edges.get_env_varR(posR).get_expansion_term(mpsR, mpoR, alpha.value()); break;
        }
        Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsL.spin_dim(), mpsL.get_chiL(), PR.dimension(1));
        Eigen::Tensor<cplx, 3> ML_P0 = mpsL.get_M_bare().concatenate(P0, 2); // mpsL is going to be optimized, zero-padded with P0
        Eigen::Tensor<cplx, 3> MR_PR = mpsR.get_M_bare().concatenate(PR, 1); // mpsR is going into the environment, enriched with PR.

        auto [U, S, V] = svd.schmidt_into_right_normalized(MR_PR, mpsR.spin_dim(), svd_cfg.value());
        mpsR.set_M(V);
        mpsR.stash_U(U, posL);
        mpsR.stash_C(S, -1.0, posL); // Set a negative truncation error to ignore it.
        mpsL.set_M(ML_P0);
        mpsL.take_stash(mpsR); // normalization of mpsL is lost here.

        {
            // Make mpsL normalized. This is strictly optional, but we do it so that normalization checks can succeed
            cplx                   norm_old = tools::common::contraction::contract_mps_norm(mpsL.get_M());
            Eigen::Tensor<cplx, 3> M_tmp    = mpsL.get_M_bare() * mpsL.get_M_bare().constant(std::pow(norm_old, -0.5)); // Rescale
            mpsL.set_M(M_tmp);
            if constexpr(settings::debug or settings::debug_expansion) {
                auto mpsL_final = state.get_multisite_mps({mpsL.get_position()});
                cplx norm_new   = tools::common::contraction::contract_mps_norm(mpsL_final);
                tools::log->debug("Normalized expanded mps {}({}): {:.16f} -> {:.16f}", mpsL.get_label(), mpsL.get_position(), std::abs(norm_old),
                                  std::abs(norm_new));
            }
        }

        tools::log->debug("Environment expansion forward pos {} | alpha {:.2e} | χ {} -> {} -> {} | svd_ε {:.2e} | χlim {}", pos_expanded, alpha.value(),
                          dimL_old[2], ML_P0.dimension(2), mpsL.get_chiR(), svd_cfg->truncation_limit, svd_cfg->rank_max);
    }
    if(state.get_direction() < 0) {
        // The expanded bond sits between mpsL and mpsR. When direction is right-to-left:
        //      * mpsL is "A(i-1)" and belongs in envL later in the optimization step
        //      * mpsR is "AC(i)" and is the active site
        auto &mpoL = model.get_mpo(posL);
        mpsL.set_M(state.get_multisite_mps({posL})); // mpsL absorbs the L(i) bond from the right so that the SVD makes sense later
        Eigen::Tensor<cplx, 3> PL;
        switch(envExpandMode) {
            case EnvExpandMode::ENE: PL = edges.get_env_eneL(posL).get_expansion_term(mpsL, mpoL, alpha.value()); break;
            case EnvExpandMode::VAR: PL = edges.get_env_varL(posL).get_expansion_term(mpsL, mpoL, alpha.value()); break;
        }
        Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsR.spin_dim(), PL.dimension(2), mpsR.get_chiR());
        Eigen::Tensor<cplx, 3> ML_PL = mpsL.get_M_bare().concatenate(PL, 2);
        Eigen::Tensor<cplx, 3> MR_P0 = mpsR.get_M_bare().concatenate(P0, 1);

        auto [U, S, V] = svd.schmidt_into_left_normalized(ML_PL, mpsL.spin_dim(), svd_cfg.value());
        // Recall: the following lines hold because we demand that mpsL is an "A" and mpsR an "AC"
        mpsL.set_M(U);
        mpsL.stash_S(S, -1.0, posR); // Set a negative truncation error to ignore it.
        mpsL.stash_V(V, posR);
        mpsR.set_M(MR_P0);
        mpsR.take_stash(mpsL); // normalization of mpsR is lost here

        {
            // Make mpsR normalized. This is strictly optional, but we do it so that normalization checks can succeed
            cplx                   norm_old = tools::common::contraction::contract_mps_norm(mpsR.get_M());
            Eigen::Tensor<cplx, 3> M_tmp    = mpsR.get_M_bare() * mpsR.get_M_bare().constant(std::pow(norm_old, -0.5)); // Rescale by the norm
            mpsR.set_M(M_tmp);
            if constexpr(settings::debug or settings::debug_expansion) {
                auto mpsR_final = state.get_multisite_mps({mpsR.get_position()});
                cplx norm_new   = tools::common::contraction::contract_mps_norm(mpsR_final);
                tools::log->debug("Normalized expanded mps {}({}): {:.16f} -> {:.16f}", mpsR.get_label(), mpsR.get_position(), std::abs(norm_old),
                                  std::abs(norm_new));
            }
        }
        tools::log->debug("Environment expansion forward pos {} | alpha {:.2e} | χ {} -> {} -> {} | ε {:.2e} | χlim {}", pos_expanded, alpha.value(),
                          dimR_old[1], MR_P0.dimension(1), mpsR.get_chiL(), svd_cfg->truncation_limit, svd_cfg->rank_max);
    }

    if(dimL_old[1] != mpsL.get_chiL()) throw except::runtime_error("mpsL changed chiL during environment expansion: {} -> {}", dimL_old, mpsL.dimensions());
    if(dimR_old[2] != mpsR.get_chiR()) throw except::runtime_error("mpsR changed chiR during environment expansion: {} -> {}", dimR_old, mpsR.dimensions());
    if constexpr(settings::debug or settings::debug_expansion) mpsL.assert_normalized();
    if constexpr(settings::debug or settings::debug_expansion) mpsR.assert_normalized();
    state.clear_cache();
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
    //    size_t posL_active      = safe_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
    //    size_t posR_active      = safe_cast<size_t>(std::clamp<long>(current_position, 0, state.get_length<long>() - 1));
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
