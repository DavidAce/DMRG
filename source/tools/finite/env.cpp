//
// Created by david on 2020-05-13.
//
#include "env.h"
#include <math/num.h>
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
#include <general/nmspc_tensor_extra.h>

void tools::finite::env::expand_subspace(class_state_finite &state, const class_model_finite &model, class_edges_finite &edges, double alpha, long chi_lim, std::optional<double> svd_threshold){
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw std::runtime_error(
            fmt::format("All lengths not equal: state {} | model {} | edges {}", state.get_length(), model.get_length(), edges.get_length()));
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw std::runtime_error(
            fmt::format("All active sites are not equal: state {} | model {} | edges {}", state.active_sites, model.active_sites, edges.active_sites));
    if(alpha < 1e-12) return;
    if(state.active_sites.empty()) throw std::runtime_error("No active sites for subspace expansion");
    using Scalar = class_state_finite::Scalar;
    // Follows the subspace expansion technique explained in https://link.aps.org/doi/10.1103/PhysRevB.91.155115
    {
        auto posL = state.active_sites.front();
        // Construct PL = alpha * eneL(i-1) * mps(i-1) * mpo(i-1)
        if(posL > 0 and posL < state.get_length<size_t>() and state.get_direction()==1){
            auto & mpsL = state.get_mps_site(posL-1); // The site to the left
            auto & mpsR = state.get_mps_site(posL); // The site to the right
            if(mpsL.get_label() == "A" or mpsL.get_label() == "AC"){
                auto chiL_old = mpsL.get_chiR();
                auto dimL_old = mpsL.dimensions();
                auto dimR_old = mpsR.dimensions();
                auto &mpoL = model.get_mpo(posL-1);
                auto &eneL = edges.get_eneL(posL-1);
                auto &varL = edges.get_varL(posL-1);
                auto PL1 = eneL.get_expansion_term(mpsL,mpoL,alpha);
                auto PL2 = varL.get_expansion_term(mpsL,mpoL,alpha);
                Eigen::Tensor<Scalar,3> PL = PL1.concatenate(PL2,2);
                auto zero = Textra::TensorConstant<Scalar>(0.0, mpsR.spin_dim(), PL.dimension(2), mpsR.get_chiR());
                auto ones = Textra::TensorConstant<Scalar>(1.0, PL.dimension(2));

                Eigen::Tensor<Scalar,3> ML_PL = Textra::asNormalized(mpsL.get_M_bare().concatenate(PL,2));
                Eigen::Tensor<Scalar,3> MR_zero = mpsR.get_M_bare().concatenate(zero, 1);

                if(mpsL.isCenter()){
                    Eigen::Tensor<Scalar,1> LC_ones = mpsL.get_LC().concatenate(ones, 0);
                    mpsL.set_LC(LC_ones,mpsL.get_truncation_error_LC());
                    tools::log->info("Subspace expansion -> mpsL({}) label {}: Set LC {} with norm {:.16f}",
                                     mpsL.get_position(),mpsL.get_label(), LC_ones.dimensions(), Textra::norm(LC_ones));
                }
                if(mpsL.get_label() == "B"){
                    Eigen::Tensor<Scalar,1> L_ones = mpsL.get_L().concatenate(ones, 0);
                    mpsL.set_L(L_ones,mpsL.get_truncation_error());
                    tools::log->info("Subspace expansion -> mpsL({}) label {}: Set L {} with norm {:.16f}",
                                     mpsL.get_position(),mpsL.get_label(), L_ones.dimensions(), Textra::norm(L_ones));
                }
                if(mpsR.get_label().find('A') != std::string::npos){
                    Eigen::Tensor<Scalar,1> L_ones = mpsR.get_L().concatenate(ones, 0);
                    mpsR.set_L(L_ones,mpsR.get_truncation_error());
                    tools::log->info("Subspace expansion -> mpsR({}) label {}: Set L {} with norm {:.16f}",
                                     mpsR.get_position(),mpsR.get_label(), L_ones.dimensions(), Textra::norm(L_ones));
                }

                chi_lim = std::min(chi_lim, mpsL.spin_dim() * std::min(mpsL.get_chiL(), mpsR.get_chiR())); // Bond dimension can't grow faster than x spin_dim
                mpsL.set_M(ML_PL);
                mpsR.set_M(MR_zero);
                state.clear_cache();

                // Make both left-normalized
                long center = state.get_position<long>();
                auto multisite_tensor = state.get_multisite_mps({posL-1});
                tools::finite::mps::merge_multisite_tensor(state,multisite_tensor,{posL-1}, center, chi_lim, svd_threshold);

                auto mpsL_new = state.get_mps_site(posL-1); // The site to the left
                auto mpsR_new = state.get_mps_site(posL); // The site to the right
                auto chiL_new = mpsL_new.get_chiR();
                auto dimL_new = mpsL_new.dimensions();
                auto dimR_new = mpsR_new.dimensions();
                tools::log->info("Subspace expansion pos L {} {} | alpha {:.2e} | ML {} -> {} | MR {} -> {} | ML_PL {} | MR_0 {} | χ {} -> {} -> {} | χlim {} ",
                                 posL-1, posL,alpha ,dimL_old, dimL_new, dimR_old, dimR_new, ML_PL.dimensions(),MR_zero.dimensions(), chiL_old, multisite_tensor.dimension(2),chiL_new, chi_lim);
                if(dimL_old[1] != dimL_new[1] ) throw std::runtime_error(fmt::format("mpsL changed chiL during left-moving expansion: {} -> {}",dimL_old, dimL_new  ));
                if(dimR_old[2] != dimR_new[2] ) throw std::runtime_error(fmt::format("mpsR changed chiR during left-moving expansion: {} -> {}",dimR_old, dimR_new  ));
            }
        }
    }

    {
        auto posR = state.active_sites.back();
        // Construct PL = alpha * eneL(i+1) * mps(i+1) * mpo(i+1)
        if(posR < state.get_length<size_t>()-1 and state.get_direction()==-1){
            auto & mpsR = state.get_mps_site(posR+1); // The the site to the right
            auto & mpsL = state.get_mps_site(posR); // The site to the left
            if(mpsR.get_label() == "B"){
                auto chiR_old = mpsR.get_chiL();
                auto dimR_old = mpsR.dimensions();
                auto dimL_old = mpsL.dimensions();
                auto &mpoR = model.get_mpo(posR+1);
                auto &eneR = edges.get_eneR(posR+1);
                auto &varR = edges.get_varR(posR+1);
                auto PR1 = eneR.get_expansion_term(mpsR,mpoR,alpha);
                auto PR2 = varR.get_expansion_term(mpsR,mpoR,alpha);
                Eigen::Tensor<Scalar,3> PR = PR1.concatenate(PR2,1);
                auto zero = Textra::TensorConstant<Scalar>(0.0, mpsL.spin_dim(), mpsL.get_chiL(), PR.dimension(1));
                auto ones = Textra::TensorConstant<Scalar>(1.0, PR.dimension(1));

                Eigen::Tensor<Scalar,3> MR_PR = mpsR.get_M_bare().concatenate(PR,1);
                Eigen::Tensor<Scalar,3> ML_zero = mpsL.get_M_bare().concatenate(zero, 2);

                if(mpsL.isCenter()){
                    Eigen::Tensor<Scalar,1> LC_ones = mpsL.get_LC().concatenate(ones, 0);
                    mpsL.set_LC(LC_ones,mpsL.get_truncation_error_LC());
                    tools::log->info("Subspace expansion <- mpsL({}) label {}: Set LC {} with norm {:.16f}",
                                     mpsL.get_position(),mpsL.get_label(), LC_ones.dimensions(), Textra::norm(LC_ones));
                }
                if(mpsL.get_label() == "B"){
                    Eigen::Tensor<Scalar,1> L_ones = mpsL.get_L().concatenate(ones, 0);
                    mpsL.set_L(L_ones,mpsL.get_truncation_error());
                    tools::log->info("Subspace expansion <- mpsL({}) label {}: Set L {} with norm {:.16f}",
                                     mpsL.get_position(),mpsL.get_label(), mpsL.get_L().dimensions(), Textra::norm(L_ones));
                }
                if(mpsR.get_label().find('A') != std::string::npos){
                    Eigen::Tensor<Scalar,1> L_ones = mpsR.get_L().concatenate(ones, 0);
                    mpsR.set_L(L_ones,mpsR.get_truncation_error());
                    tools::log->info("Subspace expansion <- mpsR({}) label {}: Set L {} with norm {:.16f}",
                                     mpsR.get_position(),mpsR.get_label(), mpsR.get_L().dimensions(), Textra::norm(L_ones));
                }
                chi_lim = std::min(chi_lim, mpsR.spin_dim() * std::min(mpsL.get_chiL(), mpsR.get_chiR())); // Bond dimension can't grow faster than x spin_dim
                mpsL.set_M(ML_zero);
                mpsR.set_M(Textra::asNormalized(MR_PR));
                state.clear_cache();


                // Make right-normalized
                Eigen::Tensor<Scalar,3> multisite_tensor;
                if(mpsL.isCenter()) multisite_tensor = Textra::asNormalized(
                                                            mpsR.get_M()
                                                           .contract(Textra::asDiagonal(mpsL.get_LC()), Textra::idx({1},{0}))
                                                           .shuffle(Textra::array3{0,2,1}));

                long center = state.get_position<long>();
//                auto multisite_tensor = state.get_multisite_mps({posR+1});
                tools::finite::mps::merge_multisite_tensor(state,multisite_tensor,{posR+1}, center, chi_lim, svd_threshold);
                auto mpsR_new = state.get_mps_site(posR+1);
                auto mpsL_new = state.get_mps_site(posR);
                auto chiR_new = mpsR_new.get_chiL();
                auto dimR_new = mpsR_new.dimensions();
                auto dimL_new = mpsL_new.dimensions();

                tools::log->info("Subspace expansion pos R {} {} | alpha {:.2e} | ML {} -> {} | MR {} -> {} | ML_0 {} | MR_PR {} | PR {} | 0 {} | χ {} -> {} -> {} | χlim {} ",
                                 posR, posR+1,alpha, dimL_old, dimL_new, dimR_old, dimR_new, ML_zero.dimensions(), MR_PR.dimensions(), PR.dimensions(), zero.dimensions(), chiR_old, multisite_tensor.dimension(1),chiR_new, chi_lim);
                if(dimL_old[1] != dimL_new[1] ) throw std::runtime_error(fmt::format("mpsL changed chiL during right-moving expansion: {} -> {}",dimL_old, dimL_new  ));
                if(dimR_old[2] != dimR_new[2] ) throw std::runtime_error(fmt::format("mpsR changed chiR during right-moving expansion: {} -> {}",dimR_old, dimR_new  ));
            }
        }
    }


    env::rebuild_edges(state,model,edges);
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
    tools::common::profile::get_default_prof()["t_env"]->tic();

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
    tools::common::profile::get_default_prof()["t_env"]->toc();
}

void tools::finite::env::rebuild_edges_var(const class_state_finite &state, const class_model_finite &model, class_edges_finite &edges) {
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw std::runtime_error(
            fmt::format("All lengths not equal: state {} | model {} | edges {}", state.get_length(), model.get_length(), edges.get_length()));
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw std::runtime_error(
            fmt::format("All active sites are not equal: state {} | model {} | edges {}", state.active_sites, model.active_sites, edges.active_sites));
    tools::common::profile::get_default_prof()["t_env"]->tic();

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
    tools::common::profile::get_default_prof()["t_env"]->toc();
}

void tools::finite::env::rebuild_edges(const class_state_finite &state, const class_model_finite &model, class_edges_finite &edges) {
    rebuild_edges_ene(state, model, edges);
    rebuild_edges_var(state, model, edges);
}
