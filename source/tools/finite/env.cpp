//
// Created by david on 2020-05-13.
//
#include "env.h"
#include <math/num.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/edges/class_env_ene.h>
#include <tensors/edges/class_env_var.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>


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

//    edges.eject_edges_stale_ene(state.get_edge_ene_status()); // Discard edges past the updated sites

    // If there are no active sites then we can build up until current position
    size_t posL_active = state.get_position();
    size_t posR_active = state.get_position();
    if(not edges.active_sites.empty()) {
        posL_active = edges.active_sites.front();
        posR_active = edges.active_sites.back();
//        edges.eject_edges_inactive_ene(); // Discard edges past the active ones
    }
    // Sometimes environments have been ejected outside of the active sites. These have to
    // be rebuilt so
//    auto [posL_ejected, posR_ejected] = edges.get_ejected_positions_ene();
//    tools::log->debug("Ejected edges: [{},{}]", posL_ejected,posR_ejected);
//    posL_active = std::min(posL_ejected,posL_active);
//    posR_active = std::min(posR_ejected,posR_active);


    tools::log->debug("Inspecting edges eneL from [{} to {}]", min_pos,posL_active);
    std::vector<size_t> ene_pos_log;
    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &ene_curr = edges.get_eneL(pos);
        if(pos == 0 and not ene_curr.has_block()){
            ene_pos_log.emplace_back(pos);
            ene_curr.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
//            state.set_edge_ene_status(pos,EdgeStatus::FRESH);
        }
        if(pos >= std::min(posL_active,state.get_length()-1)) continue;
        auto &ene_next = edges.get_eneL(pos + 1);
        auto id = ene_next.get_unique_id();
        ene_next.refresh(ene_curr, state.get_mps_site(pos), model.get_mpo(pos));
        if(id != ene_next.get_unique_id()) ene_pos_log.emplace_back(ene_next.get_position());
    }
    if(not ene_pos_log.empty()) tools::log->trace("Rebuilt L ene edges: {}", ene_pos_log);
    ene_pos_log.clear();
    tools::log->debug("Inspecting edges eneR from [{} to {}]", posR_active,max_pos);
    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); pos--) {
        auto &ene_curr = edges.get_eneR(pos);
        if(pos == state.get_length() -1 and not ene_curr.has_block()){
            ene_pos_log.emplace_back(pos);
            ene_curr.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
//            state.set_edge_ene_status(pos,EdgeStatus::FRESH);
        }
        if(pos <= std::max(posR_active,0ul)) continue;
        auto &ene_prev = edges.get_eneR(pos - 1);
        auto id = ene_prev.get_unique_id();
        ene_prev.refresh(ene_curr, state.get_mps_site(pos), model.get_mpo(pos));
        if(id != ene_prev.get_unique_id()) ene_pos_log.emplace_back(ene_prev.get_position());
//        if(not ene_prev.has_block()) {
//            ene_pos_log.emplace_back(pos - 1);
//            ene_prev = ene_curr.enlarge(state.get_mps_site(pos), model.get_mpo(pos));
//            state.set_edge_ene_status(pos-1,EdgeStatus::FRESH);
//        }
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
    size_t posL_active = state.get_position();
    size_t posR_active = state.get_position();
    if(not edges.active_sites.empty()) {
        posL_active = edges.active_sites.front();
        posR_active = edges.active_sites.back();
    }
    tools::log->debug("Inspecting edges varL from [{} to {}]", min_pos,posL_active);
    std::vector<size_t> var_pos_log;
    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &var_curr = edges.get_varL(pos);
        if(pos == 0 and not var_curr.has_block()){
            var_pos_log.emplace_back(pos);
            var_curr.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
        }
        if(pos >= std::min(posL_active,state.get_length()-1)) continue;
        auto &var_next = edges.get_varL(pos + 1);
        auto id = var_next.get_unique_id();
        var_next.refresh(var_curr, state.get_mps_site(pos), model.get_mpo(pos));
        if(id != var_next.get_unique_id()) var_pos_log.emplace_back(var_next.get_position());
    }

    if(not var_pos_log.empty()) tools::log->trace("Rebuilt L var edges: {}", var_pos_log);
    var_pos_log.clear();
    tools::log->debug("Inspecting edges varR from [{} to {}]", posR_active,max_pos);
    for(size_t pos = max_pos; pos > posR_active and pos < state.get_length(); pos--) {
        auto &var_curr = edges.get_varR(pos);
        if(pos == state.get_length() - 1 and not var_curr.has_block()){
            var_pos_log.emplace_back(pos);
            var_curr.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
        }
        if(pos <= std::max(posR_active,0ul)) continue;
        auto &var_prev = edges.get_varR(pos - 1);
        auto id = var_prev.get_unique_id();
        var_prev.refresh(var_curr,state.get_mps_site(pos), model.get_mpo(pos));
        if(id != var_prev.get_unique_id()) var_pos_log.emplace_back(var_prev.get_position());
    }
    std::reverse(var_pos_log.begin(), var_pos_log.end());
    if(not var_pos_log.empty()) tools::log->trace("Rebuilt R var edges: {}", var_pos_log);
    if(not edges.get_varL(posL_active).has_block()) throw std::logic_error(fmt::format("Left active var edge has undefined block"));
    if(not edges.get_varR(posR_active).has_block()) throw std::logic_error(fmt::format("Right active var edge has undefined block"));
    tools::common::profile::get_default_prof()["t_env"]->toc();
}



void tools::finite::env::rebuild_edges(const class_state_finite &state, const class_model_finite &model, class_edges_finite &edges) {
    rebuild_edges_ene(state,model,edges);
    rebuild_edges_var(state,model,edges);
}
