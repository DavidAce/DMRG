//
// Created by david on 2020-05-13.
//
#include "env.h"
#include <tensors/edges/class_edges_finite.h>
#include <tensors/edges/class_env_ene.h>
#include <tensors/edges/class_env_var.h>
#include <math/nmspc_math.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>

void tools::finite::env::rebuild_all_edges(const class_state_finite & state, const class_model_finite & model, class_edges_finite & edges) {
    tools::log->trace("Rebuilding all edges");
    if(not math::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw std::runtime_error(fmt::format("All lengths not equal: state {} | model {} | edges {}",state.get_length(), model.get_length(), edges.get_length() ));

    // Note that the last environment should always be empty!
    // Going left to right, that will be the back of the list.
    // Going right to left, that will be the front of the list.
    edges.eject_all_edges();
    size_t min_pos = 0;
    size_t max_pos = state.get_length()-1;
    {
        auto & ene_front = edges.get_eneL(min_pos);
        auto & var_front = edges.get_varL(min_pos);
        ene_front.set_edge_dims(state.get_mps_site(min_pos),model.get_mpo(min_pos));
        var_front.set_edge_dims(state.get_mps_site(min_pos),model.get_mpo(min_pos));
        for(size_t pos = min_pos; pos < max_pos-1; pos++){
            auto & ene_curr = edges.get_eneL(pos);
            auto & ene_next = edges.get_eneL(pos+1);
            auto & var_curr = edges.get_varL(pos);
            auto & var_next = edges.get_varL(pos+1);
            ene_next = ene_curr.enlarge(state.get_mps_site(pos),model.get_mpo(pos));
            var_next = var_curr.enlarge(state.get_mps_site(pos),model.get_mpo(pos));
        }
    }
    {
        auto & ene_back = edges.get_eneR(max_pos);
        auto & var_back = edges.get_varR(max_pos);
        ene_back.set_edge_dims(state.get_mps_site(max_pos),model.get_mpo(max_pos));
        var_back.set_edge_dims(state.get_mps_site(max_pos),model.get_mpo(max_pos));
        for(size_t pos = max_pos; pos > min_pos; pos--){
            auto & ene_curr = edges.get_eneR(pos);
            auto & ene_prev = edges.get_eneR(pos-1);
            auto & var_curr = edges.get_varR(pos);
            auto & var_prev = edges.get_varR(pos-1);
            ene_prev = ene_curr.enlarge(state.get_mps_site(pos),model.get_mpo(pos));
            var_prev = var_curr.enlarge(state.get_mps_site(pos),model.get_mpo(pos));
        }

    }
}




void tools::finite::env::rebuild_active_edges(const class_state_finite & state, const class_model_finite & model, class_edges_finite & edges){
    if(not math::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw std::runtime_error(fmt::format("All lengths not equal: state {} | model {} | edges {}",state.get_length(), model.get_length(), edges.get_length() ));
    if(not math::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw std::runtime_error(fmt::format("All active sites are not equal: state {} | model {} | edges {}",state.active_sites, model.active_sites, edges.active_sites));

    edges.eject_inactive_edges();
    size_t min_pos = 0;
    size_t max_pos = state.get_length()-1;

    // If there are no active sites then we can build up until current position
    size_t posL_active = state.get_position();
    size_t posR_active = state.get_position();
    if(not edges.active_sites.empty()){
        posL_active = edges.active_sites.front();
        posR_active = edges.active_sites.back();
    }


    for(size_t pos = min_pos; pos < posL_active; pos++){
        auto & ene_curr = edges.get_eneL(pos);
        auto & ene_next = edges.get_eneL(pos+1);
        auto & var_curr = edges.get_varL(pos);
        auto & var_next = edges.get_varL(pos+1);
        if(not ene_next.has_block()) ene_next = ene_curr.enlarge(state.get_mps_site(pos),model.get_mpo(pos));
        if(not var_next.has_block()) var_next = var_curr.enlarge(state.get_mps_site(pos),model.get_mpo(pos));
    }

    for(size_t pos = max_pos; pos > posR_active; pos--){
        auto & ene_curr = edges.get_eneR(pos);
        auto & ene_prev = edges.get_eneR(pos-1);
        auto & var_curr = edges.get_varR(pos);
        auto & var_prev = edges.get_varR(pos-1);
        if(not ene_prev.has_block()) ene_prev = ene_curr.enlarge(state.get_mps_site(pos),model.get_mpo(pos));
        if(not var_prev.has_block()) var_prev = var_curr.enlarge(state.get_mps_site(pos),model.get_mpo(pos));
    }

    if(not edges.get_eneL(posL_active).has_block()) throw std::logic_error(fmt::format("Left active ene edge has undefined block"));
    if(not edges.get_varL(posL_active).has_block()) throw std::logic_error(fmt::format("Left active var edge has undefined block"));
    if(not edges.get_eneR(posR_active).has_block()) throw std::logic_error(fmt::format("Right active ene edge has undefined block"));
    if(not edges.get_varR(posR_active).has_block()) throw std::logic_error(fmt::format("Right active var edge has undefined block"));
}
