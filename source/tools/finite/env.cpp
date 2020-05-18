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

void tools::finite::env::rebuild_edges(const class_state_finite & state, const class_model_finite & model, class_edges_finite & edges) {
    tools::log->trace("Rebuilding edges");
    if(not math::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw std::runtime_error(fmt::format("All lengths not equal: state {} | model {} | edges {}",state.get_length(), model.get_length(), edges.get_length() ));

    // Note that the last environment should always be empty!
    // Going left to right, that will be the back of the list.
    // Going right to left, that will be the front of the list.
    size_t min_pos = 0;
    size_t max_pos = state.get_length()-1;
    {
        const auto & ene_front = edges.get_ene(min_pos);
        const auto & var_front = edges.get_var(min_pos);
        ene_front.L.set_edge_dims(state.get_mps_site(min_pos),model.get_mpo(min_pos));
        var_front.L.set_edge_dims(state.get_mps_site(min_pos),model.get_mpo(min_pos));
        for(size_t pos = min_pos; pos < max_pos-1; pos++){
            const auto & ene_curr = edges.get_ene(pos);
            const auto & ene_next = edges.get_ene(pos+1);
            const auto & var_curr = edges.get_var(pos);
            const auto & var_next = edges.get_var(pos+1);
            ene_next.L.block = ene_curr.L.enlarge(state.get_mps_site(pos),model.get_mpo(pos));
            var_next.L.block = var_curr.L.enlarge(state.get_mps_site(pos),model.get_mpo(pos));
        }
    }
    {
        const auto & ene_back = edges.get_ene(max_pos);
        const auto & var_back = edges.get_var(max_pos);
        ene_back.R.set_edge_dims(state.get_mps_site(max_pos),model.get_mpo(max_pos));
        var_back.R.set_edge_dims(state.get_mps_site(max_pos),model.get_mpo(max_pos));
        for(size_t pos = max_pos; pos > min_pos; pos--){
            const auto & ene_curr = edges.get_ene(pos);
            const auto & ene_prev = edges.get_ene(pos-1);
            const auto & var_curr = edges.get_var(pos);
            const auto & var_prev = edges.get_var(pos-1);
            ene_prev.R.block = ene_curr.R.enlarge(state.get_mps_site(pos),model.get_mpo(pos));
            var_prev.R.block = var_curr.R.enlarge(state.get_mps_site(pos),model.get_mpo(pos));
        }

    }
}