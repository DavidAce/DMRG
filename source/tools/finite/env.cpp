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
    size_t min_pos = 0;
    size_t max_pos = state.get_length()-1;
    {
        edges.ene.front().L.set_edge_dims(state.get_mps(min_pos),model.get_mpo(min_pos));
        edges.var.front().L.set_edge_dims(state.get_mps(min_pos),model.get_mpo(min_pos));
        for(size_t pos = min_pos; pos < max_pos; pos++){
            auto & ene_curr = edges.get_ene(pos);
            auto & ene_next = edges.get_ene(pos+1);
            auto & var_curr = edges.get_var(pos);
            auto & var_next = edges.get_var(pos+1);
            ene_next.L.block = ene_curr.L.enlarge(state.get_mps(pos),model.get_mpo(pos));
            var_next.L.block = var_curr.L.enlarge(state.get_mps(pos),model.get_mpo(pos));
        }
    }
    {
        edges.ene.back().R.set_edge_dims(state.get_mps(max_pos),model.get_mpo(max_pos));
        edges.var.back().R.set_edge_dims(state.get_mps(max_pos),model.get_mpo(max_pos));
        for(size_t pos = max_pos; pos >= min_pos; pos--){
            auto & ene_curr = edges.get_ene(pos);
            auto & ene_prev = edges.get_ene(pos-1);
            auto & var_curr = edges.get_var(pos);
            auto & var_prev = edges.get_var(pos-1);
            ene_prev.L.block = ene_curr.L.enlarge(state.get_mps(pos),model.get_mpo(pos));
            var_prev.L.block = var_curr.L.enlarge(state.get_mps(pos),model.get_mpo(pos));
        }

    }



//    if(state.MPS_L.size() != state.MPO_L.size())
//        throw std::runtime_error(fmt::format("Size mismatch in MPSL and MPOL: {} != {}", state.MPS_L.size(), state.MPO_L.size()));
//    if(state.MPS_R.size() != state.MPO_R.size())
//        throw std::runtime_error(fmt::format("Size mismatch in MPSR and MPOR: {} != {}", state.MPS_R.size(), state.MPO_R.size()));
    // Generate new environments

    {
//        state.ENV_L.clear();
//        state.ENV2_L.clear();
//
//        auto mpsL_it = state.MPS_L.begin();
//        auto mpoL_it = state.MPO_L.begin();
//        auto ENV_L   = class_environment("L", *mpsL_it, *mpoL_it->get());     // Initialized envs
//        auto ENV2_L  = class_environment_var("L", *mpsL_it, *mpoL_it->get()); // Initialized envs
//        while(mpsL_it != state.MPS_L.end() and mpoL_it != state.MPO_L.end()) {
//            if(ENV_L.has_nan()) throw std::runtime_error("ENV_L " + std::to_string(ENV_L.get_position()) + " has NAN's");
//            state.ENV_L.emplace_back(ENV_L);
//            state.ENV2_L.emplace_back(ENV2_L);
//            if(mpsL_it->get_position() != state.ENV_L.back().get_position())
//                throw std::runtime_error(fmt::format("Size mismatch in MPSL and ENVL: {} != {}", mpsL_it->get_position(), state.ENV_L.back().get_position()));
//            if(mpsL_it->get_chiL() != state.ENV_L.back().block.dimension(0))
//                throw std::runtime_error(
//                    fmt::format("Size mismatch in MPSL and ENVL dimensions {} != {}", mpsL_it->get_chiL(), state.ENV_L.back().block.dimension(2)));
//
//            ENV_L  = ENV_L.enlarge(*mpsL_it, *mpoL_it->get());
//            ENV2_L = ENV2_L.enlarge(*mpsL_it, *mpoL_it->get());
//            mpsL_it++;
//            mpoL_it++;
//        }
    }

    {
//        state.ENV_R.clear();
//        state.ENV2_R.clear();
//
//        auto mpsR_it = state.MPS_R.rbegin();
//        auto mpoR_it = state.MPO_R.rbegin();
//        auto ENV_R   = class_environment("R", *mpsR_it, *mpoR_it->get());
//        auto ENV2_R  = class_environment_var("R", *mpsR_it, *mpoR_it->get());
//        while(mpsR_it != state.MPS_R.rend() and mpoR_it != state.MPO_R.rend()) {
//            if(ENV_R.has_nan()) throw std::runtime_error("ENV_R " + std::to_string(ENV_R.get_position()) + " has NAN's");
//            state.ENV_R.emplace_front(ENV_R);
//            state.ENV2_R.emplace_front(ENV2_R);
//            if(mpsR_it->get_position() != state.ENV_R.front().get_position())
//                throw std::runtime_error(fmt::format("Size mismatch in MPSR and ENVR: {} != {}", mpsR_it->get_position(), state.ENV_R.front().get_position()));
//            if(mpsR_it->get_chiR() != state.ENV_R.front().block.dimension(0))
//                throw std::runtime_error(
//                    fmt::format("Size mismatch in MPSR and ENVR dimensions {} != {}", mpsR_it->get_chiR(), state.ENV_R.front().block.dimension(2)));
//            ENV_R  = ENV_R.enlarge(*mpsR_it, *mpoR_it->get());
//            ENV2_R = ENV2_R.enlarge(*mpsR_it, *mpoR_it->get());
//            mpsR_it++;
//            mpoR_it++;
//        }
    }
}