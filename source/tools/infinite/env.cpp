//
// Created by david on 2019-10-24.
//

#include <tools/nmspc_tools.h>
#include <state/class_infinite_state.h>
#include <state/class_environment.h>
#include <state/class_mps_2site.h>
#include <state/class_mps_site.h>
#include <model/class_model_base.h>
void tools::infinite::env::initialize  (class_infinite_state & state){
    state.Lblock = std::make_unique<class_environment>("L");
    state.Rblock = std::make_unique<class_environment>("R");
    state.Lblock2 = std::make_unique<class_environment_var>("L");
    state.Rblock2 = std::make_unique<class_environment_var>("R");

    state.Lblock->set_edge_dims(*state.MPS->MPS_A, state.HA->MPO());
    state.Rblock->set_edge_dims(*state.MPS->MPS_B, state.HB->MPO());
    state.Lblock2->set_edge_dims(*state.MPS->MPS_A, state.HA->MPO());
    state.Rblock2->set_edge_dims(*state.MPS->MPS_B, state.HB->MPO());
}
