//
// Created by david on 2019-10-24.
//

#include "env.h"
#include <tensors/state/class_mps_2site.h>
#include <tensors/state/class_state_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/edges/class_env_ene.h>
#include <tensors/edges/class_env_var.h>


void tools::infinite::env::enlarge_edges(const class_state_infinite & state, const class_model_infinite & model, class_edges_infinite & edges){
    const auto & ene = edges.get_ene();
    ene.L = ene.L.enlarge(state.get_mps_siteA(), model.get_mpo_siteA());
    ene.R = ene.R.enlarge(state.get_mps_siteB(), model.get_mpo_siteB());

    const auto & var = edges.get_var();
    var.L = var.L.enlarge(state.get_mps_siteA(), model.get_mpo_siteA());
    var.R = var.R.enlarge(state.get_mps_siteB(), model.get_mpo_siteB());
}
