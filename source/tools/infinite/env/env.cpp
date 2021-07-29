#include "../env.h"
#include <tensors/edges/EdgesInfinite.h>
#include <tensors/model/ModelInfinite.h>
#include <tensors/site/env/EnvEne.h>
#include <tensors/site/env/EnvVar.h>
#include <tensors/state/StateInfinite.h>

void tools::infinite::env::reset_edges(const StateInfinite &state, const ModelInfinite &model, EdgesInfinite &edges) {
    edges.get_ene().L.set_edge_dims(state.get_mps_siteA(), model.get_mpo_siteA());
    edges.get_ene().R.set_edge_dims(state.get_mps_siteB(), model.get_mpo_siteB());
    edges.get_var().L.set_edge_dims(state.get_mps_siteA(), model.get_mpo_siteA());
    edges.get_var().R.set_edge_dims(state.get_mps_siteB(), model.get_mpo_siteB());
}

void tools::infinite::env::enlarge_edges(const StateInfinite &state, const ModelInfinite &model, EdgesInfinite &edges) {
    const auto &ene = edges.get_ene();
    ene.L           = ene.L.enlarge(state.get_mps_siteA(), model.get_mpo_siteA());
    ene.R           = ene.R.enlarge(state.get_mps_siteB(), model.get_mpo_siteB());

    const auto &var = edges.get_var();
    var.L           = var.L.enlarge(state.get_mps_siteA(), model.get_mpo_siteA());
    var.R           = var.R.enlarge(state.get_mps_siteB(), model.get_mpo_siteB());
}
