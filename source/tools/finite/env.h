#pragma once

class class_state_finite;
class class_model_finite;
class class_edges_finite;

namespace tools::finite::env {
    //    extern void rebuild_all_edges(const class_state_finite & state, const class_model_finite & model, class_edges_finite & edges);
    extern void assert_edges(const class_state_finite &state, const class_model_finite &model, const class_edges_finite &edges);
    extern void assert_edges_var(const class_state_finite &state, const class_model_finite &model, const class_edges_finite &edges);
    extern void assert_edges_ene(const class_state_finite &state, const class_model_finite &model, const class_edges_finite &edges);
    extern void rebuild_edges(const class_state_finite &state, const class_model_finite &model, class_edges_finite &edges);
    extern void rebuild_edges_var(const class_state_finite &state, const class_model_finite &model, class_edges_finite &edges);
    extern void rebuild_edges_ene(const class_state_finite &state, const class_model_finite &model, class_edges_finite &edges);
}