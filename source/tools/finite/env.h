#pragma once
#include <optional>
#include <vector>
class class_state_finite;
class class_model_finite;
class class_edges_finite;

namespace tools::finite::env {
    extern std::vector<size_t> expand_subspace(class_state_finite &state, const class_model_finite &model, class_edges_finite &edges,
                                               std::optional<double> alpha, long chi_lim,std::optional<double> svd_threshold = std::nullopt);

    extern void assert_edges(const class_state_finite &state, const class_model_finite &model, const class_edges_finite &edges);
    extern void assert_edges_var(const class_state_finite &state, const class_model_finite &model, const class_edges_finite &edges);
    extern void assert_edges_ene(const class_state_finite &state, const class_model_finite &model, const class_edges_finite &edges);
    extern void rebuild_edges(const class_state_finite &state, const class_model_finite &model, class_edges_finite &edges);
    extern void rebuild_edges_var(const class_state_finite &state, const class_model_finite &model, class_edges_finite &edges);
    extern void rebuild_edges_ene(const class_state_finite &state, const class_model_finite &model, class_edges_finite &edges);
}