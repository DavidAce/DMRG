#pragma once
#include <math/svd/config.h>
#include <optional>
#include <vector>
class StateFinite;
class ModelFinite;
class EdgesFinite;
enum class EnvExpandMode;
namespace tools::finite::env {
    extern std::vector<size_t> expand_environment(StateFinite &state, const ModelFinite &model, EdgesFinite &edges, std::optional<double> alpha,
                                                  EnvExpandMode envExpandMode, std::optional<svd::config> svd_cfg = std::nullopt);
    extern void                assert_edges(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges);
    extern void                assert_edges_var(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges);
    extern void                assert_edges_ene(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges);
    extern void                rebuild_edges(const StateFinite &state, const ModelFinite &model, EdgesFinite &edges);
    extern void                rebuild_edges_var(const StateFinite &state, const ModelFinite &model, EdgesFinite &edges);
    extern void                rebuild_edges_ene(const StateFinite &state, const ModelFinite &model, EdgesFinite &edges);
}