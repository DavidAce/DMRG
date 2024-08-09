#pragma once
#include "math/svd/config.h"
#include <optional>
#include <vector>
class StateFinite;
class ModelFinite;
class EdgesFinite;
class EnvEne;
class EnvVar;
enum class EnvExpandMode;
struct EnvExpansionResult;

namespace tools::finite::env {
    extern double             get_optimal_mixing_factor_ene_new(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                                const EdgesFinite &edges, EnvExpandMode envExpandMode);
    extern double             get_optimal_mixing_factor_var_new(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                                const EdgesFinite &edges, EnvExpandMode envExpandMode);
    extern double             get_optimal_mixing_factor_ene_old(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                                const EdgesFinite &edges, EnvExpandMode envExpandMode);
    extern double             get_optimal_mixing_factor_var_old(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                                const EdgesFinite &edges, EnvExpandMode envExpandMode);
    extern EnvExpansionResult expand_environment_backward(StateFinite &state, const ModelFinite &model, EdgesFinite &edges, EnvExpandMode envExpandMode,
                                                          svd::config svd_cfg);
    extern EnvExpansionResult expand_environment_forward(StateFinite &state, const ModelFinite &model, EdgesFinite &edges, EnvExpandMode envExpandMode,
                                                         svd::config svd_cfg);
    // extern ExpansionResult     expand_environment_forward_1site(StateFinite &state, const ModelFinite &model, EdgesFinite &edges, double alpha,
    //                                                       EnvExpandMode envExpandMode, svd::config svd_cfg);
    extern void assert_edges(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges);
    extern void assert_edges_var(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges);
    extern void assert_edges_ene(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges);
    extern void rebuild_edges(const StateFinite &state, const ModelFinite &model, EdgesFinite &edges);
    extern void rebuild_edges_var(const StateFinite &state, const ModelFinite &model, EdgesFinite &edges);
    extern void rebuild_edges_ene(const StateFinite &state, const ModelFinite &model, EdgesFinite &edges);
}