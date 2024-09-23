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
enum class OptAlgo;
enum class OptRitz;
enum class EnvExpandMode;
struct EnvExpansionResult;

namespace tools::finite::env {
    extern EnvExpansionResult get_optimally_mixed_block(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                        const EdgesFinite &edges, OptAlgo algo, OptRitz ritz);
    namespace internal {
        extern EnvExpansionResult get_optimally_mixed_block_H1(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                               const EdgesFinite &edges, OptRitz ritz);
        extern EnvExpansionResult get_optimally_mixed_block_H2(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                               const EdgesFinite &edges, OptRitz ritz);
        extern EnvExpansionResult get_optimally_mixed_block_VarH(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                                 const EdgesFinite &edges, OptRitz ritz);
        extern EnvExpansionResult get_optimally_mixed_block_GsiH(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                                 const EdgesFinite &edges, OptRitz ritz);
    }

    extern std::tuple<double, double, double> get_optimal_mixing_factors(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                                         const EdgesFinite &edges, EnvExpandMode envExpandMode);
    extern std::array<double, 2> get_optimal_mixing_factor_ene(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                               const EdgesFinite &edges, EnvExpandMode envExpandMode);
    extern std::array<double, 2> get_optimal_mixing_factor_var(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                               const EdgesFinite &edges, EnvExpandMode envExpandMode);
    extern double                get_optimal_mixing_factor_ene_old(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                                   const EdgesFinite &edges, EnvExpandMode envExpandMode);
    extern double                get_optimal_mixing_factor_var_old(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                                   const EdgesFinite &edges, EnvExpandMode envExpandMode);
    extern EnvExpansionResult    expand_environment_backward(StateFinite &state, const ModelFinite &model, EdgesFinite &edges, EnvExpandMode envExpandMode,
                                                             svd::config svd_cfg);
    extern EnvExpansionResult    expand_environment_reversed(StateFinite &state, const ModelFinite &model, EdgesFinite &edges, EnvExpandMode envExpandMode,
                                                             svd::config svd_cfg);
    extern EnvExpansionResult    expand_environment_forward_old(StateFinite &state, const ModelFinite &model, EdgesFinite &edges, EnvExpandMode envExpandMode,
                                                                svd::config svd_cfg);
    extern EnvExpansionResult    expand_environment_forward_nsite(StateFinite &state, const ModelFinite &model, EdgesFinite &edges, OptAlgo algo, OptRitz ritz,
                                                                  svd::config svd_cfg);
    extern EnvExpansionResult    expand_environment_forward_1site(StateFinite &state, const ModelFinite &model, EdgesFinite &edges, EnvExpandMode envExpandMode,
                                                                  svd::config svd_cfg);
    extern EnvExpansionResult expand_environment_forward_active(StateFinite &state, const ModelFinite &model, EdgesFinite &edges, EnvExpandMode envExpandMode,
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