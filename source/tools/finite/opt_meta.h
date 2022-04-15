#pragma once
#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

enum class OptSolver;
enum class OptType;
enum class OptMode;
enum class OptInit;
enum class OptWhen;
enum class OptExit;
enum class OptRitz;

namespace tools::finite::opt {
    struct OptMeta {
        OptMode               optMode;
        OptSolver             optSolver;
        OptType               optType;
        OptInit               optInit;
        OptWhen               optWhen;
        OptRitz               optRitz;
        OptExit               optExit;
        size_t                max_sites        = 2;
        size_t                min_sites        = 1;
        long                  bond_lim         = -1;
        long                  max_problem_size = 0;
        long                  problem_size     = 0;
        bool                  retry            = true;
        std::optional<double> alpha_expansion  = std::nullopt;
        std::array<long, 3>   problem_dims     = {};
        std::vector<size_t>   chosen_sites     = {};
        std::string           label;
        std::optional<bool>   compress_otf            = std::nullopt; // Compress on the fly
        std::optional<double> eigs_max_tol            = std::nullopt;
        std::optional<int>    eigs_max_iter           = std::nullopt;
        std::optional<int>    eigs_max_ncv            = std::nullopt;
        std::optional<double> eigs_grad_tol           = std::nullopt;
        std::optional<double> bfgs_grad_tol           = std::nullopt;
        std::optional<int>    bfgs_max_iter           = std::nullopt;
        std::optional<int>    bfgs_max_rank           = std::nullopt;
        std::optional<double> bfgs_func_tol           = std::nullopt;
        std::optional<bool>   bfgs_eigenvalue_scaling = std::nullopt;

        OptMeta();
        explicit OptMeta(OptRitz ritz, OptMode mode);
        [[nodiscard]] bool should_proceed(OptExit previous_exit) const;
        void               validate() const;
    };
}
