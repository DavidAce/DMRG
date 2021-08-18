#pragma once
#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

enum class OptSpace;
enum class OptType;
enum class OptMode;
enum class OptInit;
enum class OptWhen;
enum class OptExit;
enum class OptRitz;

namespace tools::finite::opt {
    struct OptMeta {
        OptMode               optMode;
        OptSpace              optSpace;
        OptType               optType;
        OptInit               optInit;
        OptWhen               optWhen;
        OptRitz               optRitz;
        OptExit               optExit;
        size_t                max_sites        = 2;
        size_t                min_sites        = 1;
        long                  max_problem_size = 0;
        long                  problem_size     = 0;
        bool                  retry            = true;
        std::optional<double> alpha_expansion  = std::nullopt;
        std::array<long, 3>   problem_dims     = {};
        std::vector<size_t>   chosen_sites     = {};
        std::string           label;
        std::optional<double> max_grad_tolerance = std::nullopt;
        std::optional<int>    ceres_max_num_iterations = std::nullopt;
        std::optional<int>    ceres_max_lbfgs_rank     = std::nullopt;
        std::optional<double> ceres_function_tolerance = std::nullopt;
        std::optional<bool>   ceres_eigenvalue_scaling = std::nullopt;

        OptMeta();
        explicit OptMeta(OptRitz ritz);
        [[nodiscard]] bool should_proceed(OptExit previous_exit) const;
    };
}
