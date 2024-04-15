#pragma once
#include "math/svd/config.h"
#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

enum class OptSolver;
enum class OptFunc;
enum class OptAlgo;
enum class OptType;
enum class OptWhen;
enum class OptExit;
enum class OptRitz;

namespace tools::finite::opt {
    struct OptMeta {
        OptFunc                    optFunc;
        OptAlgo                    optAlgo;
        OptSolver                  optSolver;
        OptType                    optType;
        OptWhen                    optWhen;
        OptRitz                    optRitz;
        OptExit                    optExit;
        size_t                     max_sites        = 2;
        size_t                     min_sites        = 1;
        long                       max_problem_size = 0;
        long                       problem_size     = 0;
        std::optional<double>      alpha_expansion  = std::nullopt;
        std::array<long, 3>        problem_dims     = {};
        std::vector<size_t>        chosen_sites     = {};
        std::string                label;
        std::optional<bool>        compress_otf  = std::nullopt; // Compress on the fly
        std::optional<double>      subspace_tol  = std::nullopt;
        std::optional<double>      eigv_target   = std::nullopt; // AKA shift
        std::optional<double>      eigs_tol      = std::nullopt;
        std::optional<int>         eigs_nev      = std::nullopt;
        std::optional<int>         eigs_ncv      = std::nullopt;
        std::optional<int>         eigs_iter_max = std::nullopt;
        std::optional<svd::config> svd_cfg       = std::nullopt;

                           OptMeta();
        explicit           OptMeta(OptRitz ritz, OptFunc mode);
        [[nodiscard]] bool should_proceed(OptExit previous_exit) const;
        void               validate() const;
    };
}
