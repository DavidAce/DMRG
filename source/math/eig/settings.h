#pragma once
#include "enums.h"
#include <array>
#include <optional>
#include <vector>
namespace eig {

    class settings {
        public:
        // Solver library
        std::optional<Lib> lib = std::nullopt;
        // Precision
        std::optional<double>    tol     = std::nullopt;
        std::optional<double>    maxTime = std::nullopt;
        std::optional<int>       maxIter = std::nullopt;
        std::optional<size_type> maxNev  = std::nullopt;
        std::optional<size_type> maxNcv  = std::nullopt;

        // Solver properties
        std::optional<Form>    form            = std::nullopt;
        std::optional<Type>    type            = std::nullopt;
        std::optional<Side>    side            = std::nullopt;
        std::optional<Ritz>    ritz            = std::nullopt;
        std::optional<Storage> storage         = std::nullopt;
        std::optional<cplx>    sigma           = std::nullopt; // Sigma value for eigenvalue shift and shift-invert mode. Finds eigenvalues at the value "sigma"
        std::optional<Shinv>   shift_invert    = std::nullopt; // To find interior solutions. https://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html
        std::optional<Vecs>    compute_eigvecs = std::nullopt;
        std::optional<Dephase> remove_phase    = std::nullopt;
        std::optional<bool>    compress        = std::nullopt;

        struct init_t {
            void *ptr = nullptr; // Points to the initial guess for an eigenvector
            long  idx = 0;       // The index of the eigenpair;
        };
        std::vector<init_t>   initial_guess = {}; // Use these as initial guess. In arpack, the first element ptr can be used as the residual pointer "residp".
        std::string           tag;                // Handy when you are using many instances
        std::optional<size_t> loglevel = std::nullopt;
        std::optional<double> logTime  = std::nullopt;

        // Primme settings
        std::optional<std::string> primme_method               = std::nullopt;
        std::optional<std::string> primme_projection           = std::nullopt;
        std::optional<int>         primme_locking              = std::nullopt;
        std::optional<int>         primme_max_inner_iterations = std::nullopt; // Strongly recommend -1 or a fixed number like 100 here.
        std::optional<double>      primme_grad_tol             = std::nullopt;
        std::optional<int32_t>     primme_grad_iter            = std::nullopt;
        std::optional<double>      primme_grad_time            = std::nullopt;

        template<typename T>
        struct primme_gradient_t {
            std::vector<T>      mpo;
            std::vector<T>      envL;
            std::vector<T>      envR;
            std::array<long, 4> shape_mpo;
            std::array<long, 3> shape_envL;
            std::array<long, 3> shape_envR;
        };
        void *primme_gradient = nullptr;

        void *primme_extra = nullptr;

        void clear();
        // Sanity checks
        void                      checkRitz();
        [[nodiscard]] std::string get_ritz_string() const;
    };
}
