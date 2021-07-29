#pragma once
#include "enums.h"
#include <optional>
#include <vector>
namespace eig {

    class settings {
        public:
        // Solver library
        std::optional<Lib>         lib           = std::nullopt;
        std::optional<std::string> primme_method = std::nullopt;
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
        int                    ncv_x_factor    = 2;
        std::vector<int>       iter_ncv_x;
        std::vector<double>    time_tol_x10;
        void                   clear();
        // Sanity checks
        void                      checkRitz();
        [[nodiscard]] std::string get_ritz_string() const;
    };
}
