#pragma once
#include "enums.h"
#include <optional>
namespace eig {
    class settings {
        public:
        // Precision
        std::optional<double>    eigThreshold = std::nullopt;
        std::optional<int>       eigMaxIter   = std::nullopt;
        std::optional<size_type> eigMaxNev    = std::nullopt;
        std::optional<size_type> eigMaxNcv    = std::nullopt;
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

        void clear();
        // Sanity checks
        void                      checkRitz();
        [[nodiscard]] std::string get_ritz_string() const;
    };
}
