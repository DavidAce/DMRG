#pragma once
#include "enums.h"
#include <optional>
namespace eig {
    class settings {
        public:
        // Precision
        double    eigThreshold = 1e-12;
        size_type eigMaxIter   = 2000;
        size_type eigMaxNev    = 1;
        size_type eigMaxNcv    = 16;
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

        // Sanity checks
        void checkRitz();
        [[nodiscard]] std::string get_ritz_string() const;
    };
}
