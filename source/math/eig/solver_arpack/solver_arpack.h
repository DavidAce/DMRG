#pragma once
#include "../enums.h"
#include "../settings.h"
#include "../sfinae.h"
#include "../solution.h"
#include <complex>
#include <memory>
#include <vector>
namespace tid {
    class ur;
}

namespace eig {
    template<typename MatrixType>
    class solver_arpack {
        private:
        int nev_internal;
        int ncv_internal;

        void eigs_sym();
        void eigs_nsym();
        void eigs_comp();

        void eigs_sym_rc();
        void eigs_nsym_rc();
        void eigs_comp_rc();

        public:
        using Scalar = typename MatrixType::Scalar;
        std::unique_ptr<tid::ur> t_tot;
        std::unique_ptr<tid::ur> t_pre;
        std::unique_ptr<tid::ur> t_mul;
        std::unique_ptr<tid::ur> t_fnd;

        void eigs();
        template<typename Derived>
        void find_solution(Derived &solver, eig::size_type nev);
        template<typename Derived>
        void find_solution_rc(Derived &solver);
        template<typename Derived>
        void copy_solution(Derived &solver);

        MatrixType    &matrix;
        eig::settings &config;
        eig::solution &result;
        Scalar        *residual = nullptr;
        solver_arpack(MatrixType &matrix_, eig::settings &config_, eig::solution &result_, Scalar *residual_);
    };

}
