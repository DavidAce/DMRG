//
// Created by david on 2018-10-30.
//

#pragma once
#include "math/eig/enums.h"
#include "math/eig/settings.h"
#include "math/eig/sfinae.h"
#include "math/eig/solution.h"
#include <complex>
#include <vector>

class class_tic_toc;

namespace eig {
    template<typename MatrixType>
    class arpack_solver {
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
        std::unique_ptr<class_tic_toc> t_sol;

        std::unique_ptr<class_tic_toc> t_tot;
        std::unique_ptr<class_tic_toc> t_pre;
        std::unique_ptr<class_tic_toc> t_mul;
        std::unique_ptr<class_tic_toc> t_fnd;

        void eigs();
        template<typename Derived>
        void find_solution(Derived &solver, eig::size_type nev);
        template<typename Derived>
        void find_solution_rc(Derived &solver);
        template<typename Derived>
        void copy_solution(Derived &solver);

        MatrixType &   matrix;
        eig::settings &config;
        eig::solution &result;
        Scalar *       residual = nullptr;
        arpack_solver(MatrixType &matrix_, eig::settings &config_, eig::solution &result_, Scalar *residual_);
    };

}
