//
// Created by david on 2018-10-30.
//

#pragma once

#include "../enums.h"
#include "../settings.h"
#include "../solution.h"
#include "../sfinae.h"
#include <complex>
#include <general/class_tic_toc.h>
#include <vector>

#define profile_arpack 0



namespace eig{
    template<typename MatrixType>
    class arpack_solver {
    private:
        int nev_internal;
        int ncv_internal;


        void eigs_sym();
        void eigs_nsym();
        void eigs_comp();

    public:
        using Scalar = typename MatrixType::Scalar;
        class_tic_toc t_sol;
        class_tic_toc t_get;
        class_tic_toc t_sub;
        class_tic_toc t_all;

        void eigs();
        template<typename Derived>
        void find_solution(Derived &solver, eig::size_type nev);

        MatrixType &   matrix;
        eig::settings &config;
        eig::solution &result;
        Scalar *       residual;
        arpack_solver(MatrixType &matrix_, eig::settings &config_, eig::solution &result_, Scalar *residual_);

        template<typename Derived>
        void copy_solution(Derived &solver) {
            // In this function we copy the solution "naively"
            eig::log->trace("Copying solution");
            using eigval_type = std::remove_pointer_t<decltype(solver.RawEigenvalues())>;
            using eigvec_type = std::remove_pointer_t<decltype(solver.RawEigenvectors())>;
            auto  eigvecsize = result.meta.rows * result.meta.cols;
            auto  eigvalsize = result.meta.cols;
            auto  eigvalsize_t = static_cast<size_t>(eigvalsize);
            auto  eigvecsize_t = static_cast<size_t>(eigvecsize);
            constexpr auto eigval_has_imag_separately = eig::sfinae::has_RawEigenvaluesImag_v<Derived>;
            constexpr auto eigval_is_cplx = std::is_same_v<cplx,eigval_type>;
            constexpr auto eigval_is_real = std::is_same_v<real,eigval_type>;
            constexpr auto eigvec_is_cplx = std::is_same_v<cplx,eigvec_type>;
            constexpr auto eigvec_is_real = std::is_same_v<real,eigvec_type>;


            if constexpr (eigval_has_imag_separately){
                eig::log->trace("Copying eigenvalues from separate real and imaginary buffers");
                result.eigvals_imag.resize(eigvalsize_t);
                result.eigvals_real.resize(eigvalsize_t);
                std::copy(solver.RawEigenvaluesImag(), solver.RawEigenvaluesImag() + eigvalsize, result.eigvals_imag.begin());
                std::copy(solver.RawEigenvaluesReal(), solver.RawEigenvaluesReal() + eigvalsize, result.eigvals_real.begin());
            }
            if constexpr(eigval_is_real){
                eig::log->trace("Copying real eigenvalues");
                result.eigvals_real.resize(eigvalsize_t);
                std::copy(solver.RawEigenvalues(), solver.RawEigenvalues() + eigvalsize, result.eigvals_real.begin());
            }
            if constexpr(eigval_is_cplx){
                eig::log->trace("Copying complex eigenvalues");
                result.eigvals_cplx.resize(eigvalsize_t);
                std::copy(solver.RawEigenvalues(), solver.RawEigenvalues() + eigvalsize, result.eigvals_cplx.begin());
            }
            if constexpr(eigvec_is_real){
                eig::log->trace("Copying real eigenvectors");
                result.eigvecsR_real.resize(eigvecsize_t);
                std::copy(solver.RawEigenvectors(), solver.RawEigenvectors() + eigvecsize, result.eigvecsR_real.begin());
            }

            if constexpr(eigvec_is_cplx){
                eig::log->trace("Copying complex eigenvectors");
                if(config.side == Side::L){
                    result.eigvecsL_cplx.resize(eigvecsize_t);
                    std::copy(solver.RawEigenvectors(), solver.RawEigenvectors() + eigvecsize, result.eigvecsL_cplx.begin());
                }else{
                    result.eigvecsR_cplx.resize(eigvecsize_t);
                    std::copy(solver.RawEigenvectors(), solver.RawEigenvectors() + eigvecsize, result.eigvecsR_cplx.begin());
                }
            }
        }

    };

}
