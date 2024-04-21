#include <complex>

#ifndef lapack_complex_float
    #define lapack_complex_float std::complex<float>
#endif
#ifndef lapack_complex_double
    #define lapack_complex_double std::complex<double>
#endif

#if defined(MKL_AVAILABLE)
    #include <mkl_lapacke.h>
#elif defined(OPENBLAS_AVAILABLE)
    #include <openblas/lapacke.h>
#else
    #include <lapacke.h>
#endif

#include "../log.h"
#include "../solver.h"
#include "math/cast.h"
#include <chrono>

int eig::solver::dgeev(real *matrix, size_type L) {
    eig::log->trace("Starting eig_dgeev (non-optimized)");
    auto  t_start      = std::chrono::high_resolution_clock::now();
    auto &eigvals_real = result.eigvals_real;
    auto &eigvals_imag = result.eigvals_imag;
    eigvals_real.resize(safe_cast<size_t>(L));
    eigvals_imag.resize(safe_cast<size_t>(L));
    std::vector<double> eigvecsR_tmp(safe_cast<size_t>(L * L));
    std::vector<double> eigvecsL_tmp(safe_cast<size_t>(L * L));

    int    info = 0;
    int    Lint = safe_cast<int>(L);
    double lwork_query;
    char   jobz = config.compute_eigvecs == Vecs::ON ? 'V' : 'N';

    info = LAPACKE_dgeev_work(LAPACK_COL_MAJOR, jobz, jobz, Lint, matrix, safe_cast<int>(L), eigvals_real.data(), eigvals_imag.data(), eigvecsL_tmp.data(),
                              Lint, eigvecsR_tmp.data(), safe_cast<int>(L), &lwork_query, -1);
    int                 lwork = (int) std::real(2.0 * lwork_query); // Make it twice as big for performance.
    std::vector<double> work(safe_cast<size_t>(lwork));
    auto                t_prep = std::chrono::high_resolution_clock::now();
    info = LAPACKE_dgeev_work(LAPACK_COL_MAJOR, jobz, jobz, Lint, matrix, safe_cast<int>(L), eigvals_real.data(), eigvals_imag.data(), eigvecsL_tmp.data(),
                              Lint, eigvecsR_tmp.data(), safe_cast<int>(L), work.data(), lwork);
    auto t_total = std::chrono::high_resolution_clock::now();
    if(info == 0) {
        result.meta.eigvecsR_found = true;
        result.meta.eigvecsL_found = true;
        result.meta.eigvals_found  = true;
        result.meta.rows           = L;
        result.meta.cols           = L;
        result.meta.nev            = Lint;
        result.meta.nev_converged  = Lint;
        result.meta.n              = L;
        result.meta.form           = Form::NSYM;
        result.meta.type           = Type::REAL;
        result.meta.time_prep      = std::chrono::duration<double>(t_prep - t_start).count();
        result.meta.time_total     = std::chrono::duration<double>(t_total - t_start).count();
    } else {
        throw std::runtime_error("LAPACK dgeev failed with error: " + std::to_string(info));
    }

    auto &eigvals  = result.eigvals_cplx;
    auto &eigvecsR = result.eigvecsR_cplx;
    auto &eigvecsL = result.eigvecsL_cplx;
    eigvals.resize(safe_cast<size_t>(L));
    eigvecsR.resize(safe_cast<size_t>(L * L));
    eigvecsL.resize(safe_cast<size_t>(L * L));

    // Copy eigenvalues
    for(size_t i = 0; i < safe_cast<size_t>(L); i++) eigvals[i] = std::complex<double>(eigvals_real[i], eigvals_imag[i]);
    // Release real/imag parts
    eigvals_real.clear();
    eigvals_imag.clear();

    // Copy eigenvectors
    auto count = 0ul;
    auto rows  = safe_cast<size_t>(L);
    auto cols  = safe_cast<size_t>(L);
    for(size_t i = 0; i < rows; i++) {
        size_t j = 0;
        while(j < cols) {
            if(eigvals_imag[j] == 0.0) {
                eigvecsR[count] = eigvecsR_tmp[i + j * rows]; //(i,j)
                eigvecsL[count] = eigvecsL_tmp[i + j * rows]; //(i,j)
                count++;
                j++;
            } else {
                eigvecsR[count] = std::complex<double>(eigvecsR_tmp[i + j * rows], eigvecsR_tmp[i + (j + 1) * rows]);
                eigvecsL[count] = std::complex<double>(eigvecsL_tmp[i + j * rows], eigvecsL_tmp[i + (j + 1) * rows]);
                count++;
                eigvecsR[count] = std::complex<double>(eigvecsR_tmp[i + j * rows], -eigvecsR_tmp[i + (j + 1) * rows]);
                eigvecsL[count] = std::complex<double>(eigvecsL_tmp[i + j * rows], -eigvecsL_tmp[i + (j + 1) * rows]);
                count++;
                j += 2;
            }
        }
    }

    return info;
}
