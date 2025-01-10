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

int eig::solver::zgeev(cx64 *matrix, size_type L) {
    eig::log->trace("Starting eig_zgeev (non-optimized");
    auto  t_start  = std::chrono::high_resolution_clock::now();
    auto &eigvals  = result.get_eigvals<Form::NSYM>();
    auto &eigvecsR = result.get_eigvecs<Form::NSYM, Type::CPLX, Side::R>();
    auto &eigvecsL = result.get_eigvecs<Form::NSYM, Type::CPLX, Side::L>();
    eigvals.resize(safe_cast<size_t>(L));
    eigvecsR.resize(safe_cast<size_t>(L * L));
    eigvecsL.resize(safe_cast<size_t>(L * L));

    // int lwork   =  2*2*L;
    // For some reason the recommended lwork from netlib doesn't work. It's better to ask lapack with a query.
    int               lrwork = safe_cast<int>(2 * L);
    int               info   = 0;
    int               Lint   = safe_cast<int>(L);
    cx64              lwork_query;
    std::vector<fp64> rwork(safe_cast<size_t>(lrwork));
    auto              matrix_ptr      = reinterpret_cast<lapack_complex_double *>(const_cast<cx64 *>(matrix));
    auto              eigvals_ptr     = reinterpret_cast<lapack_complex_double *>(eigvals.data());
    auto              eigvecsL_ptr    = reinterpret_cast<lapack_complex_double *>(eigvecsL.data());
    auto              eigvecsR_ptr    = reinterpret_cast<lapack_complex_double *>(eigvecsR.data());
    auto              lwork_query_ptr = reinterpret_cast<lapack_complex_double *>(&lwork_query);
    char              jobz            = config.compute_eigvecs == Vecs::ON ? 'V' : 'N';

    info = LAPACKE_zgeev_work(LAPACK_COL_MAJOR, jobz, jobz, safe_cast<int>(L), matrix_ptr, safe_cast<int>(L), eigvals_ptr, eigvecsL_ptr, safe_cast<int>(L),
                              eigvecsR_ptr, Lint, lwork_query_ptr, -1, rwork.data());
    int                                lwork = (int) std::real(2.0 * lwork_query); // Make it twice as big for performance.
    std::vector<lapack_complex_double> work((unsigned long) lwork);
    auto                               t_prep = std::chrono::high_resolution_clock::now();
    info = LAPACKE_zgeev_work(LAPACK_COL_MAJOR, jobz, jobz, safe_cast<int>(L), matrix_ptr, safe_cast<int>(L), eigvals_ptr, eigvecsL_ptr, safe_cast<int>(L),
                              eigvecsR_ptr, Lint, work.data(), lwork, rwork.data());
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
        result.meta.type           = Type::CPLX;
        result.meta.time_prep      = std::chrono::duration<double>(t_prep - t_start).count();
        result.meta.time_total     = std::chrono::duration<double>(t_total - t_start).count();
    } else {
        throw std::runtime_error("LAPACK zgeev failed with error: " + std::to_string(info));
    }
    return info;
}
