#include <complex>

#ifndef lapack_complex_float
    #define lapack_complex_float std::complex<float>
#endif
#ifndef lapack_complex_double
    #define lapack_complex_double std::complex<double>
#endif

#if __has_include(<mkl_lapacke.h>)
    #include <mkl_lapacke.h>
#elif __has_include(<openblas/lapacke.h>)
    #include <openblas/lapacke.h>
#else
    #include <lapacke.h>
#endif
#include "../log.h"
#include "../solver.h"
#include <chrono>

using namespace eig;

int eig::solver::dsyevx(const real *matrix, size_type L, char range, int il, int iu, double vl, double vu) {
    eig::log->trace("Starting eig dsyevx");
    auto  t_start = std::chrono::high_resolution_clock::now();
    auto &eigvals = result.get_eigvals<Form::SYMM>();
    auto &eigvecs = result.get_eigvecs<Form::SYMM, Type::REAL>();
    auto  A       = std::vector<real>(matrix, matrix + L * L);
    // Call lapack solver
    // For some reason the recommended lwork from netlib doesn't work. It's better to ask lapack with a query.
    // These nice values are inspired from armadillo. The prefactors give good performance.
    //    int lwork  = 2 * (1 + 6*L + 2*(L*L));
    //    int liwork = 3 * (3 + 5*L);
    char             jobz   = config.compute_eigvecs == Vecs::ON ? 'V' : 'N';
    int              info   = 0;
    int              m      = iu - il + 1;
    int              lda    = static_cast<int>(L);
    int              ldz    = std::max(1, lda);
    int              liwork = 5 * lda;
    double           lwork_query[1];
    std::vector<int> ifail(static_cast<unsigned long>(L));
    std::vector<int> iwork(static_cast<size_t>(liwork));

    eigvals.resize(static_cast<size_t>(m));
    eigvecs.resize(static_cast<size_t>(L * m));

    info = LAPACKE_dsyevx_work(LAPACK_COL_MAJOR, jobz, range, 'U', static_cast<int>(L), A.data(), lda, vl, vu, il, iu, 2 * LAPACKE_dlamch('S'), &m,
                               eigvals.data(), eigvecs.data(), ldz, lwork_query, -1, iwork.data(), ifail.data());

    int lwork = static_cast<int>(2 * lwork_query[0]); // Make it twice as big for performance.
    eig::log->trace(" lwork  = {}", lwork);
    eig::log->trace(" liwork = {}", liwork);

    std::vector<double> work(static_cast<size_t>(lwork));
    auto                t_prep = std::chrono::high_resolution_clock::now();
    info = LAPACKE_dsyevx_work(LAPACK_COL_MAJOR, jobz, range, 'U', static_cast<int>(L), A.data(), lda, vl, vu, il, iu, 2 * LAPACKE_dlamch('S'), &m,
                               eigvals.data(), eigvecs.data(), ldz, work.data(), lwork, iwork.data(), ifail.data());

    auto t_total = std::chrono::high_resolution_clock::now();
    if(info == 0) {
        result.meta.eigvecsR_found = true;
        result.meta.eigvals_found  = true;
        result.meta.rows           = L;
        result.meta.cols           = m;
        result.meta.nev            = m;
        result.meta.nev_converged  = m;
        result.meta.n              = L;
        result.meta.form           = Form::SYMM;
        result.meta.type           = Type::REAL;
        result.meta.time_prep      = std::chrono::duration<double>(t_prep - t_start).count();
        result.meta.time_total     = std::chrono::duration<double>(t_total - t_start).count();
    } else {
        throw std::runtime_error("LAPACK dsyevx failed with error: " + std::to_string(info));
    }
    return info;
}
