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

int eig::solver::dsyevx(const real *matrix, size_type L, char range, int il, int iu, double vl, double vu, int m) {
    eig::log->trace("Starting eig dsyevx | range {} | i [{},{}] | v [{},{}] | m {}", range, il, iu, vl, vu, m);
    auto t_start = std::chrono::high_resolution_clock::now();
    auto A       = std::vector<real>(matrix, matrix + L * L);
    char jobz    = config.compute_eigvecs == Vecs::ON ? 'V' : 'N';
    int  info    = 0;
    int  n       = static_cast<int>(L);
    int  lda     = std::max(1, n);
    int  ldz     = std::max(1, n);
    int  liwork  = std::max(1, 5 * n);
    auto lifail  = std::max(1ul, static_cast<size_t>(n));

    if(range == 'I') m = iu - il + 1;
    m = std::min(m, n);

    int              m_found = m;
    double           lwork_query[1];
    std::vector<int> iwork(static_cast<size_t>(liwork));
    std::vector<int> ifail(lifail);

    auto &eigvals = result.get_eigvals<Form::SYMM>();
    auto &eigvecs = result.get_eigvecs<Form::SYMM, Type::REAL>();
    eigvals.resize(static_cast<size_t>(ldz));
    eigvecs.resize(static_cast<size_t>(ldz) * static_cast<size_t>(ldz)); // Docs claim ldz * m, but it segfaults when 'V' finds more than m eigvals

    info = LAPACKE_dsyevx_work(LAPACK_COL_MAJOR, jobz, range, 'U', static_cast<int>(L), A.data(), lda, vl, vu, il, iu, 2 * LAPACKE_dlamch('S'), &m_found,
                               eigvals.data(), eigvecs.data(), ldz, lwork_query, -1, iwork.data(), ifail.data());

    int lwork = static_cast<int>(lwork_query[0]);
    eig::log->trace(" lwork  = {}", lwork);
    eig::log->trace(" liwork = {}", liwork);
    std::vector<double> work(static_cast<size_t>(lwork));

    auto t_prep = std::chrono::high_resolution_clock::now();

    info = LAPACKE_dsyevx_work(LAPACK_COL_MAJOR, jobz, range, 'U', static_cast<int>(L), A.data(), lda, vl, vu, il, iu, 2 * LAPACKE_dlamch('S'), &m_found,
                               eigvals.data(), eigvecs.data(), ldz, work.data(), lwork, iwork.data(), ifail.data());

    auto t_total = std::chrono::high_resolution_clock::now();
    if(info == 0) {
        eig::log->trace("Found {} eigenvalues | requested {}", m_found, m);
        eigvals.resize(static_cast<size_t>(m_found));
        eigvecs.resize(static_cast<size_t>(m_found) * static_cast<size_t>(ldz));
        result.meta.eigvecsR_found = true;
        result.meta.eigvals_found  = true;
        result.meta.rows           = L;
        result.meta.cols           = m_found;
        result.meta.nev            = m_found;
        result.meta.nev_converged  = m_found;
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
