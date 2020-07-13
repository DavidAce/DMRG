#include <complex>
//#include <complex.h>
//#undef I

#ifndef lapack_complex_float
#define lapack_complex_float  std::complex<float>
#endif
#ifndef lapack_complex_double
#define lapack_complex_double std::complex<double>
#endif
#if __has_include(<mkl_lapacke.h>)
#include <mkl_lapacke.h>
#elif __has_include(<lapacke.h>)
#include <lapacke.h>
#endif

#include "math/eig.h"
#include "solver.h"

int eig::solver::zheevd(const cplx* matrix, size_type L){
    eig::log->trace("Starting eig zheevd. Eigvecs: {}", config.compute_eigvecs.value());

    auto & eigvals = result.get_eigvals<Form::SYMM>();
    auto & eigvecs = result.get_eigvecs<Form::SYMM,Type::CPLX>();
    eigvals.resize(static_cast<size_t>(L));
    eigvecs.resize(static_cast<size_t>(L*L));
    std::copy(matrix, matrix + L*L, eigvecs.begin());

    //These nice values are inspired from armadillo. The prefactors give good performance.
    int lwork  = static_cast<int>(2 * (2*L + L*L));
    int lrwork = static_cast<int>(2 * (1 + 5*L + 2*(L*L)));
    int liwork = static_cast<int>(3 * (3 + 5*L));
    int info   = 0;
    std::vector<cplx> work    ( static_cast<size_t>(lwork) );
    std::vector<real> rwork   ( static_cast<size_t>(lrwork) );
    std::vector<int   > iwork ( static_cast<size_t>(liwork) );
    char jobz = config.compute_eigvecs == Vecs::ON ? 'V' : 'N';

    info = LAPACKE_zheevd_work(LAPACK_COL_MAJOR,jobz,'U',
                               static_cast<int>(L),
                               reinterpret_cast< lapack_complex_double *>(eigvecs.data()),
                               static_cast<int>(L),
                               eigvals.data(),
                               reinterpret_cast< lapack_complex_double *>(work.data()),
                               lwork,
                               rwork.data(),
                               lrwork,
                               iwork.data(),
                               liwork);


    if (info == 0){
        result.meta.eigvecsR_found = true;
        result.meta.eigvals_found  = true;
        result.meta.rows           = L;
        result.meta.cols           = L;
        result.meta.nev            = L;
        result.meta.n              = L;
        result.meta.form           = Form::SYMM;
        result.meta.type           = Type::CPLX;
    }else{
        throw std::runtime_error("LAPACK zheevd failed with error: " + std::to_string(info));
    }
    return info;
}

