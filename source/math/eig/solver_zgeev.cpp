
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
#elif __has_include(<openblas/lapacke.h>)
#include <openblas/lapacke.h>
#else
#include <lapacke.h>
#endif

#include "solver.h"

int eig::solver::zgeev(const cplx* matrix, size_type L){
    eig::log->trace("Starting eig_zgeev");
    auto & eigvals  = result.get_eigvals<Form::NSYM>();
    auto & eigvecsR = result.get_eigvecs<Form::NSYM,Type::CPLX, Side::R>();
    auto & eigvecsL = result.get_eigvecs<Form::NSYM,Type::CPLX, Side::L>();
    eigvals.resize(static_cast<size_t>(L));
    eigvecsR.resize(static_cast<size_t>(L*L));
    eigvecsL.resize(static_cast<size_t>(L*L));

    // int lwork   =  2*2*L;
    // For some reason the recommended lwork from netlib doesn't work. It's better to ask lapack with a query.
    int lrwork  =  static_cast<int>(2*L);
    int info   = 0;
    cplx lwork_query;
    std::vector<real> rwork  ( static_cast<size_t>(lrwork));
    auto matrix_ptr         =  reinterpret_cast< lapack_complex_double *>(const_cast<cplx *>(matrix));
    auto eigvals_ptr        =  reinterpret_cast< lapack_complex_double *>(eigvals.data());
    auto eigvecsL_ptr       =  reinterpret_cast< lapack_complex_double *>(eigvecsL.data());
    auto eigvecsR_ptr       =  reinterpret_cast< lapack_complex_double *>(eigvecsR.data());
    auto lwork_query_ptr    =  reinterpret_cast< lapack_complex_double *>(&lwork_query);
    char jobz = config.compute_eigvecs == Vecs::ON ? 'V' : 'N';

    info = LAPACKE_zgeev_work(LAPACK_COL_MAJOR,jobz,jobz,
                              static_cast<int>(L),
                              matrix_ptr,
                              static_cast<int>(L),
                              eigvals_ptr,
                              eigvecsL_ptr,
                              static_cast<int>(L),
                              eigvecsR_ptr,
                              static_cast<int>(L),
                              lwork_query_ptr,
                              -1,
                              rwork.data());
    int lwork = (int) std::real(2.0*lwork_query); //Make it twice as big for performance.
    std::vector<lapack_complex_double> work  ( (unsigned long)lwork );

    info = LAPACKE_zgeev_work(LAPACK_COL_MAJOR,jobz,jobz,
                              static_cast<int>(L),
                              matrix_ptr,
                              static_cast<int>(L),
                              eigvals_ptr,
                              eigvecsL_ptr,
                              static_cast<int>(L),
                              eigvecsR_ptr,
                              static_cast<int>(L),
                              work.data(),
                              lwork,
                              rwork.data());



    if (info == 0){
        result.meta.eigvecsR_found = true;
        result.meta.eigvecsL_found = true;
        result.meta.eigvals_found  = true;
        result.meta.rows           = L;
        result.meta.cols           = L;
        result.meta.nev            = L;
        result.meta.n              = L;
        result.meta.form           = Form::NSYM;
        result.meta.type           = Type::CPLX;
    }else{
        throw std::runtime_error("LAPACK zgeev failed with error: " + std::to_string(info));
    }
    return info;
}



