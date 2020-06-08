#if __has_include(<mkl_lapacke.h>)
#include <mkl_lapacke.h>
#elif __has_include(<lapacke.h>)
#include <lapacke.h>
#endif

#include "math/eig.h"
#include "solver.h"

using namespace eig;

int eig::solver::dsyevd(const real* matrix, size_type L){
        eig::log->trace("Starting eig dsyevd");
    auto & eigvals = result.get_eigvals<Form::SYMM>();
    auto & eigvecs = result.get_eigvecs<Form::SYMM,Type::REAL>();
    eigvals.resize(static_cast<size_t>(L));
    eigvecs.resize(static_cast<size_t>(L*L));
    std::copy(matrix, matrix + L*L, eigvecs.begin());


    // Call lapack solver
    // For some reason the recommended lwork from netlib doesn't work. It's better to ask lapack with a query.
    // These nice values are inspired from armadillo. The prefactors give good performance.
    //    int lwork  = 2 * (1 + 6*L + 2*(L*L));
    //    int liwork = 3 * (3 + 5*L);
    int info   = 0;
    char jobz = config.compute_eigvecs == Vecs::ON ? 'V' : 'N';
    double lwork_query [1];
    int    liwork_query[1];

    info = LAPACKE_dsyevd_work(LAPACK_COL_MAJOR,jobz,'U',
                               static_cast<int>(L),
                               eigvecs.data(),
                               static_cast<int>(L),
                               eigvals.data(),
                               lwork_query,
                               -1,
                               liwork_query,
                               -1);

    int lwork     = static_cast<int>(2 * lwork_query[0]); //Make it twice as big for performance.
    int liwork    = static_cast<int>(3 * liwork_query[0]); //Make it thrice as big for performance.
    eig::log->trace(" lwork  = {}", lwork);
    eig::log->trace(" liwork = {}", liwork);

    std::vector<double> work  ( static_cast<size_t>(lwork));
    std::vector<int   > iwork ( static_cast<size_t>(liwork));

    info = LAPACKE_dsyevd_work(LAPACK_COL_MAJOR,jobz,'U',
                               static_cast<int>(L),
                               eigvecs.data(),
                               static_cast<int>(L),
                               eigvals.data(),
                               work.data(),
                               lwork,
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
        result.meta.type           = Type::REAL;
    }else{
        throw std::runtime_error("LAPACK dsyevd failed with error: " + std::to_string(info));
    }
    return info;
}
