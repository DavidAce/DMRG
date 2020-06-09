//
// Created by david on 2019-08-07.
//

#include <complex.h>
#undef I
#include<complex>
#define lapack_complex_float  std::complex<float>
#define lapack_complex_double std::complex<double>

#if __has_include(<mkl_lapacke.h>)
#include <mkl_lapacke.h>
#elif __has_include(<lapacke.h>)
#include <lapacke.h>
#endif

#include <Eigen/Core>
#include <math/svd.h>

template<typename Scalar>
std::tuple<svd::solver::MatrixType<Scalar>, svd::solver::VectorType<Scalar>,svd::solver::MatrixType<Scalar> , long>
svd::solver::do_svd_lapacke(const Scalar * mat_ptr, long rows, long cols, std::optional<long> rank_max){
    if(not rank_max.has_value()) rank_max = std::min(rows,cols);
    MatrixType<Scalar> A = Eigen::Map<const MatrixType<Scalar>> (mat_ptr,rows,cols);

    if (rows <= 0)              throw std::runtime_error("SVD error: rows() == 0");
    if (cols <= 0)              throw std::runtime_error("SVD error: cols() == 0");
    if (not A.allFinite())      throw std::runtime_error("SVD error: matrix has inf's or nan's");
    if (A.isZero(0))            throw std::runtime_error("SVD error: matrix is all zeros");

    int info   = 0;
    int rowsU  = static_cast<int>(rows);
    int colsU  = static_cast<int>(std::min(rows,cols));
    int rowsVT = static_cast<int>(std::min(rows,cols));
    int colsVT = static_cast<int>(cols);
    int sizeS  = static_cast<int>(std::min(rows,cols));
    int lda    = static_cast<int>(rows);
    int ldu    = static_cast<int>(rowsU);
    int ldvt   = static_cast<int>(rowsVT);


    MatrixType<Scalar> U(rowsU,colsU);
    VectorType<double> S(sizeS);
    MatrixType<Scalar> VT(rowsVT,colsVT);
    VectorType<Scalar> work(1);

    if constexpr (std::is_same<Scalar,double>::value){
        info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', static_cast<int>(rows),static_cast<int>(cols), A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), -1);
        int lwork  = (int) work(0);
        work.resize(lwork);
        info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', static_cast<int>(rows),static_cast<int>(cols), A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), lwork);
    }
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){
        int lrwork = (int) (5 * std::min(rows,cols));
        VectorType<double> rwork(lrwork);
        auto Ap      =  reinterpret_cast< lapack_complex_double *>(A.data());
        auto Up      =  reinterpret_cast< lapack_complex_double *>(U.data());
        auto VTp     =  reinterpret_cast< lapack_complex_double *>(VT.data());
        auto Wp_qry  =  reinterpret_cast< lapack_complex_double *>(work.data());

        info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', static_cast<int>(rows),static_cast<int>(cols), Ap, lda, S.data(), Up, ldu, VTp, ldvt, Wp_qry, -1,rwork.data());
        int lwork  = (int) std::real(work(0));
        work.resize(lwork);
        auto Wp  =  reinterpret_cast< lapack_complex_double *>(work.data());
        info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', static_cast<int>(rows),static_cast<int>(cols), Ap, lda, S.data(), Up, ldu, VTp, ldvt, Wp, lwork,rwork.data());
    }

    long max_size    =  std::min(S.size(),rank_max.value());
    long rank        = (S.head(max_size).array() >= SVDThreshold).count();
    if(rank == S.size()){
        truncation_error = 0;
    }else{
        truncation_error = S.tail(S.size()-rank).norm();
    }

    if (rank <= 0
        or not U.leftCols(rank).allFinite()
        or not S.head(rank).allFinite()
        or not VT.topRows(rank).allFinite() )
    {
        std::cerr   << "SVD error \n"
                    << "  svd_threshold     = " << SVDThreshold << '\n'
                    << "  Truncation Error = " << truncation_error << '\n'
                    << "  Rank             = " << rank << '\n'
                    << "  U all finite     : " << std::boolalpha << U.leftCols(rank).allFinite() << '\n'
                    << "  S all finite     : " << std::boolalpha << S.head(rank).allFinite() << '\n'
                    << "  V all finite     : " << std::boolalpha << VT.topRows(rank).allFinite() << '\n'
                    << "  Lapacke info     = " << info << '\n';
        throw std::runtime_error("SVD lapacke error:  Erroneous results");
    }

    return std::make_tuple(
            U.leftCols(rank),
            S.head(rank),
            VT.topRows(rank),
            rank
    );
}



//! \relates svd::class_SVD
//! \brief force instantiation of do_svd_lapacke for type 'double'
template std::tuple<svd::solver::MatrixType<double>, svd::solver::VectorType<double>,svd::solver::MatrixType<double> , long>
svd::solver::do_svd_lapacke(const double *, long, long, std::optional<long>);




using cplx = std::complex<double>;
//! \relates svd::class_SVD
//! \brief force instantiation of do_svd_lapacke for type 'std::complex<double>'
template std::tuple<svd::solver::MatrixType<cplx>, svd::solver::VectorType<cplx>,svd::solver::MatrixType<cplx> , long>
svd::solver::do_svd_lapacke(const cplx *, long, long, std::optional<long>);
