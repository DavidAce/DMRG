//
// Created by david on 2019-08-07.
//

#include <complex.h>
#undef I

#ifdef MKL_AVAILABLE
#ifndef MKL_Complex8
#define MKL_Complex8 std::complex<float>
#endif
#ifndef MKL_Complex16
#define MKL_Complex16 std::complex<double>
#endif
#ifndef ComplexFloat_
#define ComplexFloat_  MKL_Complex8
#endif
#ifndef ComplexDouble_
#define ComplexDouble_ MKL_Complex16
#endif
#include <mkl_lapacke.h>
#else
#define ComplexFloat_  __complex__ float
#define ComplexDouble_ __complex__ double
#include <lapacke.h>
#endif

#include <math/class_svd_wrapper.h>
#include <Eigen/Core>
//
//template<typename Scalar>
//
//std::tuple<class_SVD::MatrixType<Scalar>, class_SVD::VectorType<Scalar>,class_SVD::MatrixType<Scalar> , long>
//do_svd_dgesvd(const Scalar *matrix, int rows, int cols){
//
//}




template<typename Scalar>
std::tuple<class_SVD::MatrixType<Scalar>, class_SVD::VectorType<Scalar>,class_SVD::MatrixType<Scalar> , long>
class_SVD::do_svd_lapacke(const Scalar * mat_ptr, long rows, long cols, long rank_max){
    MatrixType<Scalar> A = Eigen::Map<const MatrixType<Scalar>> (mat_ptr,rows,cols);

    if (rows <= 0)              throw std::runtime_error("SVD error: rows() == 0");
    if (cols <= 0)              throw std::runtime_error("SVD error: cols() == 0");
    if (not A.allFinite())      throw std::runtime_error("SVD error: matrix has inf's or nan's");
    if (A.isZero(0))            throw std::runtime_error("SVD error: matrix is all zeros");

    [[maybe_unused]] int info   = 0;
    int rowsU  = (int) rows;
    int colsU  = (int) std::min(rows,cols);
    int rowsVT = (int) std::min(rows,cols);
    int colsVT = (int) cols;
    int sizeS  = (int) std::min(rows,cols);
    int lda    = (int) rows;
    int ldu    = (int) rowsU;
    int ldvt   = (int) rowsVT;


    MatrixType<Scalar> U(rowsU,colsU);
    VectorType<double> S(sizeS);
    MatrixType<Scalar> VT(rowsVT,colsVT);
    VectorType<Scalar> work(1);

    if constexpr (std::is_same<Scalar,double>::value){
        info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rows,cols, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), -1);
        int lwork  = (int) work(0);
        work.resize(lwork);
        info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rows,cols, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), lwork);
    }
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){
        int lrwork = (int) (5 * std::min(rows,cols));
        VectorType<double> rwork(lrwork);
        auto Ap  =  reinterpret_cast< ComplexDouble_ *>(A.data());
        auto Up  =  reinterpret_cast< ComplexDouble_ *>(U.data());
        auto VTp =  reinterpret_cast< ComplexDouble_ *>(VT.data());
        info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rows,cols, Ap, lda, S.data(), Up, ldu, VTp, ldvt, work.data(), -1,rwork.data());
        int lwork  = (int) std::real(work(0));
        work.resize(lwork);
        info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rows,cols, Ap, lda, S.data(), Up, ldu, VTp, ldvt, work.data(), lwork,rwork.data());
    }

    long max_size    =  std::min(S.size(),rank_max);
    long rank        = (S.head(max_size).array() >= SVDThreshold).count();
    if(rank == S.size()){
        truncation_error = 0;
    }else{
        truncation_error = S.tail(S.size()-rank).squaredNorm();
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
                    << "  V all finite     : " << std::boolalpha << VT.topRows(rank).allFinite() << '\n';
//        return do_svd_lapacke(mat_ptr, rows,cols,rank_max);
        throw std::runtime_error("SVD lapacke error:  Erroneous results");
    }

    return std::make_tuple(
            U.leftCols(rank),
            S.head(rank),
            VT.topRows(rank),
            rank
    );
}



//! \relates class_SVD
//! \brief force instantiation of do_svd_lapacke for type 'double'
template std::tuple<class_SVD::MatrixType<double>, class_SVD::VectorType<double>,class_SVD::MatrixType<double> , long>
class_SVD::do_svd_lapacke(const double *, long, long, long);




using cplx = std::complex<double>;
//! \relates class_SVD
//! \brief force instantiation of do_svd_lapacke for type 'std::complex<double>'
template std::tuple<class_SVD::MatrixType<cplx>, class_SVD::VectorType<cplx>,class_SVD::MatrixType<cplx> , long>
class_SVD::do_svd_lapacke(const cplx *, long, long, long);


