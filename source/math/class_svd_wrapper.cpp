//
// Created by david on 2019-05-27.
//



#ifdef  EIGEN_USE_BLAS
#undef  EIGEN_USE_BLAS
#ifndef EIGEN_USE_LAPACKE_STRICT
#define EIGEN_USE_LAPACKE_STRICT
#endif
#endif

#ifdef  EIGEN_USE_MKL_ALL
#undef  EIGEN_USE_MKL_ALL
#undef  EIGEN_USE_BLAS
#ifndef EIGEN_USE_LAPACKE_STRICT
#define EIGEN_USE_LAPACKE_STRICT
#endif
#endif


#include <complex.h>
#undef I


#include <Eigen/SVD>
#include <math/class_svd_wrapper.h>

double class_SVD::get_truncation_error(){
    return truncation_error;
}

void class_SVD::setThreshold(double newThreshold) {
    SVDThreshold = newThreshold;
}

/*! \brief Performs SVD on a matrix
 *  This function is defined in cpp to avoid long compilation times when having Eigen::BDCSVD included everywhere in headers.
 *  Performs rigorous checks to ensure stability of DMRG.
*   \param mat_ptr Pointer to the matrix. Supported are double * and std::complex<double> *
*   \param rows Rows of the matrix
*   \param cols Columns of the matrix
*   \param rank_max Maximum number of singular values
*   \return The U, S, and V matrices (with S as a vector) extracted from the Eigen::BCDSVD SVD object.
*/
template<typename Scalar>
std::tuple<class_SVD::MatrixType<Scalar>, class_SVD::VectorType<Scalar>,class_SVD::MatrixType<Scalar> , long>
class_SVD::do_svd(const Scalar * mat_ptr, long rows, long cols, long rank_max){
    auto mat = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows,cols);
    if (rows <= 0)              throw std::runtime_error("SVD error: rows() == 0");
    if (cols <= 0)              throw std::runtime_error("SVD error: cols() == 0");
    if (not mat.allFinite())    throw std::runtime_error("SVD error: matrix has inf's or nan's");
    if (mat.isZero(0))          throw std::runtime_error("SVD error: matrix is all zeros");
    Eigen::BDCSVD<MatrixType<Scalar>> SVD;
    SVD.setThreshold(SVDThreshold);
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long rank = std::min(SVD.rank(),rank_max);
    truncation_error = SVD.singularValues().normalized().tail(SVD.nonzeroSingularValues()-rank).squaredNorm();

    if (SVD.rank() <= 0
    or not SVD.matrixU().leftCols(rank).allFinite()
    or not SVD.singularValues().head(rank).allFinite()
    or not SVD.matrixV().leftCols(rank).allFinite() )
    {
        std::cerr   << "SVD error \n"
                    << "  SVDThreshold     = " << SVDThreshold << '\n'
                    << "  Truncation Error = " << truncation_error << '\n'
                    << "  Rank             = " << rank << '\n'
                    << "  U all finite     : " << std::boolalpha << SVD.matrixU().leftCols(rank).allFinite() << '\n'
                    << "  S all finite     : " << std::boolalpha << SVD.singularValues().head(rank).allFinite() << '\n'
                    << "  V all finite     : " << std::boolalpha << SVD.matrixV().leftCols(rank).allFinite() << '\n';
        throw std::runtime_error("SVD error:  Erroneous results");
    }


//    std::cout << "Singular values           : " << SVD.singularValues().transpose() << std::endl;
//    std::cout << "Singular values after norm: " << SVD.singularValues().head(rank).normalized().transpose() << std::endl;
//    std::cout << "Rank                      : " << rank << std::endl;
//    std::cout << "Threshold                 : " << SVDThreshold << std::endl;
//    std::cout << "Truncation error          : " << truncation_error << std::endl;

    return std::make_tuple(
            SVD.matrixU().leftCols(rank),
            SVD.singularValues().head(rank),
            SVD.matrixV().leftCols(rank).adjoint(),
            rank
            );
}

//! \relates class_SVD
//! \brief force instantiation of do_svd for type 'double'
template std::tuple<class_SVD::MatrixType<double>, class_SVD::VectorType<double>,class_SVD::MatrixType<double> , long>
class_SVD::do_svd(const double *, long, long, long);




using cplx = std::complex<double>;
//! \relates class_SVD
//! \brief force instantiation of do_svd for type 'std::complex<double>'
template std::tuple<class_SVD::MatrixType<cplx>, class_SVD::VectorType<cplx>,class_SVD::MatrixType<cplx> , long>
class_SVD::do_svd(const cplx *, long, long, long);



