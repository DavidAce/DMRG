//
// Created by david on 2019-05-27.
//


#ifdef EIGEN_USE_BLAS
#undef EIGEN_USE_BLAS
#endif

#include <Eigen/SVD>
#include <general/class_svd_wrapper.h>

double class_SVD::get_truncation_error(){
    return truncation_error;
}

void class_SVD::setThreshold(double newThreshold) {
    SVDThreshold = newThreshold;
}


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
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-rank).squaredNorm();

    if (SVD.rank() <= 0 or not SVD.matrixU().allFinite() or not SVD.matrixV().allFinite() ){
        std::cerr << "M: \n" << mat << std::endl;
        std::cerr << "U: \n" << SVD.matrixU() << std::endl;
        std::cerr << "S: \n" << SVD.singularValues() << std::endl;
        std::cerr << "V: \n" << SVD.matrixV() << std::endl;
        std::cerr << "rank: " << rank << std::endl;
        throw std::runtime_error("SVD error:  Erroneous results");
    }

    return std::make_tuple(
            SVD.matrixU().leftCols(rank),
            SVD.singularValues().head(rank),
            SVD.matrixV().leftCols(rank).adjoint(),
            rank
            );
}

template std::tuple<class_SVD::MatrixType<double>, class_SVD::VectorType<double>,class_SVD::MatrixType<double> , long>
class_SVD::do_svd(const double *, long, long, long);

using cplx = std::complex<double>;
template std::tuple<class_SVD::MatrixType<cplx>, class_SVD::VectorType<cplx>,class_SVD::MatrixType<cplx> , long>
class_SVD::do_svd(const cplx *, long, long, long);



