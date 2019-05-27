//
// Created by david on 2019-05-27.
//
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
    if (rows <= 0)  {throw std::runtime_error("SVD error: rows() == 0");}
    if (cols <= 0)  {throw std::runtime_error("SVD error: cols() == 0");}
    auto mat = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows,cols);
    Eigen::BDCSVD<MatrixType<Scalar>> SVD;
    SVD.setThreshold(SVDThreshold);
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long rank = std::min(SVD.rank(),rank_max);
    truncation_error = 1.0 - SVD.singularValues().head(rank).squaredNorm();
    if (SVD.rank() <= 0 or mat.isZero(0)){
        std::cout << "SVD error: Rank is zero or invalid matrix" << std::endl;
        std::cerr << "M: \n" << mat << std::endl;
        std::cerr << "U: \n" << SVD.matrixU() << std::endl;
        std::cerr << "S: \n" << SVD.singularValues() << std::endl;
        std::cerr << "V: \n" << SVD.matrixV() << std::endl;
        std::cerr << "rank: " << rank << std::endl;
        throw std::runtime_error("SVD error:  Rank is zero or invalid matrix");
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





