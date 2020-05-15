//
// Created by david on 2019-05-27.
//



//#define EIGEN_USE_MKL_ALL
//#define EIGEN_USE_BLAS
//#define EIGEN_USE_LAPACKE_STRICT

//#define EIGEN_DONT_VECTORIZE
//#define EIGEN_DONT_ALIGN_STATICALLY

//#define EIGEN_MAX_ALIGN_BYTES 0
//#define EIGEN_ENABLE_AVX512
//#define EIGEN_UNALIGNED_VECTORIZE 0
//#define EIGEN_DONT_ALIGN
//#define EIGEN_MALLOC_ALREADY_ALIGNED 0


#include <complex.h>
#undef I

//For svd debugging
//#include <h5pp/h5pp.h>
//#include <config/nmspc_settings.h>


#include <Eigen/SVD>
#include <Eigen/QR>
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
class_SVD::do_svd(const Scalar * mat_ptr, long rows, long cols, std::optional<long> rank_max){
    if (use_lapacke) return do_svd_lapacke(mat_ptr, rows,cols,rank_max);
    if(not rank_max.has_value()) rank_max = std::min(rows,cols);

    MatrixType<Scalar> mat = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows,cols);
    if (rows <= 0)              throw std::runtime_error("SVD error: rows() == 0");
    if (cols <= 0)              throw std::runtime_error("SVD error: cols() == 0");
    if (not mat.allFinite())    throw std::runtime_error("SVD error: matrix has inf's or nan's");
    if (mat.isZero(0))          throw std::runtime_error("SVD error: matrix is all zeros");

    Eigen::BDCSVD<MatrixType<Scalar>> SVD;
    SVD.setThreshold(SVDThreshold);
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long max_size =  std::min(SVD.singularValues().size(),rank_max.value());
    long rank     = (SVD.singularValues().head(max_size).array() >= SVDThreshold).count();
    if(rank == SVD.singularValues().size()){
        truncation_error = 0;
    }else{
        truncation_error = SVD.singularValues().tail(SVD.singularValues().size() - rank).norm();
    }

    if (SVD.rank() <= 0
    or not SVD.matrixU().leftCols(rank).allFinite()
    or not SVD.singularValues().head(rank).allFinite()
    or not SVD.matrixV().leftCols(rank).allFinite() )
    {
        std::cerr   << "SVD error \n"
                    << "  svd_threshold     = " << SVDThreshold << '\n'
                    << "  Truncation Error = " << truncation_error << '\n'
                    << "  Rank             = " << rank << '\n'
                    << "  U all finite     : " << std::boolalpha << SVD.matrixU().leftCols(rank).allFinite() << '\n'
                    << "  S all finite     : " << std::boolalpha << SVD.singularValues().head(rank).allFinite() << '\n'
                    << "  V all finite     : " << std::boolalpha << SVD.matrixV().leftCols(rank).allFinite() << '\n'
                    << "Trying SVD with LAPACKE instead \n";
        return do_svd_lapacke(mat_ptr, rows,cols,rank_max);
    }

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
class_SVD::do_svd(const double *, long, long, std::optional<long>);




using cplx = std::complex<double>;
//! \relates class_SVD
//! \brief force instantiation of do_svd for type 'std::complex<double>'
template std::tuple<class_SVD::MatrixType<cplx>, class_SVD::VectorType<cplx>,class_SVD::MatrixType<cplx> , long>
class_SVD::do_svd(const cplx *, long, long, std::optional<long>);





template<typename Scalar>
Eigen::Tensor<Scalar, 2>
class_SVD::pseudo_inverse(const Eigen::Tensor<Scalar, 2> &tensor){
    if (tensor.dimension(0) <= 0)  {throw std::runtime_error("pseudo_inverse error: Dimension is zero: tensor.dimension(0)");}
    if (tensor.dimension(1) <= 0)  {throw std::runtime_error("pseudo_inverse error: Dimension is zero: tensor.dimension(1)");}
    Eigen::Map<const MatrixType<Scalar>> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    return Textra::MatrixTensorMap(mat.completeOrthogonalDecomposition().pseudoInverse() );
}



//! \relates class_SVD
//! \brief force instantiation of pseudo_inverse for type 'double'
template Eigen::Tensor<double, 2> class_SVD::pseudo_inverse(const Eigen::Tensor<double, 2> &tensor);
//! \relates class_SVD
//! \brief force instantiation of pseudo_inverse for type 'std::complex<double>'
template Eigen::Tensor<cplx, 2>   class_SVD::pseudo_inverse(const Eigen::Tensor<cplx  , 2> &tensor);