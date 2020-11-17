//
// Created by david on 2019-05-27.
//

#include <complex.h>
#undef I

#include <Eigen/QR>
#include <Eigen/SVD>
#include <math/svd.h>
svd::solver::solver(size_t logLevel) {
    if(not svd::log) { tools::Logger::setLogger(svd::log, "svd", logLevel); }
}

double svd::solver::get_truncation_error() { return truncation_error; }

void svd::solver::setThreshold(double newThreshold, std::optional<double> overrideThreshold) {
    threshold = newThreshold;
    if(overrideThreshold) threshold = overrideThreshold.value();
}
void svd::solver::setSwitchSize(size_t newSwitchSize, std::optional<size_t> overrideSwitchSize) {
    switchsize = newSwitchSize;
    if(overrideSwitchSize) switchsize = overrideSwitchSize.value();
}

/*! \brief Performs SVD on a matrix
 *  This function is defined in cpp to avoid long compilation times when having Eigen::BDCSVD included everywhere in headers.
 *  Performs rigorous checks to ensure stability of DMRG.
 *  In some cases Eigen::BCDSVD/JacobiSVD will fail with segfault. Here we use a patched version of Eigen that throws an error
 *  instead so we get a chance to catch it and use lapack svd instead.
 *   \param mat_ptr Pointer to the matrix. Supported are double * and std::complex<double> *
 *   \param rows Rows of the matrix
 *   \param cols Columns of the matrix
 *   \param rank_max Maximum number of singular values
 *   \return The U, S, and V matrices (with S as a vector) extracted from the Eigen::BCDSVD SVD object.
 */
template<typename Scalar>
std::tuple<svd::solver::MatrixType<Scalar>, svd::solver::VectorType<Scalar>, svd::solver::MatrixType<Scalar>, long>
    svd::solver::do_svd(const Scalar *mat_ptr, long rows, long cols, std::optional<long> rank_max) {
    if(use_lapacke) return do_svd_lapacke(mat_ptr, rows, cols, rank_max);
    if(not rank_max.has_value()) rank_max = std::min(rows, cols);

    svd::log->trace("Starting SVD with Eigen");
    Eigen::Map<const MatrixType<Scalar>> mat(mat_ptr, rows, cols);

    if(rows <= 0) throw std::runtime_error(fmt::format("SVD error: rows = {}", rows));
    if(cols <= 0) throw std::runtime_error(fmt::format("SVD error: cols = {}", cols));

#ifndef NDEBUG
    // These are more expensive debugging operations
    if(not mat.allFinite()) throw std::runtime_error("SVD error: matrix has inf's or nan's");
    if(mat.isZero(0)) throw std::runtime_error("SVD error: matrix is all zeros");
    if(mat.isZero(1e-12)) svd::log->warn("Lapacke SVD Warning\n\t Given matrix elements are all close to zero (prec 1e-12)");
#endif

    Eigen::BDCSVD<MatrixType<Scalar>> SVD;

    // Setup the SVD solver
    if(switchsize) SVD.setSwitchSize(static_cast<int>(switchsize.value()));
    else{
        switchsize = 16;
        SVD.setSwitchSize(static_cast<int>(switchsize.value()));
    }
    if(threshold) SVD.setThreshold(threshold.value());
    else
        threshold = SVD.threshold();
    bool use_jacobi = std::min(rows, cols) < static_cast<long>(switchsize.value());
    unsigned int bitfield = Eigen::ComputeThinU | Eigen::ComputeThinV;
    if(use_jacobi){
        // We only use Jacobi for precision. So we use all the precision we can get.
        bitfield = Eigen::ComputeFullU | Eigen::ComputeFullV | Eigen::FullPivHouseholderQRPreconditioner;
        svd::log->trace("Using JacobiSVD with flags Eigen::ComputeFullU | Eigen::ComputeFullV | Eigen::FullPivHouseholderQRPreconditioner");

    }else{
        svd::log->trace("Using BDCSVD with flags Eigen::ComputeThinU | Eigen::ComputeThinV");
    }

    // Run the svd
    svd::log->trace("Running SVD with threshold {:.4e} | switchsize {} | size {}", threshold.value(), switchsize.value(), rank_max.value());
    SVD.compute(mat, bitfield);
    long max_size = std::min(SVD.singularValues().size(), rank_max.value());
    long rank     = (SVD.singularValues().head(max_size).array() >= threshold.value()).count();
    svd::log->trace("Truncation singular values");
    if(rank == SVD.singularValues().size()) {
        truncation_error = 0;
    } else {
        truncation_error = SVD.singularValues().tail(SVD.singularValues().size() - rank).norm();
    }

    if(SVD.rank() <= 0 or rank == 0 or not SVD.matrixU().leftCols(rank).allFinite() or not SVD.singularValues().head(rank).allFinite() or
       not SVD.matrixV().leftCols(rank).allFinite()) {
        svd::log->warn("Eigen SVD error \n"
                       "  svd_threshold    = {:.4e}\n"
                       "  Truncation Error = {:.4e}\n"
                       "  Rank             = {}\n"
                       "  U all finite     : {}\n"
                       "  S all finite     : {}\n"
                       "  V all finite     : {}\n"
                       "Trying SVD with LAPACKE instead \n",
                       threshold.value(), truncation_error, rank, SVD.matrixU().leftCols(rank).allFinite(), SVD.singularValues().head(rank).allFinite(),
                       SVD.matrixV().leftCols(rank).allFinite());
        return do_svd_lapacke(mat_ptr, rows, cols, rank_max);
    }
    svd::log->trace("SVD with Eigen finished successfully");

    return std::make_tuple(SVD.matrixU().leftCols(rank), SVD.singularValues().head(rank), SVD.matrixV().leftCols(rank).adjoint(), rank);
}

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'double'
template std::tuple<svd::solver::MatrixType<double>, svd::solver::VectorType<double>, svd::solver::MatrixType<double>, long>
    svd::solver::do_svd(const double *, long, long, std::optional<long>);

using cplx = std::complex<double>;
//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'std::complex<double>'
template std::tuple<svd::solver::MatrixType<cplx>, svd::solver::VectorType<cplx>, svd::solver::MatrixType<cplx>, long>
    svd::solver::do_svd(const cplx *, long, long, std::optional<long>);

template<typename Scalar>
Eigen::Tensor<Scalar, 2> svd::solver::pseudo_inverse(const Eigen::Tensor<Scalar, 2> &tensor) {
    if(tensor.dimension(0) <= 0) { throw std::runtime_error("pseudo_inverse error: Dimension is zero: tensor.dimension(0)"); }
    if(tensor.dimension(1) <= 0) { throw std::runtime_error("pseudo_inverse error: Dimension is zero: tensor.dimension(1)"); }
    Eigen::Map<const MatrixType<Scalar>> mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    return Textra::MatrixToTensor(mat.completeOrthogonalDecomposition().pseudoInverse());
}

//! \relates svd::class_SVD
//! \brief force instantiation of pseudo_inverse for type 'double'
template Eigen::Tensor<double, 2> svd::solver::pseudo_inverse(const Eigen::Tensor<double, 2> &tensor);
//! \relates svd::class_SVD
//! \brief force instantiation of pseudo_inverse for type 'std::complex<double>'
template Eigen::Tensor<cplx, 2> svd::solver::pseudo_inverse(const Eigen::Tensor<cplx, 2> &tensor);