//
// Created by david on 2019-05-27.
//

#include <complex.h>
#undef I

#include <Eigen/QR>
#include <Eigen/SVD>
#include <general/class_tic_toc.h>
#include <iostream>
#include <math/svd.h>

std::optional<long long> svd::solver::count = 0;

svd::solver::solver() {
    setLogLevel(2);
    t_wrk = std::make_unique<class_tic_toc>();
    t_adj = std::make_unique<class_tic_toc>();
    t_jac = std::make_unique<class_tic_toc>();
    t_svd = std::make_unique<class_tic_toc>();
    if(not count) count = 0;
}
void svd::solver::copy_settings(const svd::settings &svd_settings) {
    if(svd_settings.threshold) threshold = svd_settings.threshold.value();
    if(svd_settings.switchsize) switchsize = svd_settings.switchsize.value();
    if(svd_settings.loglevel) setLogLevel(svd_settings.loglevel.value());
    if(svd_settings.use_bdc) use_bdc = svd_settings.use_bdc.value();
    if(svd_settings.use_lapacke) use_lapacke = svd_settings.use_lapacke.value();
    if(svd_settings.profile and svd_settings.profile.value()) enableProfiling();
}

svd::solver::solver(const svd::settings &svd_settings) : solver() { copy_settings(svd_settings); }

svd::solver::solver(std::optional<svd::settings> svd_settings) : solver() {
    if(svd_settings) copy_settings(svd_settings.value());
}

void svd::solver::enableProfiling() {
    t_wrk->set_properties(true, 5, "work");
    t_adj->set_properties(true, 5, "adjoint");
    t_jac->set_properties(true, 5, "jacobi");
    t_svd->set_properties(true, 5, "bdcsvd");
}
void svd::solver::disableProfiling() {
    t_wrk->set_properties(false, 0, "");
    t_adj->set_properties(false, 0, "");
    t_jac->set_properties(false, 0, "");
    t_svd->set_properties(false, 0, "");
}

void svd::solver::setLogLevel(size_t logLevel) {
    if(not svd::log)
        tools::Logger::setLogger(svd::log, "svd", logLevel, true);
    else
        tools::Logger::setLogLevel(svd::log, logLevel);
    svd::log->set_pattern("[%Y-%m-%d %H:%M:%S.%e][%n]%^[%=8l]%$ %v");
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
    SVD.setSwitchSize(static_cast<int>(switchsize));
    SVD.setThreshold(threshold);
    bool use_jacobi = std::min(rows, cols) < static_cast<long>(switchsize);
    svd::log->trace("Running SVD with threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, rank_max.value());
    if(use_jacobi) {
        // We only use Jacobi for precision. So we use all the precision we can get.
        svd::log->debug("Running Eigen::JacobiSVD threshold {:.4e} | switchsize {} | rank_max {}", threshold, switchsize, rank_max.value());
        // Run the svd
        t_jac->tic();
        SVD.compute(mat, Eigen::ComputeFullU | Eigen::ComputeFullV | Eigen::FullPivHouseholderQRPreconditioner);
        t_jac->toc();
    } else {
        svd::log->debug("Running Eigen::BDCSVD threshold {:.4e} | switchsize {} | rank_max {}", threshold, switchsize, rank_max.value());
        // Run the svd
        t_svd->tic();
        SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
        t_svd->toc();
    }
    if(count) count.value()++;
    long max_size = std::min(SVD.singularValues().size(), rank_max.value());
    long rank     = (SVD.singularValues().head(max_size).array() >= threshold).count();
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
                       threshold, truncation_error, rank, SVD.matrixU().leftCols(rank).allFinite(), SVD.singularValues().head(rank).allFinite(),
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
    return Textra::TensorCast(mat.completeOrthogonalDecomposition().pseudoInverse());
}

//! \relates svd::class_SVD
//! \brief force instantiation of pseudo_inverse for type 'double'
template Eigen::Tensor<double, 2> svd::solver::pseudo_inverse(const Eigen::Tensor<double, 2> &tensor);
//! \relates svd::class_SVD
//! \brief force instantiation of pseudo_inverse for type 'std::complex<double>'
template Eigen::Tensor<cplx, 2> svd::solver::pseudo_inverse(const Eigen::Tensor<cplx, 2> &tensor);