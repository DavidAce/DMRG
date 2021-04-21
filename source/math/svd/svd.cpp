//
// Created by david on 2019-05-27.
//

#include <Eigen/QR>
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
    if(use_lapacke) {
        try {
            return do_svd_lapacke(mat_ptr, rows, cols, rank_max);
        } catch(const std::exception &ex) {
            tools::log->warn("Lapacke failed to perform SVD: {} | Trying Eigen", ex.what());
            return do_svd_eigen(mat_ptr, rows, cols, rank_max);
        }
    }else {
        try {
            return do_svd_eigen(mat_ptr, rows, cols, rank_max);
        } catch(const std::exception &ex) {
            tools::log->warn("Eigen failed to perform SVD: {} | Trying Lapacke", ex.what());
            return do_svd_lapacke(mat_ptr, rows, cols, rank_max);
        }
    }
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