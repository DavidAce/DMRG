#include <complex.h>
#undef I

#include "rsvd/Constants.hpp"
#include "rsvd/ErrorEstimators.hpp"
#include "rsvd/RandomizedSvd.hpp"
#include <Eigen/Dense>
#include <math/svd.h>
#include <tid/tid.h>
#include <math/rnd.h>

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
    svd::solver::do_svd_rsvd(const Scalar *mat_ptr, long rows, long cols, std::optional<long> rank_max) {
    auto t_eigen = tid::tic_scope("rsvd");
    if(not rank_max.has_value()) rank_max = std::min(rows, cols);

    svd::log->trace("Starting SVD with RSVD");
    MatrixType<Scalar> mat = Eigen::Map<const MatrixType<Scalar>> (mat_ptr, rows, cols);

    if(rows <= 0) throw std::runtime_error(fmt::format("SVD error: rows = {}", rows));
    if(cols <= 0) throw std::runtime_error(fmt::format("SVD error: cols = {}", cols));

#if !defined(NDEBUG)
    // These are more expensive debugging operations
    if(not mat.allFinite()) throw std::runtime_error("SVD error: matrix has inf's or nan's");
    if(mat.isZero(0)) throw std::runtime_error("SVD error: matrix is all zeros");
    if(mat.isZero(1e-12)) svd::log->warn("Lapacke SVD Warning\n\t Given matrix elements are all close to zero (prec 1e-12)");
#endif



    // Randomized SVD
//    std::mt19937_64 randomEngine{};
//    randomEngine.seed(777);
//    rnd::

    Rsvd::RandomizedSvd<MatrixType<Scalar>, pcg64, Rsvd::SubspaceIterationConditioner::Lu> SVD(rnd::internal::rng);

    svd::log->debug("Running RSVD threshold {:.4e} | rank_max {}", threshold, rank_max.value());
    // Run the svd
    SVD.compute(mat, rank_max.value(), 10, 4U);

    if(count) count.value()++;
    long max_size = std::min(SVD.singularValues().size(), rank_max.value());
    long rank     = (SVD.singularValues().topRows(max_size).array().real() >= threshold).count();
    svd::log->trace("Truncating singular values");
    if(rank == SVD.singularValues().size()) {
        truncation_error = 0;
    } else {
        truncation_error = SVD.singularValues().bottomRows(SVD.singularValues().size() - rank).norm();
    }

    if(rank == 0 or not SVD.matrixU().leftCols(rank).allFinite() or not SVD.singularValues().topRows(rank).allFinite() or
       not SVD.matrixV().leftCols(rank).allFinite()) {
        throw std::runtime_error(fmt::format(FMT_STRING("RSVD SVD error \n"
                                                        "  svd_threshold    = {:.4e}\n"
                                                        "  Truncation Error = {:.4e}\n"
                                                        "  Rank             = {}\n"
                                                        "  Dims             = ({}, {})\n"
                                                        "  A all finite     : {}\n"
                                                        "  U all finite     : {}\n"
                                                        "  S all finite     : {}\n"
                                                        "  V all finite     : {}\n"),
                                             threshold, truncation_error, rank, rows, cols, mat.allFinite(), SVD.matrixU().leftCols(rank).allFinite(),
                                             SVD.singularValues().topRows(rank).allFinite(), SVD.matrixV().leftCols(rank).allFinite()));
    }
    svd::log->trace("SVD with Eigen finished successfully");

    return std::make_tuple(SVD.matrixU().leftCols(rank), SVD.singularValues().topRows(rank), SVD.matrixV().leftCols(rank).adjoint(), rank);
}

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'double'
template std::tuple<svd::solver::MatrixType<double>, svd::solver::VectorType<double>, svd::solver::MatrixType<double>, long>
    svd::solver::do_svd_rsvd(const double *, long, long, std::optional<long>);

using cplx = std::complex<double>;
//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'std::complex<double>'
template std::tuple<svd::solver::MatrixType<cplx>, svd::solver::VectorType<cplx>, svd::solver::MatrixType<cplx>, long>
    svd::solver::do_svd_rsvd(const cplx *, long, long, std::optional<long>);