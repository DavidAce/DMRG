#include <math/svd.h>

#if defined(DMRG_ENABLE_RSVD)
    #include "debug/exceptions.h"
    #include "math/rnd.h"
    #include "rsvd/Constants.hpp"
    #include "rsvd/ErrorEstimators.hpp"
    #include "rsvd/RandomizedSvd.hpp"
    #include "tid/tid.h"
    #include <Eigen/Dense>
#endif

/*! \brief Performs randomized SVD on a matrix
 */
template<typename Scalar>
std::tuple<svd::solver::MatrixType<Scalar>, svd::solver::VectorType<Scalar>, svd::solver::MatrixType<Scalar>, long>
    svd::solver::do_svd_rsvd([[maybe_unused]] const Scalar *mat_ptr, [[maybe_unused]] long rows, [[maybe_unused]] long cols,
                             [[maybe_unused]] std::optional<long> rank_max) {
#if !defined(DMRG_ENABLE_RSVD)
    throw std::runtime_error("Define DMRG_ENABLE_RSVD to use rsvd");
#else
    auto t_rsvd = tid::tic_scope("rsvd");
    if(not rank_max.has_value())
        rank_max = std::min(rows, cols);
    else
        rank_max = std::min({rank_max.value(), rows, cols});

    svd::log->trace("Starting SVD with RSVD | {} x {} |threshold {:.4e} | rank_max {}", rows, cols, threshold, rank_max.value());
    MatrixType<Scalar> mat = Eigen::Map<const MatrixType<Scalar> >(mat_ptr, rows, cols);

    if(rows <= 0) throw except::runtime_error("SVD error: rows = {}", rows);
    if(cols <= 0) throw except::runtime_error("SVD error: cols = {}", cols);

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

    svd::log->debug("Running RSVD | {} x {} |threshold {:.4e} | rank_max {}", rows, cols, threshold, rank_max.value());
    // Run the svd
    //    SVD.compute(mat, rank_max.value(), 10, 4U);
    SVD.compute(mat, rank_max.value());

    if(count) count.value()++;
    long rank     = SVD.rank();
    long max_size = std::min(rank, rank_max.value());
    // Truncation error needs normalized singular values
    std::tie(rank, truncation_error) = get_rank_by_truncation_error(SVD.singularValues().head(max_size).normalized());

    if(rank == 0 or not SVD.matrixU().leftCols(rank).allFinite() or not SVD.singularValues().topRows(rank).allFinite() or
       not SVD.matrixV().leftCols(rank).allFinite()) {
        throw except::runtime_error("RSVD SVD error \n"
                                    "  svd_threshold    = {:.4e}\n"
                                    "  Truncation Error = {:.4e}\n"
                                    "  Rank             = {}\n"
                                    "  Dims             = ({}, {})\n"
                                    "  A all finite     : {}\n"
                                    "  U all finite     : {}\n"
                                    "  S all finite     : {}\n"
                                    "  V all finite     : {}\n",
                                    threshold, truncation_error, rank, rows, cols, mat.allFinite(), SVD.matrixU().leftCols(rank).allFinite(),
                                    SVD.singularValues().topRows(rank).allFinite(), SVD.matrixV().leftCols(rank).allFinite());
    }
    svd::log->trace("SVD with RND SVD finished successfully | threshold {:<8.2e} | rank {:<4} | rank_max {:<4} | {:>4} x {:<4} | trunc {:8.2e}, time {:8.2e}",
                    threshold, rank, rank_max.value(), rows, cols, truncation_error, t_rsvd->get_last_interval());
    // Not all calls to do_svd need normalized S, so we do not normalize here!
    return std::make_tuple(SVD.matrixU().leftCols(rank), SVD.singularValues().topRows(rank), SVD.matrixV().leftCols(rank).adjoint(), rank);
#endif
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