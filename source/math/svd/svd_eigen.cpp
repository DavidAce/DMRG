#include <complex.h>
#undef I

#include "../svd.h"
#include "debug/exceptions.h"
#include "tid/tid.h"
#include <Eigen/QR>
#include <Eigen/SVD>

//#if !defined(NDEBUG)
#include <h5pp/h5pp.h>
//#endif
//

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
    svd::solver::do_svd_eigen(const Scalar *mat_ptr, long rows, long cols, std::optional<long> rank_max) {
    auto t_eigen = tid::tic_scope("eigen", tid::extra);
    if(not rank_max.has_value()) rank_max = std::min(rows, cols);

    svd::log->trace("Starting SVD with Eigen");
    Eigen::Map<const MatrixType<Scalar>> mat(mat_ptr, rows, cols);

    if(rows <= 0) throw except::runtime_error("SVD error: rows = {}", rows);
    if(cols <= 0) throw except::runtime_error("SVD error: cols = {}", cols);

#if !defined(NDEBUG)
    // These are more expensive debugging operations
    if(not mat.allFinite()) throw std::runtime_error("SVD error: matrix has inf's or nan's");
    if(mat.isZero(0)) throw std::runtime_error("SVD error: matrix is all zeros");
    if(mat.isZero(1e-12)) svd::log->warn("Lapacke SVD Warning\n\t Given matrix elements are all close to zero (prec 1e-12)");
#endif
    std::vector<std::pair<std::string, std::string>> details;
    if(save_fail or save_result) details = {{"Eigen Version", fmt::format("{}.{}.{}", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION)}};

    Eigen::BDCSVD<MatrixType<Scalar>> SVD;

    // Setup the SVD solver
    SVD.setSwitchSize(static_cast<int>(switchsize_bdc));
    SVD.setThreshold(threshold);
    // Add suffix for more detailed breakdown of matrix sizes
    auto t_suffix = benchmark ? fmt::format("{}", num::next_multiple<long>(rank_max.value(), 5l)) : "";

    bool use_jacobi = std::min(rows, cols) < static_cast<long>(switchsize_bdc);
    if(use_jacobi) {
        // We only use Jacobi for precision. So we use all the precision we can get.
        svd::log->debug("Running Eigen::JacobiSVD threshold {:.4e} | switchsize bdc {} | rank_max {}", threshold, switchsize_bdc, rank_max.value());
        // Run the svd
        auto t_jcb = tid::tic_token(fmt::format("jcb{}", t_suffix), tid::detailed);
        SVD.compute(mat, Eigen::ComputeFullU | Eigen::ComputeFullV | Eigen::FullPivHouseholderQRPreconditioner);
    } else {
        svd::log->debug("Running Eigen::BDCSVD threshold {:.4e} | switchsize bdc {} | rank_max {}", threshold, switchsize_bdc, rank_max.value());
        // Run the svd
        auto t_bdc = tid::tic_token(fmt::format("bdc{}", t_suffix), tid::detailed);
        SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    }
    if(count) count.value()++;
    long max_size = std::min(SVD.singularValues().size(), rank_max.value());
    long rank     = (SVD.singularValues().head(max_size).array() >= threshold).count();
    svd::log->trace("Truncating singular values");
    if(rank == SVD.singularValues().size()) {
        truncation_error = 0;
    } else {
        truncation_error = SVD.singularValues().tail(SVD.singularValues().size() - rank).norm();
    }

    auto [rank_tr, truncation_error_tr] = truncation_error_limited_rank(SVD.singularValues().head(max_size));
    if(rank_tr < rank) {
        rank             = rank_tr;
        truncation_error = truncation_error_tr;
    }

    bool U_finite   = SVD.matrixU().leftCols(rank).allFinite();
    bool S_finite   = SVD.singularValues().head(rank).allFinite();
    bool V_finite   = SVD.matrixV().leftCols(rank).allFinite();
    bool S_positive = (SVD.singularValues().head(rank).array() >= 0).all();

    if(SVD.rank() <= 0 or rank == 0 or not U_finite or not S_finite or not S_positive or not V_finite) {
        if(not mat.allFinite()) {
            print_matrix(mat.data(), mat.rows(), mat.cols());
            svd::log->critical("Eigen SVD error: matrix has inf's or nan's");
        }
        if(mat.isZero(1e-12)) {
            print_matrix(mat.data(), mat.rows(), mat.cols(), 16);
            svd::log->critical("Eigen SVD error: matrix is all zeros");
        }
        if(not S_positive) {
            print_vector(SVD.singularValues().head(rank).data(), rank, 16);
            svd::log->critical("Eigen SVD error: S is not positive");
        }
        //#if !defined(NDEBUG)
        if(save_fail)
            save_svd<Scalar>(mat, SVD.matrixU().leftCols(rank), SVD.singularValues().head(rank), SVD.matrixV().leftCols(rank).adjoint(), rank_max.value(),
                             "Eigen", details);

        //#endif

        throw except::runtime_error("Eigen SVD error \n"
                                    "  svd_threshold    = {:.4e}\n"
                                    "  Truncation Error = {:.4e}\n"
                                    "  Rank             = {}\n"
                                    "  Rank max         = {}\n"
                                    "  Dims             = ({}, {})\n"
                                    "  A all finite     : {}\n"
                                    "  U all finite     : {}\n"
                                    "  S all finite     : {}\n"
                                    "  S all positive   : {}\n"
                                    "  V all finite     : {}\n",
                                    threshold, truncation_error, rank, rank_max.value(), rows, cols, mat.allFinite(), U_finite, S_finite, S_positive, V_finite);
    }

    if(save_result)
        save_svd<Scalar>(mat, SVD.matrixU().leftCols(rank), SVD.singularValues().head(rank), SVD.matrixV().leftCols(rank).adjoint(), rank_max.value(), "Eigen",
                         details);
    svd::log->trace("SVD with Eigen finished successfully");

    return std::make_tuple(SVD.matrixU().leftCols(rank), SVD.singularValues().head(rank), SVD.matrixV().leftCols(rank).adjoint(), rank);
}

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'double'
template std::tuple<svd::solver::MatrixType<double>, svd::solver::VectorType<double>, svd::solver::MatrixType<double>, long>
    svd::solver::do_svd_eigen(const double *, long, long, std::optional<long>);

using cplx = std::complex<double>;
//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'std::complex<double>'
template std::tuple<svd::solver::MatrixType<cplx>, svd::solver::VectorType<cplx>, svd::solver::MatrixType<cplx>, long>
    svd::solver::do_svd_eigen(const cplx *, long, long, std::optional<long>);