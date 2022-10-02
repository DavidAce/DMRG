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
 *   \return The U, S, and V matrices (with S as a vector) extracted from the Eigen::BCDSVD SVD object.
 */
template<typename Scalar>
std::tuple<svd::solver::MatrixType<Scalar>, svd::solver::VectorType<Scalar>, svd::solver::MatrixType<Scalar>>
    svd::solver::do_svd_eigen(const Scalar *mat_ptr, long rows, long cols) const {
    auto t_eigen = tid::tic_scope("eigen", tid::higher);
    svd::log->trace("Starting SVD with Eigen");
    auto                                 minRC = std::min(rows, cols);
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
    // Add suffix for more detailed breakdown of matrix sizes
    auto t_suffix = benchmark ? fmt::format("{}", num::next_multiple<long>(minRC, 5l)) : "";
    auto svd_info = fmt::format("| {} x {} | rank_max {} | truncation limit {:.4e} | switchsize bdc {}", rows, cols, rank_max, truncation_lim, switchsize_bdc);
    bool use_jacobi = minRC < static_cast<long>(switchsize_bdc);
    if(use_jacobi) {
        // We only use Jacobi for precision. So we use all the precision we can get.
        svd::log->debug("Running Eigen::JacobiSVD {}", svd_info);
        // Run the svd
        auto t_jcb = tid::tic_token(fmt::format("jcb{}", t_suffix), tid::highest);
        SVD.compute(mat, Eigen::ComputeFullU | Eigen::ComputeFullV | Eigen::FullPivHouseholderQRPreconditioner);
    } else {
        svd::log->debug("Running Eigen::BDCSVD {}", svd_info);
        // Run the svd
        auto t_bdc = tid::tic_token(fmt::format("bdc{}", t_suffix), tid::highest);
        SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    }

    long rank_lim   = rank_max > 0 ? std::min(minRC, rank_max) : minRC;
    rank            = SVD.nonzeroSingularValues();
    long max_size   = std::min(rank, rank_lim);
    bool U_finite   = SVD.matrixU().leftCols(max_size).allFinite();
    bool S_finite   = SVD.singularValues().head(max_size).allFinite();
    bool V_finite   = SVD.matrixV().leftCols(max_size).allFinite();
    bool S_positive = (SVD.singularValues().head(max_size).array() >= 0).all();

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
        if(save_fail) save_svd<Scalar>(mat, SVD.matrixU(), SVD.singularValues(), SVD.matrixV().adjoint(), "Eigen", details);

        //#endif

        throw except::runtime_error("Eigen SVD error \n"
                                    "  Rank max       = {}\n"
                                    "  Dims             = ({}, {})\n"
                                    "  A all finite     : {}\n"
                                    "  U all finite     : {}\n"
                                    "  S all finite     : {}\n"
                                    "  S all positive   : {}\n"
                                    "  V all finite     : {}\n",
                                    rank_max, rows, cols, mat.allFinite(), U_finite, S_finite, S_positive, V_finite);
    }

    // Truncation error needs normalized singular values
    std::tie(rank, truncation_error) = get_rank_from_truncation_error(SVD.singularValues().head(max_size).normalized());

    if(save_result)
        save_svd<Scalar>(mat, SVD.matrixU().leftCols(rank), SVD.singularValues().head(rank), SVD.matrixV().leftCols(rank).adjoint(), "Eigen", details);
    svd::log->trace("SVD with Eigen finished successfully");
    // Not all calls to do_svd need normalized S, so we do not normalize here!
    return std::make_tuple(SVD.matrixU().leftCols(rank), SVD.singularValues().head(rank), SVD.matrixV().leftCols(rank).adjoint());
}

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'double'
template std::tuple<svd::solver::MatrixType<double>, svd::solver::VectorType<double>, svd::solver::MatrixType<double>>
    svd::solver::do_svd_eigen(const double *, long, long) const;

using cplx = std::complex<double>;
//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'std::complex<double>'
template std::tuple<svd::solver::MatrixType<cplx>, svd::solver::VectorType<cplx>, svd::solver::MatrixType<cplx>> svd::solver::do_svd_eigen(const cplx *, long,
                                                                                                                                           long) const;