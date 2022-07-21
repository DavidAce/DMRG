#include "../svd.h"
#include "debug/exceptions.h"
#include "tid/tid.h"
#include <complex>

#ifndef lapack_complex_float
    #define lapack_complex_float std::complex<float>
#endif
#ifndef lapack_complex_double
    #define lapack_complex_double std::complex<double>
#endif

// complex must be included before lapacke!
#if __has_include(<mkl_lapacke.h>)
    #include <mkl_lapacke.h>
#elif __has_include(<openblas/lapacke.h>)
    #include <openblas/cblas.h>
    #include <openblas/lapacke.h>
    #include <openblas_config.h>
#else
    #include <lapacke.h>
#endif

namespace svd {
    static constexpr bool use_jsv = true;
    namespace internal {
        // These are workspace arrays used by LAPACK which can be reused for the duration of the program.
        // Call clear() to recover the memory space
        std::vector<int>                  iwork(1ul, 0);
        std::vector<std::complex<double>> cwork(1ul, 0);
        std::vector<double>               rwork(1ul, 0);
        void                              clear_lapack() {
                                         iwork = std::vector<int>(1ul, 0);
                                         cwork = std::vector<std::complex<double>>(1ul, 0);
                                         rwork = std::vector<double>(1ul, 0);
                                         iwork.shrink_to_fit();
                                         cwork.shrink_to_fit();
                                         rwork.shrink_to_fit();
        }
    }
}

template<typename Scalar>
std::tuple<svd::solver::MatrixType<Scalar>, svd::solver::VectorType<Scalar>, svd::solver::MatrixType<Scalar>>
    svd::solver::do_svd_lapacke(const Scalar *mat_ptr, long rows, long cols) const {
    // Setup useful sizes
    int rowsA = static_cast<int>(rows);
    int colsA = static_cast<int>(cols);
    int sizeS = std::min(rowsA, colsA);

    // Setup the SVD solver
    bool use_jacobi = static_cast<size_t>(sizeS) < switchsize_bdc;

    if(use_jacobi and rows < cols) {
        // The jacobi routine needs a tall matrix
        auto t_adj = tid::tic_token("adjoint", tid::detailed);
        svd::log->trace("Transposing {}x{} into a tall matrix {}x{}", rows, cols, cols, rows);
        MatrixType<Scalar> A = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows, cols);
        A.adjointInPlace(); // Adjoint directly on a map seems to give a bug?
        // Sanity checks
        if(A.rows() <= 0) throw std::runtime_error("SVD error: rows() == 0");
        if(A.cols() <= 0) throw std::runtime_error("SVD error: cols() == 0");

        t_adj.toc();
        auto [U, S, VT] = do_svd_lapacke(A.data(), A.rows(), A.cols());
        if(U.rows() != A.rows()) throw except::logic_error("U.rows():{} != A.rows():{}", U.rows(), A.rows());
        if(VT.cols() != A.cols()) throw except::logic_error("VT.cols():{} != A.cols():{}", VT.cols(), A.cols());
        return std::make_tuple(VT.adjoint(), S, U.adjoint());
    }
    auto t_lpk = tid::tic_scope("lapacke", tid::extra);

    // Sanity checks
    if(rows <= 0) throw std::runtime_error("SVD error: rows() <= 0");
    if(cols <= 0) throw std::runtime_error("SVD error: cols() <= 0");

    MatrixType<Scalar>                               A = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows, cols); // gets destroyed in some routines
    MatrixType<Scalar>                               A_original;
    std::vector<std::pair<std::string, std::string>> details;
    if(save_fail or save_result) {
        A_original = A;
#if defined(OPENBLAS_AVAILABLE)
        details = {{"library", "OpenBLAS"},
                   {"OPENBLAS_VERSION", OPENBLAS_VERSION},
                   {"openblas_num_threads", std::to_string(openblas_get_num_threads())},
                   {"openblas_parallel_mode", std::to_string(openblas_get_parallel())},
                   {"openblas_corename", openblas_get_corename()},
                   {"openblas_config", openblas_get_config()},
                   {"OPENBLAS_GEMM_MULTITHREAD_THRESHOLD", std::to_string(OPENBLAS_GEMM_MULTITHREAD_THRESHOLD)}};
#endif
#if defined(MKL_AVAILABLE)
        MKLVersion Version;
        mkl_get_version(&Version);
        details = {{"library", "Intel MKL"}, {"Intel-MKL-Version", fmt::format("{}.{}.{}", Version.MajorVersion, Version.MinorVersion, Version.UpdateVersion)}};
#endif
    }
    // Add suffix for more detailed breakdown of matrix sizes
    auto t_suffix = benchmark ? fmt::format("{}", num::next_multiple<int>(sizeS, 5)) : "";

    // Initialize containers
    MatrixType<Scalar> U;
    VectorType<double> S;
    MatrixType<Scalar> V;
    MatrixType<Scalar> VT;
    svd::log->trace("Starting SVD with lapacke | rows {} | cols {}", rows, cols);

    int info   = 0;
    int rowsU  = rowsA;
    int colsU  = std::min(rowsA, colsA);
    int rowsVT = std::min(rowsA, colsA);
    int colsVT = colsA;
    int rowsV  = colsA;
    int colsV  = std::min(rowsA, colsA);
    int lda    = rowsA;
    int ldu    = rowsU;
    int ldvt   = rowsVT;
    int ldv    = rowsV;

    int mx = std::max(rowsA, colsA);
    int mn = std::min(rowsA, colsA);

    try {
        // More sanity checks
        if(not A.allFinite()) {
            print_matrix(A.data(), A.rows(), A.cols());
            throw std::runtime_error("A has inf's or nan's");
        }
        if(A.isZero(1e-12)) {
            print_matrix(A.data(), A.rows(), A.cols(), 16);
            throw std::runtime_error("A is a zero matrix");
        }
        using namespace svd::internal;
        if constexpr(std::is_same<Scalar, double>::value) {
            if(use_jacobi) {
                if constexpr(use_jsv) {
                    svd::log->debug("Running Lapacke Jacobi SVD | truncation limit {:.4e} | switchsize bdc {} | size {}", truncation_lim, switchsize_bdc,
                                    sizeS);
                    // http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga8767bfcf983f8dc6ef2842029ab25599.html#ga8767bfcf983f8dc6ef2842029ab25599
                    // For this routine we need rows > cols
                    S.resize(sizeS);
                    U.resize(rowsU, colsU);
                    V.resize(rowsV, colsV); // Local matrix gets transposed after computation

                    int lrwork = std::max(2 * rowsA + colsA, 6 * colsA + 2 * colsA * colsA);
                    int liwork = rowsA * 3 * colsA;

                    rwork.resize(static_cast<size_t>(lrwork));
                    iwork.resize(static_cast<size_t>(liwork));

                    svd::log->trace("Running dgesvj");
                    auto t_dgejsv = tid::tic_token(fmt::format("dgejsv{}", t_suffix), tid::detailed);
                    info          = LAPACKE_dgejsv_work(LAPACK_COL_MAJOR, 'F' /* 'R' may also work well */, 'U', 'V', 'N' /* 'R' kills small columns of A */,
                                                        'T' /* T/N:  T will transpose if faster. Ignored if A is rectangular */,
                                                        'N' /* P/N: P will use perturbation to drown denormalized numbers */, rowsA, colsA, A.data(), lda, S.data(),
                                                        U.data(), ldu, V.data(), ldv, rwork.data(), lrwork, iwork.data());
                    t_dgejsv.toc();
                    if(info < 0) throw except::runtime_error("Lapacke SVD dgejsv error: parameter {} is invalid", info);
                    if(info > 0) throw except::runtime_error("Lapacke SVD dgejsv error: could not converge: info {}", info);
                    VT = V.adjoint();
                    V.resize(0, 0);
                } else {
                    svd::log->debug("Running Lapacke Jacobi SVD | truncation limit {:.4e} | switchsize bdc {} | size {}", truncation_lim, switchsize_bdc,
                                    sizeS);
                    // For this routine we need rows > cols
                    S.resize(sizeS);
                    V.resize(rowsV, colsV); // Local matrix gets transposed after computation

                    int lrwork = std::max(6, rowsA + colsA);
                    rwork.resize(static_cast<size_t>(lrwork));

                    svd::log->trace("Running dgesvj");
                    auto t_dgesvj = tid::tic_token(fmt::format("dgesvj{}", t_suffix), tid::detailed);
                    info =
                        LAPACKE_dgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, A.data(), lda, S.data(), ldv, V.data(), ldv, rwork.data(), lrwork);
                    t_dgesvj.toc();
                    if(info < 0) throw except::runtime_error("Lapacke SVD dgesvj error: parameter {} is invalid", info);
                    if(info > 0) throw except::runtime_error("Lapacke SVD dgesvj error: could not converge: info {}", info);
                    U  = A;
                    VT = V.adjoint();
                }
            } else if(use_bdc) {
                svd::log->debug("Running Lapacke BDC SVD | truncation limit {:.4e} | switchsize bdc {} | size {}", truncation_lim, switchsize_bdc, sizeS);
                int liwork = std::max(1, 8 * mn);
                int lrwork = std::max(1, mn * (6 + 4 * mn) + mx);
                iwork.resize(static_cast<size_t>(liwork));
                rwork.resize(static_cast<size_t>(lrwork));

                U.resize(rowsU, colsU);
                S.resize(sizeS);
                VT.resize(rowsVT, colsVT);

                svd::log->trace("Querying dgesdd");
                info = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, rwork.data(), -1,
                                           iwork.data());
                if(info < 0) throw except::runtime_error("Lapacke SVD dgesdd error: parameter {} is invalid", info);
                if(info > 0) throw except::runtime_error("Lapacke SVD dgesdd error: could not converge: info {}", info);

                lrwork = static_cast<int>(rwork[0]);
                rwork.resize(static_cast<size_t>(lrwork));

                svd::log->trace("Running dgesdd");
                auto t_sdd = tid::tic_token(fmt::format("dgesdd{}", t_suffix), tid::detailed);
                info = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, rwork.data(), lrwork,
                                           iwork.data());
                if(info < 0) throw except::runtime_error("Lapacke SVD dgesdd error: parameter {} is invalid", info);
                if(info > 0) throw except::runtime_error("Lapacke SVD dgesdd error: could not converge: info {}", info);

            } else {
                svd::log->debug("Running Lapacke SVD | truncation limit {:.4e} | switchsize bdc {} | size {}", truncation_lim, switchsize_bdc, sizeS);

                U.resize(rowsU, colsU);
                S.resize(sizeS);
                VT.resize(rowsVT, colsVT);

                rwork.resize(1ul);

                svd::log->trace("Querying dgesvd");
                info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, rwork.data(), -1);
                if(info < 0) throw except::runtime_error("Lapacke SVD dgesvd error: parameter {} is invalid", info);
                if(info > 0) throw except::runtime_error("Lapacke SVD dgesvd error: could not converge: info {}", info);

                int lrwork = static_cast<int>(rwork[0]);
                rwork.resize(static_cast<size_t>(lrwork));

                svd::log->trace("Running dgesvd");
                auto t_dgesvd = tid::tic_token(fmt::format("dgesvd{}", t_suffix), tid::detailed);
                info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, rwork.data(),
                                           lrwork);
                if(info < 0) throw except::runtime_error("Lapacke SVD dgesvd error: parameter {} is invalid", info);
                if(info > 0) throw except::runtime_error("Lapacke SVD dgesvd error: could not converge: info {}", info);
            }
        } else if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
            if(use_jacobi) {
                if constexpr(use_jsv) {
                    svd::log->debug("Running Lapacke Jacobi SVD | truncation limit {:.4e} | switchsize bdc {} | size {}", truncation_lim, switchsize_bdc,
                                    sizeS);
                    // http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga8767bfcf983f8dc6ef2842029ab25599.html#ga8767bfcf983f8dc6ef2842029ab25599
                    // For this routine we need rows > cols
                    S.resize(sizeS);
                    U.resize(rowsU, colsU);
                    V.resize(rowsV, colsV); // Local matrix gets transposed after computation

                    int lrwork = std::max(7, 2 * rowsA);
                    int lcwork = std::max(2, lrwork);
                    int liwork = 4; // Minimum size is 4
                    rwork.resize(static_cast<size_t>(lrwork));
                    cwork.resize(static_cast<size_t>(lcwork));
                    iwork.resize(static_cast<size_t>(liwork));

                    svd::log->trace("Querying zgejsv");
                    auto t_zgejsv = tid::tic_token(fmt::format("zgejsv{}", t_suffix), tid::detailed);
                    info          = LAPACKE_zgejsv_work(LAPACK_COL_MAJOR, 'F' /* 'R' may also work well */, 'U', 'V', 'N' /* 'R' kills small columns of A */,
                                                        'T' /* T/N:  T will transpose if faster. Ignored if A is rectangular */,
                                                        'N' /* P/N: P will use perturbation to drown denormalized numbers */, rowsA, colsA, A.data(), lda, S.data(),
                                                        U.data(), ldu, V.data(), ldv, cwork.data(), -1, rwork.data(), -1, iwork.data());

                    if(info < 0) throw except::runtime_error("zgejsv error: parameter {} is invalid", info);
                    if(info > 0) throw except::runtime_error("zgejsv error: could not converge: info {}", info);

                    lcwork = static_cast<int>(std::real(cwork[0]));
                    lrwork = static_cast<int>(rwork[0]);
                    liwork = static_cast<int>(iwork[0]);
                    cwork.resize(static_cast<size_t>(lcwork));
                    rwork.resize(static_cast<size_t>(lrwork));
                    iwork.resize(static_cast<size_t>(liwork));

                    svd::log->trace("Running zgejsv");
                    info = LAPACKE_zgejsv_work(LAPACK_COL_MAJOR, 'F' /* 'R' may also work well */, 'U', 'V', 'N' /* 'R' kills small columns of A */,
                                               'T' /* T/N:  T will transpose if faster. Ignored if A is rectangular */,
                                               'N' /* P/N: P will use perturbation to drown denormalized numbers */, rowsA, colsA, A.data(), lda, S.data(),
                                               U.data(), ldu, V.data(), ldv, cwork.data(), lcwork, rwork.data(), lrwork, iwork.data());

                    if(info < 0) throw except::runtime_error("zgejsv error: parameter {} is invalid", info);
                    if(info > 0) throw except::runtime_error("zgejsv error: could not converge: info {}", info);
                    VT = V.adjoint();
                    V.resize(0, 0);
                } else {
                    svd::log->debug("Running Lapacke Jacobi SVD | truncation limit {:.4e} | switchsize bdc {} | size {}", truncation_lim, switchsize_bdc,
                                    sizeS);
                    S.resize(sizeS);
                    V.resize(rowsV, colsV); // Local matrix gets transposed after computation

                    int lcwork = std::max(6, rowsA + colsA);
                    int lrwork = std::max(6, colsA);
                    cwork.resize(static_cast<size_t>(lcwork));
                    rwork.resize(static_cast<size_t>(lrwork));

                    svd::log->trace("Querying zgesvj");
                    auto t_zgesvj = tid::tic_token(fmt::format("zgesvj{}", t_suffix), tid::detailed);
                    info = LAPACKE_zgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, A.data(), lda, S.data(), ldv, V.data(), ldv, cwork.data(), -1,
                                               rwork.data(), -1);

                    if(info < 0) throw except::runtime_error("zgesvj error: parameter {} is invalid", info);
                    if(info > 0) throw except::runtime_error("zgesvj error: could not converge: info {}", info);

                    lcwork = static_cast<int>(std::real(cwork[0]));
                    lrwork = static_cast<int>(rwork[0]);
                    cwork.resize(static_cast<size_t>(lcwork));
                    rwork.resize(static_cast<size_t>(lrwork));

                    svd::log->trace("Running zgesvj | cwork {} | rwork {}", lcwork, lrwork);
                    info = LAPACKE_zgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, A.data(), lda, S.data(), ldv, V.data(), ldv, cwork.data(), lcwork,
                                               rwork.data(), lrwork);
                    if(info < 0) throw except::runtime_error("zgesvj error: parameter {} is invalid", info);
                    if(info > 0) throw except::runtime_error("zgesvj error: could not converge: info {}", info);
                    U  = A;
                    VT = V.adjoint();
                    V.resize(0, 0);
                }

            } else if(use_bdc) {
                svd::log->debug("Running Lapacke BDC SVD | truncation limit {:.4e} | switchsize bdc {} | size {}", truncation_lim, switchsize_bdc, sizeS);
                U.resize(rowsU, colsU);
                S.resize(sizeS);
                VT.resize(rowsVT, colsVT);

                int lcwork = std::max(1, mn * mn + 3 * mn);
                int lrwork = std::max(1, mn * std::max(5 * mn + 7, 2 * mx + 2 * mn + 1));
                int liwork = std::max(1, 8 * mn);

                cwork.resize(static_cast<size_t>(lcwork));
                rwork.resize(static_cast<size_t>(lrwork));
                iwork.resize(static_cast<size_t>(liwork));

                svd::log->trace("Querying zgesdd");
                auto t_zgesdd = tid::tic_token(fmt::format("zgesdd{}", t_suffix), tid::detailed);
                /* clang-format off */
                info = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, cwork.data(), -1, rwork.data(), iwork.data());
                if(info < 0) throw except::runtime_error("Lapacke SVD zgesdd error: parameter {} is invalid", info);
                if(info > 0) throw except::runtime_error("Lapacke SVD zgesdd error: could not converge: info {}", info);

                lcwork = static_cast<int>(std::real(cwork[0]));
                cwork.resize(static_cast<size_t>(lcwork));

                svd::log->trace("Running zgesdd");
                info = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, cwork.data(), lcwork, rwork.data(), iwork.data());
                if(info < 0) throw except::runtime_error("Lapacke SVD zgesdd error: parameter {} is invalid", info);
                if(info > 0) throw except::runtime_error("Lapacke SVD zgesdd error: could not converge: info {}", info);
                /* clang-format on */

            } else {
                svd::log->debug("Running Lapacke SVD | truncation limit {:.4e} | switchsize bdc {} | size {}", truncation_lim, switchsize_bdc, sizeS);

                U.resize(rowsU, colsU);
                S.resize(sizeS);
                VT.resize(rowsVT, colsVT);

                int lrwork = std::max(1, 5 * mn);
                int lcwork = std::max(1, 2 * mn + mx);
                rwork.resize(static_cast<size_t>(lrwork));
                cwork.resize(static_cast<size_t>(lcwork));

                svd::log->trace("Querying zgesvd");
                auto t_zgesvd = tid::tic_token(fmt::format("zgesvd{}", t_suffix), tid::detailed);
                /* clang-format off */
                info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, cwork.data(), -1, rwork.data());
                if(info < 0) throw except::runtime_error("Lapacke SVD zgesvd error: parameter {} is invalid", info);
                if(info > 0) throw except::runtime_error("Lapacke SVD zgesvd error: could not converge: info {}", info);

                lcwork = static_cast<int>(std::real(cwork[0]));
                cwork.resize(static_cast<size_t>(lcwork));

                svd::log->trace("Running zgesvd");
                info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, cwork.data(), lcwork, rwork.data());
                if(info < 0) throw except::runtime_error("Lapacke SVD zgesvd error: parameter {} is invalid", info);
                if(info > 0) throw except::runtime_error("Lapacke SVD zgesvd error: could not converge: info {}", info);
                /* clang-format on */
            }
        }

        svd::log->trace("Truncating singular values");
        auto max_size                    = S.nonZeros();
        std::tie(rank, truncation_error) = get_rank_from_truncation_error(S.head(max_size).normalized()); // Truncation error needs normalized singular values

        // Do the truncation
        U  = U.leftCols(rank).eval();
        S  = S.head(rank).eval(); // Not all calls to do_svd need normalized S, so we do not normalize here!
        VT = VT.topRows(rank).eval();

        // Sanity checks
        if(not U.allFinite()) {
            print_matrix(U.data(), U.rows(), U.cols());
            throw std::runtime_error("U has inf's or nan's");
        }

        if(not VT.allFinite()) {
            print_matrix(VT.data(), VT.rows(), VT.cols());
            throw std::runtime_error("VT has inf's or nan's");
        }
        if(not S.allFinite()) {
            print_vector(S.data(), rank, 16);
            throw std::runtime_error("S has inf's or nan's");
        }
        if(not(S.array() >= 0).all()) {
            print_vector(S.data(), rank, 16);
            throw std::runtime_error("S is not positive");
        }

    } catch(const std::exception &ex) {
        //#if !defined(NDEBUG)
        if(save_fail) { save_svd<Scalar>(A_original, U, S, VT, "lapacke", details); }
        throw except::runtime_error("Lapacke SVD error \n"
                                    "  Truncation Error = {:.4e}\n"
                                    "  Rank             = {}\n"
                                    "  Dims             = ({}, {})\n"
                                    "  Lapacke info     : {}\n"
                                    "  Error message    : {}\n",
                                    truncation_error, rank, rows, cols, info, ex.what());
    }
    if(save_result) { save_svd<Scalar>(A_original, U, S, VT, "lapacke", details); }

    svd::log->trace(
        "SVD with Lapacke finished successfully | truncation limit {:<8.2e} | rank {:<4} | rank_max {:<4} | {:>4} x {:<4} | trunc {:8.2e}, time {:8.2e}",
        truncation_lim, rank, rank_max, rows, cols, truncation_error, t_lpk->get_last_interval());
    return std::make_tuple(U, S, VT);
}
using real = double;
using cplx = std::complex<double>;

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd_lapacke for type 'double'
template std::tuple<svd::solver::MatrixType<real>, svd::solver::VectorType<real>, svd::solver::MatrixType<real>> svd::solver::do_svd_lapacke(const real *, long,
                                                                                                                                             long) const;

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd_lapacke for type 'std::complex<double>'
template std::tuple<svd::solver::MatrixType<cplx>, svd::solver::VectorType<cplx>, svd::solver::MatrixType<cplx>> svd::solver::do_svd_lapacke(const cplx *, long,
                                                                                                                                             long) const;
