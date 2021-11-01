#include "../svd.h"
#include <complex>
#include <tid/tid.h>

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
    #include <openblas/lapacke.h>
#else
    #include <lapacke.h>
#endif

//#if !defined(NDEBUG)
#include <h5pp/h5pp.h>
#include <math/linalg/matrix.h>
#include <math/linalg/tensor.h>
//#endif
//

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
std::tuple<svd::solver::MatrixType<Scalar>, svd::solver::VectorType<Scalar>, svd::solver::MatrixType<Scalar>, long>
    svd::solver::do_svd_lapacke(const Scalar *mat_ptr, long rows, long cols, std::optional<long> rank_max) {
    // Setup useful sizes
    int rowsA = static_cast<int>(rows);
    int colsA = static_cast<int>(cols);
    int sizeS = std::min(rowsA, colsA);
    if(not rank_max.has_value()) rank_max = std::min(rows, cols);
    if(rank_max.value() == 0) throw std::logic_error("rank_max == 0");
    // Setup the SVD solver
    bool use_jacobi = static_cast<size_t>(sizeS) < switchsize;

    if(use_jacobi and rows < cols) {
        // The jacobi routine needs a tall matrix
        auto t_adj = tid::tic_token("adjoint");
        svd::log->trace("Transposing {}x{} into a tall matrix {}x{}", rows, cols, cols, rows);
        MatrixType<Scalar> A = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows, cols);
        A.adjointInPlace(); // Adjoint directly on a map seems to give a bug?
        // Sanity checks
        if(A.rows() <= 0) throw std::runtime_error("SVD error: rows() == 0");
        if(A.cols() <= 0) throw std::runtime_error("SVD error: cols() == 0");

        t_adj.toc();
        auto [U, S, VT, rank] = do_svd_lapacke(A.data(), A.rows(), A.cols(), std::max(A.rows(), A.cols()));
        long max_size         = std::min(S.size(), rank_max.value());
        rank                  = (S.head(max_size).real().array() >= threshold).count();
        if(U.rows() != A.rows()) throw std::logic_error(fmt::format("U.rows():{} != A.rows():{}", U.rows(), A.rows()));
        if(VT.cols() != A.cols()) throw std::logic_error(fmt::format("VT.cols():{} != A.cols():{}", VT.cols(), A.cols()));
        return std::make_tuple(VT.adjoint().leftCols(rank), S.head(rank), U.adjoint().topRows(rank), rank);
    }
    auto t_lpk = tid::tic_scope("lapacke");

    // Sanity checks
    if(rows <= 0) throw std::runtime_error("SVD error: rows() <= 0");
    if(cols <= 0) throw std::runtime_error("SVD error: cols() <= 0");

    MatrixType<Scalar> A = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows, cols); // gets destroyed in some routines
    MatrixType<Scalar> A_original;
    if(save_fail) A_original = A;

    // Initialize containers
    MatrixType<Scalar> U;
    VectorType<double> S;
    MatrixType<Scalar> V;
    MatrixType<Scalar> VT;
    long               rank;
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
                    svd::log->debug("Running Lapacke Jacobi SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);
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
                    auto t_dgejsv = tid::tic_token("dgejsv");
                    info          = LAPACKE_dgejsv_work(LAPACK_COL_MAJOR, 'F' /* 'R' may also work well */, 'U', 'V', 'N' /* 'R' kills small columns of A */,
                                                        'T' /* T/N:  T will transpose if faster. Ignored if A is rectangular */,
                                                        'N' /* P/N: P will use perturbation to drown denormalized numbers */, rowsA, colsA, A.data(), lda, S.data(),
                                                        U.data(), ldu, V.data(), ldv, rwork.data(), lrwork, iwork.data());
                    t_dgejsv.toc();
                    if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD dgejsv error: parameter {} is invalid", info));
                    if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD dgejsv error: could not converge: info {}", info));
                    VT = V.adjoint();
                    V.resize(0, 0);
                } else {
                    svd::log->debug("Running Lapacke Jacobi SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);
                    // For this routine we need rows > cols
                    S.resize(sizeS);
                    V.resize(rowsV, colsV); // Local matrix gets transposed after computation

                    int lrwork = std::max(6, rowsA + colsA);
                    rwork.resize(static_cast<size_t>(lrwork));

                    svd::log->trace("Running dgesvj");
                    auto t_dgesvj = tid::tic_token("dgesvj");
                    info =
                        LAPACKE_dgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, A.data(), lda, S.data(), ldv, V.data(), ldv, rwork.data(), lrwork);
                    t_dgesvj.toc();
                    if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvj error: parameter {} is invalid", info));
                    if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvj error: could not converge: info {}", info));
                    U  = A;
                    VT = V.adjoint();
                }
            } else if(use_bdc) {
                svd::log->debug("Running Lapacke BDC SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);
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
                if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesdd error: parameter {} is invalid", info));
                if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesdd error: could not converge: info {}", info));

                lrwork = static_cast<int>(rwork[0]);
                rwork.resize(static_cast<size_t>(lrwork));

                svd::log->trace("Running dgesdd");
                auto t_sdd = tid::tic_token("dgesdd");
                info = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, rwork.data(), lrwork,
                                           iwork.data());
                if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesdd error: parameter {} is invalid", info));
                if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesdd error: could not converge: info {}", info));

            } else {
                svd::log->debug("Running Lapacke SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);

                U.resize(rowsU, colsU);
                S.resize(sizeS);
                VT.resize(rowsVT, colsVT);

                rwork.resize(1ul);

                svd::log->trace("Querying dgesvd");
                info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, rwork.data(), -1);
                if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: parameter {} is invalid", info));
                if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: could not converge: info {}", info));

                int lrwork = static_cast<int>(rwork[0]);
                rwork.resize(static_cast<size_t>(lrwork));

                svd::log->trace("Running dgesvd");
                auto t_dgesvd = tid::tic_token("dgesvd");
                info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, rwork.data(),
                                           lrwork);
                if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: parameter {} is invalid", info));
                if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: could not converge: info {}", info));
            }
        } else if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
            if(use_jacobi) {
                if constexpr(use_jsv) {
                    svd::log->debug("Running Lapacke Jacobi SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);
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
                    auto t_zgejsv = tid::tic_token("zgejsv");
                    info          = LAPACKE_zgejsv_work(LAPACK_COL_MAJOR, 'F' /* 'R' may also work well */, 'U', 'V', 'N' /* 'R' kills small columns of A */,
                                                        'T' /* T/N:  T will transpose if faster. Ignored if A is rectangular */,
                                                        'N' /* P/N: P will use perturbation to drown denormalized numbers */, rowsA, colsA, A.data(), lda, S.data(),
                                                        U.data(), ldu, V.data(), ldv, cwork.data(), -1, rwork.data(), -1, iwork.data());

                    if(info < 0) throw std::runtime_error(fmt::format("zgejsv error: parameter {} is invalid", info));
                    if(info > 0) throw std::runtime_error(fmt::format("zgejsv error: could not converge: info {}", info));

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

                    if(info < 0) throw std::runtime_error(fmt::format("zgejsv error: parameter {} is invalid", info));
                    if(info > 0) throw std::runtime_error(fmt::format("zgejsv error: could not converge: info {}", info));
                    VT = V.adjoint();
                    V.resize(0, 0);
                } else {
                    svd::log->debug("Running Lapacke Jacobi SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);
                    S.resize(sizeS);
                    V.resize(rowsV, colsV); // Local matrix gets transposed after computation

                    int lcwork = std::max(6, rowsA + colsA);
                    int lrwork = std::max(6, colsA);
                    cwork.resize(static_cast<size_t>(lcwork));
                    rwork.resize(static_cast<size_t>(lrwork));

                    svd::log->trace("Querying zgesvj");
                    auto t_zgesvj = tid::tic_token("zgesvj");
                    info = LAPACKE_zgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, A.data(), lda, S.data(), ldv, V.data(), ldv, cwork.data(), -1,
                                               rwork.data(), -1);

                    if(info < 0) throw std::runtime_error(fmt::format("zgesvj error: parameter {} is invalid", info));
                    if(info > 0) throw std::runtime_error(fmt::format("zgesvj error: could not converge: info {}", info));

                    lcwork = static_cast<int>(std::real(cwork[0]));
                    lrwork = static_cast<int>(rwork[0]);
                    cwork.resize(static_cast<size_t>(lcwork));
                    rwork.resize(static_cast<size_t>(lrwork));

                    svd::log->trace("Running zgesvj | cwork {} | rwork {}", lcwork, lrwork);
                    info = LAPACKE_zgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, A.data(), lda, S.data(), ldv, V.data(), ldv, cwork.data(), lcwork,
                                               rwork.data(), lrwork);
                    if(info < 0) throw std::runtime_error(fmt::format("zgesvj error: parameter {} is invalid", info));
                    if(info > 0) throw std::runtime_error(fmt::format("zgesvj error: could not converge: info {}", info));
                    U  = A;
                    VT = V.adjoint();
                    V.resize(0, 0);
                }

            } else if(use_bdc) {
                svd::log->debug("Running Lapacke BDC SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);
                U.resize(rowsU, colsU);
                S.resize(sizeS);
                VT.resize(rowsVT, colsVT);

                int liwork = std::max(1, 8 * mn);
                int lrwork = std::max(1, mn * std::max(5 * mn + 7, 2 * mx + 2 * mn + 1));
                int lcwork = std::max(1, mn * mn + 3 * mn);

                iwork.resize(static_cast<size_t>(liwork));
                rwork.resize(static_cast<size_t>(lrwork));
                cwork.resize(static_cast<size_t>(lcwork));

                svd::log->trace("Querying zgesdd");

                auto t_zgesdd = tid::tic_token("zgesdd");
                /* clang-format off */
                info = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, cwork.data(), -1, rwork.data(), iwork.data());
                if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesdd error: parameter {} is invalid", info));
                if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesdd error: could not converge: info {}", info));

                lcwork = static_cast<int>(std::real(cwork[0]));
                cwork.resize(static_cast<size_t>(lcwork));

                svd::log->trace("Running zgesdd");
                info = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, cwork.data(), lcwork, rwork.data(), iwork.data());
                if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesdd error: parameter {} is invalid", info));
                if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesdd error: could not converge: info {}", info));
                /* clang-format on */

            } else {
                svd::log->debug("Running Lapacke SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);

                U.resize(rowsU, colsU);
                S.resize(sizeS);
                VT.resize(rowsVT, colsVT);

                int lrwork = std::max(1, 5 * mn);
                int lcwork = std::max(1, 2 * mn + mx);
                rwork.resize(static_cast<size_t>(lrwork));
                cwork.resize(static_cast<size_t>(lcwork));

                svd::log->trace("Querying zgesvd");
                auto t_zgesvd = tid::tic_token("zgesvd");
                /* clang-format off */
                info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, cwork.data(), -1, rwork.data());
                if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesvd error: parameter {} is invalid", info));
                if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesvd error: could not converge: info {}", info));

                lcwork = static_cast<int>(std::real(cwork[0]));
                cwork.resize(static_cast<size_t>(lcwork));

                svd::log->trace("Running zgesvd");
                info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, cwork.data(), lcwork, rwork.data());
                if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesvd error: parameter {} is invalid", info));
                if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesvd error: could not converge: info {}", info));
                /* clang-format on */
            }
        }

        svd::log->trace("Truncating singular values");
        if(count) count.value()++;
        long max_size = std::min(S.size(), rank_max.value());
        rank          = (S.head(max_size).array() >= threshold).count();
        if(rank == S.size()) {
            truncation_error = 0;
        } else {
            truncation_error = S.tail(S.size() - rank).norm();
        }

        // Do the truncation
        U  = U.leftCols(rank).eval();
        S  = S.head(rank).eval();
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
        if(save_fail) {
            auto file       = h5pp::File("svd-failed.h5", h5pp::FilePermission::READWRITE);
            auto group_num  = 0;
            auto group_name = fmt::format("svd_lapacke_{}", group_num);
            while(file.linkExists(group_name)) group_name = fmt::format("svd_lapacke_{}", ++group_num);
            file.writeDataset(A_original, fmt::format("{}/A", group_name));
            file.writeDataset(U, fmt::format("{}/U", group_name));
            file.writeDataset(S, fmt::format("{}/S", group_name));
            file.writeDataset(VT, fmt::format("{}/V", group_name));
            file.writeAttribute(rows, "rows", group_name);
            file.writeAttribute(cols, "cols", group_name);
            file.writeAttribute(rank, "rank", group_name);
            file.writeAttribute(rank_max.value(), "rank_max", group_name);
            file.writeAttribute(info, "info", group_name);
            file.writeAttribute(use_bdc, "use_bdc", group_name);
            file.writeAttribute(threshold, "threshold", group_name);
            file.writeAttribute(switchsize, "switchsize", group_name);
#if defined(OPENBLAS_AVAILABLE)
            file.writeAttribute(OPENBLAS_VERSION, "OPENBLAS_VERSION", group_name);
            file.writeAttribute(openblas_get_num_threads(), "openblas_get_num_threads", group_name);
            file.writeAttribute(openblas_get_parallel(), "openblas_parallel_mode", group_name);
            file.writeAttribute(openblas_get_corename(), "openblas_get_corename", group_name);
            file.writeAttribute(openblas_get_config(), "openblas_get_config()", group_name);
            file.writeAttribute(OPENBLAS_GEMM_MULTITHREAD_THRESHOLD, "OPENBLAS_GEMM_MULTITHREAD_THRESHOLD", group_name);
#endif

#if defined(MKL_AVAILABLE)
            MKLVersion Version;
            mkl_get_version(&Version);
            file.writeAttribute(Version.MajorVersion, "Intel-MKL-MajorVersion", group_name);
            file.writeAttribute(Version.MinorVersion, "Intel-MKL-MinorVersion", group_name);
            file.writeAttribute(Version.UpdateVersion, "Intel-MKL-UpdateVersion", group_name);
#endif
        }
        throw std::runtime_error(fmt::format(FMT_STRING("Lapacke SVD error \n"
                                                        "  svd_threshold    = {:.4e}\n"
                                                        "  Truncation Error = {:.4e}\n"
                                                        "  Rank             = {}\n"
                                                        "  Dims             = ({}, {})\n"
                                                        "  Lapacke info     : {}\n"
                                                        "  Error message    : {}\n"),
                                             threshold, truncation_error, rank, rows, cols, info, ex.what()));
    }

    svd::log->trace("SVD with lapacke finished successfully. info = {}", info);
    return std::make_tuple(U, S, VT, rank);
}

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd_lapacke for type 'double'
template std::tuple<svd::solver::MatrixType<double>, svd::solver::VectorType<double>, svd::solver::MatrixType<double>, long>
    svd::solver::do_svd_lapacke(const double *, long, long, std::optional<long>);

using cplx = std::complex<double>;
//! \relates svd::class_SVD
//! \brief force instantiation of do_svd_lapacke for type 'std::complex<double>'
template std::tuple<svd::solver::MatrixType<cplx>, svd::solver::VectorType<cplx>, svd::solver::MatrixType<cplx>, long>
    svd::solver::do_svd_lapacke(const cplx *, long, long, std::optional<long>);
