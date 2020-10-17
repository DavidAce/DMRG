//
// Created by david on 2019-08-07.
//

#include <complex.h>
#undef I
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
    #include <openblas/lapacke.h>
#else
    #include <lapacke.h>
#endif

#include <Eigen/Core>
#include <math/svd.h>
template<typename Scalar>
std::tuple<svd::solver::MatrixType<Scalar>, svd::solver::VectorType<Scalar>, svd::solver::MatrixType<Scalar>, long>
    svd::solver::do_svd_lapacke(const Scalar *mat_ptr, long rows, long cols, std::optional<long> rank_max) {
    int rowsA = static_cast<int>(rows);
    int colsA = static_cast<int>(cols);
    int sizeS = std::min(rowsA, colsA);

    // Setup the SVD solver
    double lapacke_svd_threshold  = static_cast<double>(sizeS) * std::numeric_limits<double>::epsilon(); // Same as Eigen
    size_t lapacke_svd_switchsize = 16ul;                                                                // Same as Eigen
    if(threshold) lapacke_svd_threshold = threshold.value();
    if(switchsize) lapacke_svd_switchsize = switchsize.value();
    bool use_jacobi = static_cast<size_t>(sizeS) <= lapacke_svd_switchsize;
    if(use_jacobi and rows < cols) {
        // The jacobi routine needs a tall matrix
        svd::log->trace("Transposing into tall matrix");
        MatrixType<Scalar> A = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows, cols).adjoint();
        auto [U, S, V, rank] = do_svd_lapacke(A.data(), A.rows(), A.cols(), std::max(A.rows(), A.cols()));
        long max_size        = std::min(S.size(), std::min(rows, cols));
        rank                 = (S.head(max_size).real().array() >= lapacke_svd_threshold).count();
        return std::make_tuple(V.adjoint().leftCols(rank), S.head(rank), U.adjoint().topRows(rank), rank);
    }

    svd::log->trace("Starting SVD with lapacke");

    if(not rank_max.has_value()) rank_max = std::min(rows, cols);
    if(rows <= 0) throw std::runtime_error("SVD error: rows() == 0");
    if(cols <= 0) throw std::runtime_error("SVD error: cols() == 0");

    MatrixType<Scalar> A = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows, cols);
    if(not A.allFinite()) throw std::runtime_error("SVD error: matrix has inf's or nan's");
    if(A.isZero(1e-12)) throw std::runtime_error("SVD error: matrix is all zeros");
#ifndef NDEBUG
    // These are more expensive debugging operations
    if(not A.allFinite()) throw std::runtime_error("SVD error: matrix has inf's or nan's");
    if(A.isZero(0)) throw std::runtime_error("SVD error: matrix is all zeros");
    if(A.isZero(1e-12)) svd::log->warn("Lapacke SVD Warning\n\t Given matrix elements are all close to zero (prec 1e-12)");
#endif

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

    MatrixType<Scalar> U;
    VectorType<double> S;
    MatrixType<Scalar> VT;

    if constexpr(std::is_same<Scalar, double>::value) {
        if(use_jacobi) {
            svd::log->trace("Running Jacobi SVD with threshold {:.4e} | switchsize {} | size {}", lapacke_svd_threshold, lapacke_svd_switchsize, sizeS);
            // http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga8767bfcf983f8dc6ef2842029ab25599.html#ga8767bfcf983f8dc6ef2842029ab25599
            // For this routine we need rows > cols
            svd::log->trace("Resizing S and V containers");
            S.resize(sizeS);
            MatrixType<Scalar> V(rowsV, colsV); // Local matrix gets transposed after computation
            svd::log->trace("Resizing work arrays");
            int                lwork = std::max(6, rowsA + colsA);
            VectorType<Scalar> work(lwork);
            work.setConstant(0);
            work(0) = 1;
            svd::log->trace("Running dgejsv");
            info = LAPACKE_dgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, A.data(), lda, S.data(), ldv, V.data(), ldv, work.data(), lwork);
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));
            long max_size = std::min(S.size(), rank_max.value());
            long rank     = (S.head(max_size).array() >= lapacke_svd_threshold).count();
            U             = A.leftCols(rank);
            VT            = V.adjoint().topRows(rank);
        } else {
            svd::log->trace("Running BDC SVD with threshold {:.4e} | switchsize {} | size {}", lapacke_svd_threshold, lapacke_svd_switchsize, sizeS);
            svd::log->trace("Resizing USV^T containers");
            U.resize(rowsU, colsU);
            S.resize(sizeS);
            VT.resize(rowsVT, colsVT);
            VectorType<Scalar> work(1);
            svd::log->trace("Querying dgesvd");
            info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), -1);
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));
            svd::log->trace("Resizing work arrays");
            int lwork = static_cast<int>(work(0));
            work.resize(lwork);
            svd::log->trace("Running dgesvd");
            info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), lwork);
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));
        }
    }
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
        if(use_jacobi) {
            svd::log->trace("Running Jacobi SVD with threshold {:.4e} | switchsize {} | size {}", lapacke_svd_threshold, lapacke_svd_switchsize, sizeS);
            svd::log->trace("Resizing S and V containers");
            S.resize(sizeS);
            MatrixType<Scalar> V(rowsV, colsV); // Local matrix gets transposed after computation
            VectorType<Scalar> work(1);
            VectorType<double> rwork(1);
            auto               Ap = reinterpret_cast<lapack_complex_double *>(A.data());
            auto               Vp = reinterpret_cast<lapack_complex_double *>(V.data());
            auto               Wp = reinterpret_cast<lapack_complex_double *>(work.data());

            svd::log->trace("Querying zgesvj");
            info = LAPACKE_zgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, Ap, lda, S.data(), ldv, Vp, ldv, Wp, -1, rwork.data(), -1);
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));

            svd::log->trace("Resizing work array");
            int lwork  = static_cast<int>(std::real(work(0)));
            int lrwork = static_cast<int>(rwork(0));
            work.resize(lwork);
            rwork.resize(lrwork);
            Wp = reinterpret_cast<lapack_complex_double *>(work.data());

            svd::log->trace("Running zgesvj");
            info = LAPACKE_zgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, Ap, lda, S.data(), ldv, Vp, ldv, Wp, lwork, rwork.data(), lrwork);

            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));
            long max_size = std::min(S.size(), rank_max.value());
            long rank     = (S.head(max_size).array() >= lapacke_svd_threshold).count();
            U             = A.leftCols(rank);
            VT            = V.adjoint().topRows(rank);

        } else {
            svd::log->trace("Running BDC SVD with threshold {:.4e} | switchsize {} | size {}", lapacke_svd_threshold, lapacke_svd_switchsize, sizeS);
            svd::log->trace("Resizing USV^T containers");
            int lrwork = 5 * std::min(rowsA, colsA);
            U.resize(rowsU, colsU);
            S.resize(sizeS);
            VT.resize(rowsVT, colsVT);
            VectorType<Scalar> work(1);
            VectorType<double> rwork(lrwork);
            auto               Ap  = reinterpret_cast<lapack_complex_double *>(A.data());
            auto               Up  = reinterpret_cast<lapack_complex_double *>(U.data());
            auto               VTp = reinterpret_cast<lapack_complex_double *>(VT.data());
            auto               Wp  = reinterpret_cast<lapack_complex_double *>(work.data());

            svd::log->trace("Querying zgesvd");
            info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, Ap, lda, S.data(), Up, ldu, VTp, ldvt, Wp, -1, rwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));

            svd::log->trace("Resizing work arrays");
            int lwork = static_cast<int>(std::real(work(0)));
            work.resize(lwork);
            Wp = reinterpret_cast<lapack_complex_double *>(work.data()); // Update the pointer if reallocated

            svd::log->trace("Running zgesvd");
            info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, Ap, lda, S.data(), Up, ldu, VTp, ldvt, Wp, lwork, rwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));
        }
    }
    svd::log->trace("Truncation singular values");

    long max_size = std::min(S.size(), rank_max.value());
    long rank     = (S.head(max_size).array() >= lapacke_svd_threshold).count();
    if(rank == S.size()) {
        truncation_error = 0;
    } else {
        truncation_error = S.tail(S.size() - rank).norm();
    }

    if(rank <= 0 or not U.leftCols(rank).allFinite() or not S.head(rank).allFinite() or not VT.topRows(rank).allFinite()) {
        svd::log->warn("Lapacke SVD error \n"
                       "  svd_threshold    = {:.4e}\n"
                       "  Truncation Error = {:.4e}\n"
                       "  Rank             = {}\n"
                       "  U all finite     : {}\n"
                       "  S all finite     : {}\n"
                       "  V all finite     : {}\n"
                       "  Lapacke info     : {}\n",
                       lapacke_svd_threshold, truncation_error, rank, U.leftCols(rank).allFinite(), S.head(rank).allFinite(), VT.topRows(rank).allFinite(),
                       info);
        if(not use_lapacke) throw std::runtime_error("Lapacke SVD error:  Wrong results");
    }
    svd::log->trace("SVD with lapacke finished successfully. info = {}", info);

    return std::make_tuple(U.leftCols(rank), S.head(rank), VT.topRows(rank), rank);
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
