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

namespace svd {
    template<typename Scalar>
    void print_matrix_lapacke(const Scalar *mat_ptr, long rows, long cols, long dec = 8) {
        auto A = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(mat_ptr, rows, cols);
        svd::log->warn("Error printing matrix of dimensions {}x{}\n", rows, cols);
        for(long r = 0; r < A.rows(); r++) {
            if constexpr(std::is_same_v<Scalar, std::complex<double>>)
                for(long c = 0; c < A.cols(); c++) fmt::print("({1:.{0}f},{2:+.{0}f}) ", dec, std::real(A(r, c)), std::imag(A(r, c)));
            else
                for(long c = 0; c < A.cols(); c++) fmt::print("{1:.{0}f} ", dec, A(r, c));
            fmt::print("\n");
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

    // Setup the SVD solver
    double lapacke_svd_threshold  = static_cast<double>(sizeS) * std::numeric_limits<double>::epsilon();
    size_t lapacke_svd_switchsize = 16ul; // Same as Eigen
    if(threshold) lapacke_svd_threshold = threshold.value();
    if(switchsize) lapacke_svd_switchsize = switchsize.value();
    bool use_jacobi = static_cast<size_t>(sizeS) < lapacke_svd_switchsize;
    if(use_jacobi and rows < cols) {
        // The jacobi routine needs a tall matrix
        svd::log->trace("Transposing {}x{} into tall matrix {}x{}", rows, cols, cols, rows);
        MatrixType<Scalar> A = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows, cols);
        A.adjointInPlace(); // Adjoint directly on a map seems to give a bug?
        svd::log->trace("Sanity checks... ");
        // Sanity checks
        if(A.rows() <= 0) throw std::runtime_error("SVD error: rows() == 0");
        if(A.cols() <= 0) throw std::runtime_error("SVD error: cols() == 0");
        if(not A.allFinite()) {
            print_matrix_lapacke(A.data(), A.rows(), A.cols());
            throw std::runtime_error("SVD error: matrix has inf's or nan's");
        }
        if(A.isZero(1e-12)) {
            print_matrix_lapacke(A.data(), A.rows(), A.cols(), 16);
            throw std::runtime_error("SVD error: matrix is all zeros");
        }
        svd::log->trace("Sanity checks... OK. A: {}x{}", A.rows(), A.cols());
        auto [U, S, VT, rank] = do_svd_lapacke(A.data(), A.rows(), A.cols(), std::max(A.rows(), A.cols()));
        long max_size         = std::min(S.size(), rank_max.value());
        rank                  = (S.head(max_size).real().array() >= lapacke_svd_threshold).count();
        if(U.rows() != A.rows()) throw std::logic_error(fmt::format("U.rows():{} != A.rows():{}", U.rows(), A.rows()));
        if(VT.cols() != A.cols()) throw std::logic_error(fmt::format("VT.cols():{} != A.cols():{}", VT.cols(), A.cols()));
        return std::make_tuple(VT.adjoint().leftCols(rank), S.head(rank), U.adjoint().topRows(rank), rank);
    }

    // Sanity checks
    if(rows <= 0) throw std::runtime_error("SVD error: rows() == 0");
    if(cols <= 0) throw std::runtime_error("SVD error: cols() == 0");

    MatrixType<Scalar> A = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows, cols);
    if(not A.allFinite()) {
        print_matrix_lapacke(mat_ptr, rows, cols);
        throw std::runtime_error("SVD error: matrix has inf's or nan's");
    }
    if(A.isZero(1e-12)) {
        print_matrix_lapacke(mat_ptr, rows, cols, 16);
        throw std::runtime_error("SVD error: matrix is all zeros");
    }

    svd::log->trace("Starting SVD with lapacke");

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
            int lwork = std::max(6, rowsA + colsA);
            std::vector<Scalar> work(lwork);
//            work.setConstant(0);
//            work[0] = 1;
            S.resize(sizeS);
            MatrixType<Scalar> V(rowsV, colsV); // Local matrix gets transposed after computation

            svd::log->trace("Running dgejsv");
            info = LAPACKE_dgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, A.data(), lda, S.data(), ldv, V.data(), ldv, work.data(), lwork);
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));
            long max_size = std::min(S.size(), rank_max.value());
            long rank     = (S.head(max_size).array() >= lapacke_svd_threshold).count();
            U             = A.leftCols(rank);
            VT            = V.adjoint().topRows(rank);
        }
        else if (use_bdc) {
            svd::log->trace("Running BDC SVD with threshold {:.4e} | switchsize {} | size {}", lapacke_svd_threshold, lapacke_svd_switchsize, sizeS);
            int liwork = std::max(1, 8 * std::min(rowsA, colsA));
            std::vector<Scalar> work(1);
            std::vector<int>    iwork(liwork);

            U.resize(rowsU, colsU);
            S.resize(sizeS);
            VT.resize(rowsVT, colsVT);

            svd::log->trace("Querying dgesvd");
            info = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), -1, iwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));

            int lwork = static_cast<int>(work[0]);
            work.resize(lwork);

            svd::log->trace("Running dgesvd");
            info = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), lwork, iwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));
        }
        else {
            svd::log->trace("Running BDC SVD with threshold {:.4e} | switchsize {} | size {}", lapacke_svd_threshold, lapacke_svd_switchsize, sizeS);
            std::vector<Scalar> work(1);

            U.resize(rowsU, colsU);
            S.resize(sizeS);
            VT.resize(rowsVT, colsVT);

            svd::log->trace("Querying dgesvd");
            info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), -1);
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));

            int lwork = static_cast<int>(work[0]);
            work.resize(lwork);

            svd::log->trace("Running dgesvd");
            info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), lwork);
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));
        }
    }
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
        if(use_jacobi) {
            svd::log->trace("Running Jacobi SVD with threshold {:.4e} | switchsize {} | size {}", lapacke_svd_threshold, lapacke_svd_switchsize, sizeS);
            std::vector<Scalar> work(1);
            std::vector<double> rwork(1);

            S.resize(sizeS);
            MatrixType<Scalar> V(rowsV, colsV); // Local matrix gets transposed after computation

            auto               Ap = reinterpret_cast<lapack_complex_double *>(A.data());
            auto               Vp = reinterpret_cast<lapack_complex_double *>(V.data());
            auto               Wp = reinterpret_cast<lapack_complex_double *>(work.data());

            svd::log->trace("Querying zgesvj");
            info = LAPACKE_zgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, Ap, lda, S.data(), ldv, Vp, ldv, Wp, -1, rwork.data(), -1);
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));

            int lwork  = static_cast<int>(std::real(work[0]));
            int lrwork = static_cast<int>(rwork[0]);
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

        } else if (use_bdc)
        {
            svd::log->trace("Running BDC SVD with threshold {:.4e} | switchsize {} | size {}", lapacke_svd_threshold, lapacke_svd_switchsize, sizeS);
            int mx = std::max(rowsA,colsA);
            int mn = std::min(rowsA,colsA);
            int lrwork = std::max(1, mn * std::max(5*mn + 7, 2*mx + 2*mn + 1 ));
            int liwork = std::max(1, 8 * std::min(rowsA, colsA));
            std::vector<int>    iwork(static_cast<size_t>(liwork));
            std::vector<double> rwork(static_cast<size_t>(lrwork));
            std::vector<Scalar> work(1);

            U.resize(rowsU, colsU);
            S.resize(sizeS);
            VT.resize(rowsVT, colsVT);

            auto               Ap  = reinterpret_cast<lapack_complex_double *>(A.data());
            auto               Up  = reinterpret_cast<lapack_complex_double *>(U.data());
            auto               VTp = reinterpret_cast<lapack_complex_double *>(VT.data());
            auto               Wp  = reinterpret_cast<lapack_complex_double *>(work.data());

            svd::log->trace("Querying zgesdd");
            info = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, Ap, lda, S.data(), Up, ldu, VTp, ldvt, Wp, -1, rwork.data(), iwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));

            int lwork = static_cast<int>(std::real(work[0]));
            work.resize(static_cast<size_t>(lwork));
            Wp = reinterpret_cast<lapack_complex_double *>(work.data()); // Update the pointer if reallocated

            svd::log->trace("Running zgesdd");
            info = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, Ap, lda, S.data(), Up, ldu, VTp, ldvt, Wp, lwork, rwork.data(),iwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));
        }
        else
        {
            svd::log->trace("Running SVD with threshold {:.4e} | switchsize {} | size {}", lapacke_svd_threshold, lapacke_svd_switchsize, sizeS);
            int lrwork = 5 * std::min(rowsA, colsA);
            std::vector<Scalar> work(1);
            std::vector<double> rwork(lrwork);

            U.resize(rowsU, colsU);
            S.resize(sizeS);
            VT.resize(rowsVT, colsVT);

            auto               Ap  = reinterpret_cast<lapack_complex_double *>(A.data());
            auto               Up  = reinterpret_cast<lapack_complex_double *>(U.data());
            auto               VTp = reinterpret_cast<lapack_complex_double *>(VT.data());
            auto               Wp  = reinterpret_cast<lapack_complex_double *>(work.data());

            svd::log->trace("Querying zgesvd");
            info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, Ap, lda, S.data(), Up, ldu, VTp, ldvt, Wp, -1, rwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));

            int lwork = static_cast<int>(std::real(work[0]));
            work.resize(lwork);
            Wp = reinterpret_cast<lapack_complex_double *>(work.data()); // Update the pointer if reallocated

            svd::log->trace("Running zgesvd");
            info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, Ap, lda, S.data(), Up, ldu, VTp, ldvt, Wp, lwork, rwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", -info));
        }
    }
    svd::log->trace("Truncating singular values");

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
