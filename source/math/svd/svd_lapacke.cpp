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
#include <tid/tid.h>

//#if !defined(NDEBUG)
#include <h5pp/h5pp.h>
#include <math/linalg/matrix.h>
#include <math/linalg/tensor.h>
//#endif
//

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

    MatrixType<Scalar> A = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows, cols);

#if !defined(NDEBUG)
    {
        auto t_debug = tid::tic_scope("debug");
        if(not A.allFinite()) {
            print_matrix_lapacke(mat_ptr, rows, cols);
            throw std::runtime_error("SVD error: matrix has inf's or nan's");
        }
        if(A.isZero(1e-12)) {
            print_matrix_lapacke(mat_ptr, rows, cols, 16);
            throw std::runtime_error("SVD error: matrix is all zeros");
        }
    }
#endif

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

    MatrixType<Scalar> U;
    VectorType<double> S;
    MatrixType<Scalar> VT;

    if constexpr(std::is_same<Scalar, double>::value) {
        if(use_jacobi) {
            svd::log->debug("Running Lapacke Jacobi SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);
            // http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga8767bfcf983f8dc6ef2842029ab25599.html#ga8767bfcf983f8dc6ef2842029ab25599
            // For this routine we need rows > cols
            int                 lwork = std::max(6, rowsA + colsA);
            std::vector<Scalar> work(static_cast<size_t>(lwork));
            //            work.setConstant(0);
            //            work[0] = 1;
            S.resize(sizeS);
            MatrixType<Scalar> V(rowsV, colsV); // Local matrix gets transposed after computation

            svd::log->trace("Running dgejsv");
            auto t_dgesvj = tid::tic_token("dgesvj");
            info          = LAPACKE_dgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, A.data(), lda, S.data(), ldv, V.data(), ldv, work.data(), lwork);
            t_dgesvj.toc();
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: parameter {} is invalid", info));
            if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: could not converge: info {}", info));
            long max_size = std::min(S.size(), rank_max.value());
            long rank     = (S.head(max_size).array() >= threshold).count();
            U             = A.leftCols(rank);
            VT            = V.adjoint().topRows(rank);
        } else if(use_bdc) {
            svd::log->debug("Running Lapacke BDC SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);
            if(threshold < 1e-10) throw std::runtime_error(fmt::format("threshold: {:.4e} < 1e-10", threshold));
            int                 liwork = std::max(1, 8 * std::min(rowsA, colsA));
            std::vector<Scalar> work(1);
            std::vector<int>    iwork(static_cast<size_t>(liwork));

            U.resize(rowsU, colsU);
            S.resize(sizeS);
            VT.resize(rowsVT, colsVT);

            svd::log->trace("Querying dgesdd");
            info = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), -1,
                                       iwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: parameter {} is invalid", info));
            if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: could not converge: info {}", info));

            int lwork = static_cast<int>(work[0]);
            work.resize(static_cast<size_t>(lwork));

            svd::log->trace("Running dgesdd");
            auto t_sdd = tid::tic_token("dgesdd");
            info       = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), lwork,
                                             iwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: parameter {} is invalid", info));
            if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: could not converge: info {}", info));

        } else {
            svd::log->debug("Running Lapacke SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);
            std::vector<Scalar> work(1);

            U.resize(rowsU, colsU);
            S.resize(sizeS);
            VT.resize(rowsVT, colsVT);

            svd::log->trace("Querying dgesvd");
            info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), -1);
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: parameter {} is invalid", info));
            if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: could not converge: info {}", info));

            int lwork = static_cast<int>(work[0]);
            work.resize(static_cast<size_t>(lwork));

            svd::log->trace("Running dgesvd");
            auto t_dgesvd = tid::tic_token("dgesvd");
            info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, work.data(), lwork);
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: parameter {} is invalid", info));
            if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD dgesvd error: could not converge: info {}", info));
        }
    }
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
        if(use_jacobi) {
            svd::log->debug("Running Lapacke Jacobi SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);
            std::vector<Scalar> cwork(1);
            std::vector<double> rwork(1);
            std::vector<int>    iwork(1);

            S.resize(sizeS);
            U.resize(rowsU, colsU);             // Local matrix gets transposed after computation
            MatrixType<Scalar> V(rowsV, colsV); // Local matrix gets transposed after computation

            auto Ap     = reinterpret_cast<lapack_complex_double *>(A.data());
            auto Vp     = reinterpret_cast<lapack_complex_double *>(V.data());
            auto pcwork = reinterpret_cast<lapack_complex_double *>(cwork.data());
            rwork       = {1.0};
            svd::log->trace("Querying zgesvj");
            info = LAPACKE_zgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, Ap, lda, S.data(), ldv, Vp, ldv, pcwork, -1, rwork.data(), -1);

            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD error: parameter {} is invalid", info));
            if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD error: could not converge: info {}", info));

            int lcwork = static_cast<int>(std::real(cwork[0]));
            int lrwork = static_cast<int>(rwork[0]);
            cwork.resize(static_cast<size_t>(lcwork));
            rwork.resize(static_cast<size_t>(lrwork));
            pcwork = reinterpret_cast<lapack_complex_double *>(cwork.data());

            svd::log->trace("Running zgesvj | cwork {} | rwork {}", lcwork, lrwork);
            auto t_zgesvj = tid::tic_token("zgesvj");
            info = LAPACKE_zgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, Ap, lda, S.data(), ldv, Vp, ldv, pcwork, lcwork, rwork.data(), lrwork);
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesvj error: parameter {} is invalid", info));
            if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesvj error: could not converge: info {}", info));

            long max_size = std::min(S.size(), rank_max.value());
            long rank     = (S.head(max_size).array() >= threshold).count();
            U             = A.leftCols(rank);
            VT            = V.adjoint().topRows(rank);
            svd::log->trace("info {} | rwork {} {} {} {} | rank {}", info, rwork[0], rwork[1], rwork[2], rwork[3], rank);

        } else if(use_bdc) {
            svd::log->debug("Running Lapacke BDC SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);
            int                 mx     = std::max(rowsA, colsA);
            int                 mn     = std::min(rowsA, colsA);
            int                 lrwork = std::max(1, mn * std::max(5 * mn + 7, 2 * mx + 2 * mn + 1));
            int                 liwork = std::max(1, 8 * std::min(rowsA, colsA));
            std::vector<int>    iwork(static_cast<size_t>(liwork));
            std::vector<double> rwork(static_cast<size_t>(lrwork));
            std::vector<Scalar> work(1);

            U.resize(rowsU, colsU);
            S.resize(sizeS);
            VT.resize(rowsVT, colsVT);

            auto Ap  = reinterpret_cast<lapack_complex_double *>(A.data());
            auto Up  = reinterpret_cast<lapack_complex_double *>(U.data());
            auto VTp = reinterpret_cast<lapack_complex_double *>(VT.data());
            auto Wp  = reinterpret_cast<lapack_complex_double *>(work.data());

            svd::log->trace("Querying zgesdd");
            info = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, Ap, lda, S.data(), Up, ldu, VTp, ldvt, Wp, -1, rwork.data(), iwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesdd error: parameter {} is invalid", info));
            if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesdd error: could not converge: info {}", info));
            int lwork = static_cast<int>(std::real(work[0]));
            work.resize(static_cast<size_t>(lwork));
            Wp = reinterpret_cast<lapack_complex_double *>(work.data()); // Update the pointer if reallocated
            svd::log->trace("Running zgesdd");
            auto t_zgesdd = tid::tic_token("zgesdd");
            info = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, Ap, lda, S.data(), Up, ldu, VTp, ldvt, Wp, lwork, rwork.data(), iwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesdd error: parameter {} is invalid", info));
            if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesdd error: could not converge: info {}", info));
        } else {
            svd::log->debug("Running Lapacke SVD | threshold {:.4e} | switchsize {} | size {}", threshold, switchsize, sizeS);
            int                 lrwork = 5 * std::min(rowsA, colsA);
            std::vector<Scalar> work(1);
            std::vector<double> rwork(static_cast<size_t>(lrwork));

            U.resize(rowsU, colsU);
            S.resize(sizeS);
            VT.resize(rowsVT, colsVT);

            auto Ap  = reinterpret_cast<lapack_complex_double *>(A.data());
            auto Up  = reinterpret_cast<lapack_complex_double *>(U.data());
            auto VTp = reinterpret_cast<lapack_complex_double *>(VT.data());
            auto Wp  = reinterpret_cast<lapack_complex_double *>(work.data());

            svd::log->trace("Querying zgesvd");
            info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, Ap, lda, S.data(), Up, ldu, VTp, ldvt, Wp, -1, rwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesvd error: parameter {} is invalid", info));
            if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesvd error: could not converge: info {}", info));
            int lwork = static_cast<int>(std::real(work[0]));
            work.resize(static_cast<size_t>(lwork));
            Wp = reinterpret_cast<lapack_complex_double *>(work.data()); // Update the pointer if reallocated

            svd::log->trace("Running zgesvd");
            auto t_zgesvd = tid::tic_token("zgesvd");
            info          = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, Ap, lda, S.data(), Up, ldu, VTp, ldvt, Wp, lwork, rwork.data());
            if(info < 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesvd error: parameter {} is invalid", info));
            if(info > 0) throw std::runtime_error(fmt::format("Lapacke SVD zgesvd error: could not converge: info {}", info));
        }
    }
    svd::log->trace("Truncating singular values");
    if(count) count.value()++;
    long max_size = std::min(S.size(), rank_max.value());
    long rank     = (S.head(max_size).array() >= threshold).count();
    if(rank == S.size()) {
        truncation_error = 0;
    } else {
        truncation_error = S.tail(S.size() - rank).norm();
    }

    bool U_finite   = U.leftCols(rank).allFinite();
    bool S_finite   = S.head(rank).allFinite();
    bool V_finite   = VT.topRows(rank).allFinite();
    bool S_positive = (S.head(rank).array() >= 0).all();

    if(rank <= 0 or rank > rows or info != 0 or not U_finite or not S_finite or not S_positive or not V_finite) {
        if(not A.allFinite()) {
            print_matrix(A.data(), A.rows(), A.cols());
            svd::log->critical("Lapacke SVD error: matrix has inf's or nan's");
        }
        if(A.isZero(1e-12)) {
            print_matrix(A.data(), A.rows(), A.cols(), 16);
            svd::log->critical("Lapacke SVD error: matrix is all zeros");
        }
        if(not S_positive) {
            print_vector(S.data(), rank, 16);
            svd::log->critical("Lapacke SVD error: S is not positive");
        }
        //#if !defined(NDEBUG)
        if(save_fail) {
            auto file       = h5pp::File("svd-failed.h5", h5pp::FilePermission::READWRITE);
            auto group_num  = 0;
            auto group_name = fmt::format("svd_lapacke_{}", group_num);
            while(file.linkExists(group_name)) group_name = fmt::format("svd_lapacke_{}", ++group_num);
            file.writeDataset(A, fmt::format("{}/A", group_name));
            file.writeDataset(U.leftCols(rank), fmt::format("{}/U", group_name));
            file.writeDataset(S.head(rank), fmt::format("{}/S", group_name));
            file.writeDataset(VT.topRows(rank), fmt::format("{}/V", group_name));
            file.writeAttribute(rows, "rows", group_name);
            file.writeAttribute(cols, "cols", group_name);
            file.writeAttribute(rank, "rank", group_name);
            file.writeAttribute(rank_max.value(), "rank_max", group_name);
            file.writeAttribute(info, "info", group_name);
            file.writeAttribute(use_bdc, "use_bdc", group_name);
            file.writeAttribute(threshold, "threshold", group_name);
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

        //#endif
        throw std::runtime_error(fmt::format(FMT_STRING("Lapacke SVD error \n"
                                                        "  svd_threshold    = {:.4e}\n"
                                                        "  Truncation Error = {:.4e}\n"
                                                        "  Rank             = {}\n"
                                                        "  Dims             = ({}, {})\n"
                                                        "  A all finite     : {}\n"
                                                        "  U all finite     : {}\n"
                                                        "  S all finite     : {}\n"
                                                        "  S all positive   : {}\n"
                                                        "  V all finite     : {}\n"
                                                        "  Lapacke info     : {}\n"),
                                             threshold, truncation_error, rank, rows, cols, A.allFinite(), U_finite, S_finite, S_positive, V_finite, info));
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
