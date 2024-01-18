#include "../svd.h"
#include "debug/exceptions.h"
#include "tid/tid.h"
#include <complex>
#include <csignal>
#include <fmt/ranges.h>

#ifndef lapack_complex_float
    #define lapack_complex_float std::complex<float>
#endif
#ifndef lapack_complex_double
    #define lapack_complex_double std::complex<double>
#endif

// complex must be included before lapacke!
#if defined(MKL_AVAILABLE)
    #include <mkl_lapacke.h>
#elif defined(OPENBLAS_AVAILABLE)
    #include <openblas/lapacke.h>
#elif defined(FLEXIBLAS_AVAILABLE)
    #include <flexiblas/lapacke.h>
#else
    #include <lapacke.h>
#endif

namespace svd {
#if defined(NDEBUG)
    static constexpr bool ndebug = true;
#else
    static constexpr bool ndebug = false;
#endif
}

namespace svd {

    namespace internal {
        // These are workspace arrays used by LAPACK which can be reused for the duration of the program.
        // Call clear() to recover the memory space
        std::vector<int>                  iwork;
        std::vector<std::complex<double>> cwork;
        std::vector<double>               rwork;
        void                              clear_lapack() {
            iwork = std::vector<int>();
            cwork = std::vector<std::complex<double>>();
            rwork = std::vector<double>();
            iwork.shrink_to_fit();
            cwork.shrink_to_fit();
            rwork.shrink_to_fit();
        }
    }
}

template<typename Scalar>
std::tuple<svd::MatrixType<Scalar>, svd::VectorType<Scalar>, svd::MatrixType<Scalar>> svd::solver::do_svd_lapacke(const Scalar *mat_ptr, long rows,
                                                                                                                  long cols) const {
    // Setup useful sizes
    int rowsA = static_cast<int>(rows);
    int colsA = static_cast<int>(cols);
    int sizeS = std::min(rowsA, colsA);

    if(rows < cols and (svd_rtn == rtn::gejsv or svd_rtn == rtn::gesvj)) {
        // The jacobi routines needs a tall matrix
        auto t_adj = tid::tic_token("adjoint", tid::highest);
        svd::log->trace("Transposing {}x{} into a tall matrix {}x{}", rows, cols, cols, rows);
        MatrixType<Scalar> A = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows, cols);
        A.adjointInPlace(); // Adjoint directly on a map seems to give a bug?
        // Sanity checks
        if(A.rows() <= 0) throw std::logic_error("SVD error: rows() == 0");
        if(A.cols() <= 0) throw std::logic_error("SVD error: cols() == 0");

        t_adj.toc();
        auto [U, S, VT] = do_svd_lapacke(A.data(), A.rows(), A.cols());
        if(U.rows() != A.rows()) throw except::logic_error("U.rows():{} != A.rows():{}", U.rows(), A.rows());
        if(VT.cols() != A.cols()) throw except::logic_error("VT.cols():{} != A.cols():{}", VT.cols(), A.cols());
        return std::make_tuple(VT.adjoint(), S, U.adjoint());
    }
    auto t_lpk = tid::tic_scope("lapacke", tid::highest);

    // Sanity checks
    if(rows <= 0) throw std::logic_error("SVD error: rows() <= 0");
    if(cols <= 0) throw std::logic_error("SVD error: cols() <= 0");

    MatrixType<Scalar> A = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows, cols); // gets destroyed in some routines
    if(svd_save != svd::save::NONE) save_svd(A);
    if(svd_save == svd::save::FAIL) saveMetaData.A = A;
    //    saveMetaData.svd_is_running = true; // TODO: REMOVE THIS LINE! We don't really want to save it every time!!
    //    saveMetaData.svd_save = save::ALL;  // TODO: REMOVE THIS LINE! We don't really want to save it every time!!
    //    saveMetaData.A = A; // TODO: REMOVE THIS LINE! We don't really want to save it every time!!
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

    int         mx = std::max(rowsA, colsA);
    int         mn = std::min(rowsA, colsA);
    std::string errmsg;
    try {
        // Sanity checks
        if constexpr(!ndebug)
            if(A.isZero(1e-16)) svd::log->warn("Lapacke SVD: A is a zero matrix");
        if(not A.allFinite()) {
            print_matrix(A.data(), A.rows(), A.cols());
            throw std::runtime_error("A has inf's or nan's");
        }
        using namespace svd::internal;
        if constexpr(std::is_same<Scalar, double>::value) {
            auto t_svd = tid::tic_token(fmt::format("d{}{}", enum2sv(svd_rtn), t_suffix), tid::highest);
            if constexpr(!ndebug)
                svd::log->debug("Running Lapacke d{} | truncation limit {:.4e} | switchsize bdc {} | size {}", enum2sv(svd_rtn), truncation_lim,
                                switchsize_gesdd, sizeS);
            switch(svd_rtn) {
                case rtn::gesvd: {
                    U.resize(rowsU, colsU);
                    S.resize(sizeS);
                    VT.resize(rowsVT, colsVT);
                    rwork.resize(1ul);
                    info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, rwork.data(),
                                               -1);
                    if(info != 0) break;
                    int lrwork = static_cast<int>(rwork[0]);
                    rwork.resize(static_cast<size_t>(std::max(1, lrwork)));

                    info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, rwork.data(),
                                               lrwork);
                    break;
                }
                case rtn::gesvj: {
                    // For this routine we need rows > cols
                    S.resize(sizeS);
                    V.resize(rowsV, colsV); // Local matrix gets transposed after computation

                    int lrwork = std::max(6, rowsA + colsA);
                    rwork.resize(static_cast<size_t>(lrwork));
                    info =
                        LAPACKE_dgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, A.data(), lda, S.data(), ldv, V.data(), ldv, rwork.data(), lrwork);
                    if(info < 0) throw except::runtime_error("Lapacke SVD dgesvj error: parameter {} is invalid", info);
                    if(info > 0) throw except::runtime_error("Lapacke SVD dgesvj error: could not converge: info {}", info);
                    U  = A;
                    VT = V.adjoint();
                    break;
                }
                case rtn::gejsv: {
                    // http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga8767bfcf983f8dc6ef2842029ab25599.html#ga8767bfcf983f8dc6ef2842029ab25599
                    // For this routine we need rows > cols
                    S.resize(sizeS);
                    U.resize(rowsU, colsU);
                    V.resize(rowsV, colsV); // Local matrix gets transposed after computation

                    int lrwork = std::max(2 * rowsA + colsA, 6 * colsA + 2 * colsA * colsA);
                    int liwork = std::max(3, rowsA + 3 * colsA);

                    rwork.resize(static_cast<size_t>(lrwork));
                    iwork.resize(static_cast<size_t>(liwork));

                    info = LAPACKE_dgejsv_work(LAPACK_COL_MAJOR, 'F' /* 'R' may also work well */, 'U', 'V', 'N' /* 'R' kills small columns of A */,
                                               'T' /* T/N:  T will transpose if faster. Ignored if A is rectangular */,
                                               'N' /* P/N: P will use perturbation to drown denormalized numbers */, rowsA, colsA, A.data(), lda, S.data(),
                                               U.data(), ldu, V.data(), ldv, rwork.data(), lrwork, iwork.data());
                    if(info < 0) throw except::runtime_error("Lapacke SVD dgejsv error: parameter {} is invalid", info);
                    if(info > 0) throw except::runtime_error("Lapacke SVD dgejsv error: could not converge: info {}", info);
                    VT = V.adjoint();
                    V.resize(0, 0);
                    break;
                }
                case rtn::gesdd: {
                    U.resize(rowsU, colsU);
                    S.resize(sizeS);
                    VT.resize(rowsVT, colsVT);

                    int lrwork = std::max(1, mn * (6 + 4 * mn) + mx);
                    int liwork = std::max(1, 8 * mn);
                    rwork.resize(static_cast<size_t>(lrwork));
                    iwork.resize(static_cast<size_t>(liwork));

                    info = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, rwork.data(), -1,
                                               iwork.data());
                    if(info < 0) throw except::runtime_error("Lapacke SVD dgesdd error: parameter {} is invalid", info);
                    if(info > 0) throw except::runtime_error("Lapacke SVD dgesdd error: could not converge: info {}", info);

                    lrwork = static_cast<int>(rwork[0]);
                    rwork.resize(static_cast<size_t>(std::max(1, lrwork)));

                    auto t_sdd = tid::tic_token(fmt::format("dgesdd{}", t_suffix), tid::highest);
                    info       = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, rwork.data(),
                                                     lrwork, iwork.data());
                    if(info < 0) throw except::runtime_error("Lapacke SVD dgesdd error: parameter {} is invalid", info);
                    if(info > 0) throw except::runtime_error("Lapacke SVD dgesdd error: could not converge: info {}", info);
                    break;
                }
                case rtn::gesvdx: {
                    U.resize(rowsU, colsU);
                    S.resize(sizeS);
                    VT.resize(rowsVT, colsVT);

                    double vl    = std::min(1e-10, truncation_lim / 5.0);
                    double vu    = std::max(1e+10, truncation_lim / 5.0);
                    int    il    = 1;
                    int    iu    = std::min<int>(sizeS, static_cast<int>(rank_max));
                    char   range = rank_max < sizeS ? 'I' : 'V';

                    if(svdx_select.has_value()) {
                        if(std::holds_alternative<svdx_indices_t>(svdx_select.value())) {
                            auto sel = std::get<svdx_indices_t>(svdx_select.value());
                            iu       = std::min<int>(static_cast<int>(sel.iu), static_cast<int>(rank_max));
                            il       = std::min<int>(static_cast<int>(sel.il), static_cast<int>(rank_max));
                            range    = 'I';
                        } else if(std::holds_alternative<svdx_values_t>(svdx_select.value())) {
                            auto sel = std::get<svdx_values_t>(svdx_select.value());
                            if(std::isfinite(sel.vl)) vl = sel.vl;
                            if(std::isfinite(sel.vu)) vu = sel.vu;
                            range = 'V';
                        }
                    }

                    int ns     = 0;
                    int lrwork = std::max(1, mn * 2 + mx);
                    int liwork = std::max(1, 12 * mn);
                    rwork.resize(static_cast<size_t>(lrwork));
                    iwork.resize(static_cast<size_t>(liwork));

                    info = LAPACKE_dgesvdx_work(LAPACK_COL_MAJOR, 'V', 'V', range, rowsA, colsA, A.data(), lda, vl, vu, il, iu, &ns, S.data(), U.data(), ldu,
                                                VT.data(), ldvt, rwork.data(), -1, iwork.data());
                    if(info != 0) break;

                    lrwork = static_cast<int>(std::real(rwork[0]));
                    rwork.resize(static_cast<size_t>(std::max(1, lrwork)));

                    info = LAPACKE_dgesvdx_work(LAPACK_COL_MAJOR, 'V', 'V', 'V', rowsA, colsA, A.data(), lda, vl, vu, il, iu, &ns, S.data(), U.data(), ldu,
                                                VT.data(), ldvt, rwork.data(), lrwork, iwork.data());

                    // Select the computed region
                    U  = U.leftCols(ns).eval();
                    S  = S.head(ns).eval(); // Not all calls to do_svd need normalized S, so we do not normalize here!
                    VT = VT.topRows(ns).eval();
                    break;
                }
                default: throw std::logic_error("invalid case for enum svd::rtn");
            }
            if(info < 0) throw except::runtime_error("Lapacke SVD d{} error: parameter {} is invalid", enum2sv(svd_rtn), info);
            if(info > 0) throw except::runtime_error("Lapacke SVD d{} error: could not converge: info {}", enum2sv(svd_rtn), info);
        } else if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
            auto t_svd = tid::tic_token(fmt::format("z{}{}", enum2sv(svd_rtn), t_suffix), tid::highest);
            if constexpr(!ndebug)
                svd::log->debug("Running Lapacke z{} | truncation limit {:.4e} | switchsize bdc {} | size {}", enum2sv(svd_rtn), truncation_lim,
                                switchsize_gesdd, sizeS);
            switch(svd_rtn) {
                case rtn::gesvd: {
                    U.resize(rowsU, colsU);
                    S.resize(sizeS);
                    VT.resize(rowsVT, colsVT);

                    int lcwork = 1;
                    int lrwork = std::max(1, 5 * mn);
                    cwork.resize(static_cast<size_t>(lcwork));
                    rwork.resize(static_cast<size_t>(lrwork));

                    info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, cwork.data(),
                                               -1, rwork.data());
                    if(info != 0) break;
                    lcwork = static_cast<int>(std::real(cwork[0]));
                    cwork.resize(static_cast<size_t>(std::max(1, lcwork)));
                    info = LAPACKE_zgesvd_work(LAPACK_COL_MAJOR, 'S', 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, cwork.data(),
                                               lcwork, rwork.data());
                    break;
                }
                case rtn::gesvj: {
                    S.resize(sizeS);
                    V.resize(rowsV, colsV); // Local matrix gets transposed after computation

                    cwork.resize(1);
                    rwork.resize(6ul);

                    info = LAPACKE_zgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, A.data(), lda, S.data(), ldv, V.data(), ldv, cwork.data(), -1,
                                               rwork.data(), -1);
                    if(info != 0) break;
                    int lcwork = static_cast<int>(std::real(cwork[0]));
                    int lrwork = static_cast<int>(rwork[0]);
                    cwork.resize(static_cast<size_t>(std::max(1, lcwork)));
                    rwork.resize(static_cast<size_t>(std::max(6, lrwork)));

                    info = LAPACKE_zgesvj_work(LAPACK_COL_MAJOR, 'G', 'U', 'V', rowsA, colsA, A.data(), lda, S.data(), ldv, V.data(), ldv, cwork.data(), lcwork,
                                               rwork.data(), lrwork);
                    if(info != 0) break;
                    U  = A;
                    VT = V.adjoint();
                    V.resize(0, 0);
                    break;
                }
                case rtn::gejsv: {
                    // http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga8767bfcf983f8dc6ef2842029ab25599.html#ga8767bfcf983f8dc6ef2842029ab25599
                    // For this routine we need rows > cols
                    S.resize(sizeS);
                    U.resize(rowsU, colsU);
                    V.resize(rowsV, colsV); // Local matrix gets transposed after computation
                    cwork.resize(std::max(2ul, static_cast<size_t>(5 * rowsA + 2 * rowsA * rowsA)));
                    rwork.resize(std::max(7ul, static_cast<size_t>(2 * colsA)));
                    iwork.resize(std::max(4ul, static_cast<size_t>(2 * rowsA + colsA)));

                    info = LAPACKE_zgejsv_work(LAPACK_COL_MAJOR, 'F' /* 'R' may also work well */, 'U', 'V', 'N' /* 'R' kills small columns of A */,
                                               'T' /* T/N:  T will transpose if faster. Ignored if A is rectangular */,
                                               'N' /* P/N: P will use perturbation to drown denormalized numbers */, rowsA, colsA, A.data(), lda, S.data(),
                                               U.data(), ldu, V.data(), ldv, cwork.data(), -1, rwork.data(), -1, iwork.data());
                    if(info != 0) break;

                    int lcwork = static_cast<int>(std::real(cwork[0]));
                    int lrwork = static_cast<int>(rwork[0]);
                    int liwork = static_cast<int>(iwork[0]);
                    cwork.resize(static_cast<size_t>(std::max(2, lcwork)));
                    rwork.resize(static_cast<size_t>(std::max(7, lrwork)));
                    iwork.resize(static_cast<size_t>(std::max({4, 2 * rowsA + colsA, liwork})));

                    info = LAPACKE_zgejsv_work(LAPACK_COL_MAJOR, 'F' /* 'R' may also work well */, 'U', 'V', 'N' /* 'R' kills small columns of A */,
                                               'T' /* T/N:  T will transpose if faster. Ignored if A is rectangular */,
                                               'N' /* P/N: P will use perturbation to drown denormalized numbers */, rowsA, colsA, A.data(), lda, S.data(),
                                               U.data(), ldu, V.data(), ldv, cwork.data(), lcwork, rwork.data(), lrwork, iwork.data());
                    if(info != 0) break;
                    VT = V.adjoint();
                    V.resize(0, 0);
                    break;
                }
                case rtn::gesdd: {
                    U.resize(rowsU, colsU);
                    S.resize(sizeS);
                    VT.resize(rowsVT, colsVT);
                    int lcwork = 1;
                    int lrwork = std::max(1, mn * std::max(5 * mn + 7, 2 * mx + 2 * mn + 1));
                    int liwork = std::max(1, 8 * mn);

                    cwork.resize(static_cast<size_t>(lcwork));
                    rwork.resize(static_cast<size_t>(lrwork));
                    iwork.resize(static_cast<size_t>(liwork));
                    info = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, cwork.data(), -1,
                                               rwork.data(), iwork.data());
                    if(info != 0) break;

                    lcwork = static_cast<int>(std::real(cwork[0]));
                    cwork.resize(static_cast<size_t>(std::max(1, lcwork)));

                    info = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR, 'S', rowsA, colsA, A.data(), lda, S.data(), U.data(), ldu, VT.data(), ldvt, cwork.data(),
                                               lcwork, rwork.data(), iwork.data());
                    break;
                }
                case rtn::gesvdx: {
                    U.resize(rowsU, colsU);
                    S.resize(sizeS);
                    VT.resize(rowsVT, colsVT);

                    double vl    = std::min(1e10, truncation_lim / 5.0);
                    double vu    = std::max(1e10, truncation_lim / 5.0);
                    int    il    = 1;
                    int    iu    = std::min<int>(sizeS, static_cast<int>(rank_max));
                    char   range = rank_max < sizeS ? 'I' : 'V';

                    if(svdx_select.has_value()) {
                        if(std::holds_alternative<svdx_indices_t>(svdx_select.value())) {
                            auto sel = std::get<svdx_indices_t>(svdx_select.value());
                            iu       = std::min<int>(static_cast<int>(sel.iu), static_cast<int>(rank_max));
                            il       = std::min<int>(static_cast<int>(sel.il), static_cast<int>(rank_max));
                            range    = 'I';
                        } else if(std::holds_alternative<svdx_values_t>(svdx_select.value())) {
                            auto sel = std::get<svdx_values_t>(svdx_select.value());
                            if(std::isfinite(sel.vl)) vl = sel.vl;
                            if(std::isfinite(sel.vu)) vu = sel.vu;
                            range = 'V';
                        }
                    }
                    int ns     = 0;
                    int lcwork = std::max(1, mn * 2 + mx);
                    int lrwork = mn * (mn * 2 + 15 * mn);
                    int liwork = 12 * mn;

                    cwork.resize(static_cast<size_t>(lcwork));
                    rwork.resize(static_cast<size_t>(lrwork));
                    iwork.resize(static_cast<size_t>(liwork));

                    info = LAPACKE_zgesvdx_work(LAPACK_COL_MAJOR, 'V', 'V', range, rowsA, colsA, A.data(), lda, vl, vu, il, iu, &ns, S.data(), U.data(), ldu,
                                                VT.data(), ldvt, cwork.data(), -1, rwork.data(), iwork.data());
                    if(info != 0) break;

                    lcwork = static_cast<int>(std::real(cwork[0]));
                    cwork.resize(static_cast<size_t>(std::max(1, lcwork)));

                    info = LAPACKE_zgesvdx_work(LAPACK_COL_MAJOR, 'V', 'V', 'V', rowsA, colsA, A.data(), lda, vl, vu, il, iu, &ns, S.data(), U.data(), ldu,
                                                VT.data(), ldvt, cwork.data(), lcwork, rwork.data(), iwork.data());

                    // Select the computed region
                    U  = U.leftCols(ns).eval();
                    S  = S.head(ns).eval(); // Not all calls to do_svd need normalized S, so we do not normalize here!
                    VT = VT.topRows(ns).eval();
                    break;
                }
                default: throw std::logic_error("invalid case for enum svd::rtn");
            }
            if(info < 0) throw except::runtime_error("z{} error: parameter {} is invalid", enum2sv(svd_rtn), info);
            if(info > 0) throw except::runtime_error("z{} error: could not converge: info {}", enum2sv(svd_rtn), info);
        }

        svd::log->trace("Truncating singular values");
        auto max_size                    = S.nonZeros();
        std::tie(rank, truncation_error) = get_rank_from_truncation_error(S.head(max_size));
        // Do the truncation
        U  = U.leftCols(rank).eval();
        S  = S.head(rank).eval(); // Not all calls to do_svd need normalized S, so we do not normalize here!
        VT = VT.topRows(rank).eval();

        // Sanity checks
        if(not U.allFinite()) {
            print_matrix(U.data(), U.rows(), U.cols());
            throw except::runtime_error("U has inf's or nan's");
        }

        if(not VT.allFinite()) {
            print_matrix(VT.data(), VT.rows(), VT.cols());
            throw except::runtime_error("VT has inf's or nan's");
        }
        if(not S.allFinite()) {
            print_vector(S.data(), rank, 16);
            throw except::runtime_error("S has inf's or nan's");
        }
        if(not(S.array() >= 0).all()) {
            print_vector(S.data(), rank, 16);
            throw except::runtime_error("S is not positive");
        }
    } catch(const except::runtime_error &ex) {
        // #if !defined(NDEBUG)
        if(svd_save == svd::save::FAIL) {
            saveMetaData.U                = U;
            saveMetaData.S                = S;
            saveMetaData.VT               = VT;
            saveMetaData.rank             = rank;
            saveMetaData.truncation_error = truncation_error;
            saveMetaData.info             = info;
        }

        save_svd(); // Used on failure only if svd_save == svd::save::FAIL
        throw except::runtime_error("Lapacke SVD error \n"
                                    "  Singular values  = {::.5e}\n"
                                    "  Truncation Error = {:.4e}\n"
                                    "  Rank             = {}\n"
                                    "  Dims             = ({}, {})\n"
                                    "  Lapacke info     : {}\n"
                                    "  Error message    : {}\n",
                                    S, truncation_error, rank, rows, cols, info, ex.what());
    }
    // TODO: REMOVE THE SCOPE BELOW!
    //    if(std::min({rows, cols}) >= 2000) {
    //        svd::log->info("Dims ({},{}) | rank {} | discarded {} | norm {}", rows, cols, rank, std::min({rows, cols}) - rank, S.norm());
    ////        saveMetaData.U                = U;
    //        saveMetaData.S                = S;
    ////        saveMetaData.VT               = VT;
    //        saveMetaData.rank             = rank;
    //        saveMetaData.truncation_error = truncation_error;
    //        saveMetaData.info             = info;
    //        save_svd();
    //    }
    save_svd<Scalar>(U, S, VT, info);
    saveMetaData = svd::internal::SaveMetaData{}; // Clear
    svd::log->trace(
        "SVD with Lapacke finished successfully | truncation limit {:<8.2e} | rank {:<4} | rank_max {:<4} | {:>4} x {:<4} | trunc {:8.2e}, time {:8.2e}",
        truncation_lim, rank, rank_max, rows, cols, truncation_error, t_lpk->get_last_interval());
    return std::make_tuple(U, S, VT);
}
using real = double;
using cplx = std::complex<double>;

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd_lapacke for type 'double'
template std::tuple<svd::MatrixType<real>, svd::VectorType<real>, svd::MatrixType<real>> svd::solver::do_svd_lapacke(const real *, long, long) const;

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd_lapacke for type 'std::complex<double>'
template std::tuple<svd::MatrixType<cplx>, svd::VectorType<cplx>, svd::MatrixType<cplx>> svd::solver::do_svd_lapacke(const cplx *, long, long) const;
