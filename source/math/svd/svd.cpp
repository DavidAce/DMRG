#include "math/svd.h"
#include "tid/tid.h"
#include <Eigen/QR>

#if defined(_OPENMP)

    #include <omp.h>

#endif

svd::solver::solver() { setLogLevel(2); }

svd::config::config(long rank_max_) : rank_max(rank_max_) {}

svd::config::config(double truncation_lim_) : truncation_limit(truncation_lim_) {}

svd::config::config(long rank_max_, double truncation_lim_) : rank_max(rank_max_), truncation_limit(truncation_lim_) {}

svd::config::config(std::optional<long> rank_max_) : rank_max(rank_max_) {}

svd::config::config(std::optional<double> truncation_lim_) : truncation_limit(truncation_lim_) {}

svd::config::config(std::optional<long> rank_max_, std::optional<double> truncation_lim_) : rank_max(rank_max_), truncation_limit(truncation_lim_) {}

std::string svd::config::to_string() {
    /* clang-format off */
    std::string msg;
    if (rank_max) msg.append(fmt::format(" | rank_max {}", rank_max.value()));
    if (rank_min) msg.append(fmt::format(" | rank_min {}", rank_min.value()));
    if (truncation_limit) msg.append(fmt::format(" | truncation_lim {:.2e}", truncation_limit.value()));
    if (switchsize_gejsv) msg.append(fmt::format(" | switchsize_gejsv {}", switchsize_gejsv.value()));
    if (switchsize_gesvd) msg.append(fmt::format(" | switchsize_gesvd {}", switchsize_gesvd.value()));
    if (switchsize_gesdd) msg.append(fmt::format(" | switchsize_gesdd {}", switchsize_gesdd.value()));
    if (svdx_select and std::holds_alternative<svdx_indices_t>(svdx_select.value())) {
        auto sel = std::get<svdx_indices_t>(svdx_select.value());
        msg.append(fmt::format(" | svdx_select {}-{}", sel.il, sel.iu));
    }
    if (svdx_select and std::holds_alternative<svdx_values_t>(svdx_select.value())) {
        auto sel = std::get<svdx_values_t>(svdx_select.value());
        msg.append(fmt::format(" | svdx_select {}-{}", sel.vl, sel.vu));
    }
    if (loglevel) msg.append(fmt::format(" | loglevel {}", loglevel.value()));
    if (svd_lib) msg.append(fmt::format(" | svd_lib {}", enum2sv(svd_lib.value())));
    if (svd_rtn) msg.append(fmt::format(" | svd_rtn {}", enum2sv(svd_rtn.value())));
    if (svd_save) msg.append(fmt::format(" | svd_save {}", enum2sv(svd_save.value())));
    if (benchmark) msg.append(fmt::format(" | benchmark {}", benchmark.value()));
    return msg.empty() ? msg : "svd settings" + msg;
    /* clang-format on */
}

void svd::solver::copy_config(const svd::config &svd_cfg) {
    if(svd_cfg.rank_max) rank_max = svd_cfg.rank_max.value();
    if(svd_cfg.rank_min) rank_min = svd_cfg.rank_min.value();
    if(svd_cfg.truncation_limit) truncation_lim = svd_cfg.truncation_limit.value();
    if(svd_cfg.switchsize_gejsv) switchsize_gejsv = svd_cfg.switchsize_gejsv.value();
    if(svd_cfg.switchsize_gesvd) switchsize_gesvd = svd_cfg.switchsize_gesvd.value();
    if(svd_cfg.switchsize_gesdd) switchsize_gesdd = svd_cfg.switchsize_gesdd.value();
    if(svd_cfg.svdx_select) svdx_select = svd_cfg.svdx_select.value();
    if(svd_cfg.loglevel) setLogLevel(svd_cfg.loglevel.value());
    if(svd_cfg.svd_lib) svd_lib = svd_cfg.svd_lib.value();
    if(svd_cfg.svd_rtn) svd_rtn = svd_cfg.svd_rtn.value();
    if(svd_cfg.svd_save) svd_save = svd_cfg.svd_save.value();
    if(svd_cfg.benchmark) benchmark = svd_cfg.benchmark.value();
}

svd::solver::solver(const svd::config &svd_cfg) : solver() { copy_config(svd_cfg); }

svd::solver::solver(std::optional<svd::config> svd_cfg) : solver() {
    if(svd_cfg) copy_config(svd_cfg.value());
}

void svd::solver::set_config(const svd::config &svd_cfg) { copy_config(svd_cfg); }

void svd::solver::set_config(std::optional<svd::config> svd_cfg) {
    if(svd_cfg) copy_config(svd_cfg.value());
}

void svd::solver::setLogLevel(int logLevel) {
    if(!log) {
        std::string name = "svd";
#if defined(_OPENMP)
        if(omp_in_parallel()) name = fmt::format("svd-{}", omp_get_thread_num());
#endif
        log = spdlog::get(name);
        if(!log) {
            log = spdlog::stdout_color_mt(name, spdlog::color_mode::always);
            log->set_pattern("[%Y-%m-%d %H:%M:%S.%e][%n]%^[%=8l]%$ %v");
        }
    } else {
        if(logLevel != log->level()) { log->set_level(static_cast<spdlog::level::level_enum>(logLevel)); }
    }
}

double svd::solver::get_truncation_error() const { return truncation_error; }

long svd::solver::get_rank() const { return rank; }

long long svd::solver::get_count() { return count; }

/*! \brief Performs SVD on a matrix
 *  This function is defined in cpp to avoid long compilation times when having Eigen::BDCSVD included everywhere in headers.
 *  Performs rigorous checks to ensure stability of DMRG.
 *  In some cases Eigen::BCDSVD/JacobiSVD will fail with segfault. Here we use a patched version of Eigen that throws an error
 *  instead so we get a chance to catch it and use lapack svd instead.
 *   \param mat_ptr Pointer to the matrix. Supported are double * and std::complex<double> *
 *   \param rows Rows of the matrix
 *   \param cols Columns of the matrix
 *   \param svd_cfg Optional overrides to default svd configuration
 *   \return The U, S, and V matrices (with S as a vector) extracted from the Eigen::BCDSVD SVD object.
 */
template<typename Scalar>
std::tuple<svd::MatrixType<Scalar>, svd::VectorType<Scalar>, svd::MatrixType<Scalar>> svd::solver::do_svd_ptr(const Scalar *mat_ptr, long rows, long cols,
                                                                                                              const svd::config &svd_cfg) {
    //    auto t_svd = tid::tic_scope("svd", tid::level::highest);

    copy_config(svd_cfg);
    auto sizeS    = std::min(rows, cols);
    long rank_lim = rank_max > 0 ? std::min(sizeS, rank_max) : sizeS;

    // Resolve geauto
    if(svd_cfg.svd_rtn == svd::rtn::geauto or svd_rtn == svd::rtn::geauto) {
        svd_rtn = svd::rtn::gesvj;
        if(switchsize_gejsv != -1ul and std::cmp_greater_equal(sizeS, switchsize_gejsv)) svd_rtn = svd::rtn::gejsv;
        if(switchsize_gesvd != -1ul and std::cmp_greater_equal(sizeS, switchsize_gesvd)) svd_rtn = svd::rtn::gesvd;
        if(switchsize_gesdd != -1ul and std::cmp_greater_equal(sizeS, switchsize_gesdd)) svd_rtn = svd::rtn::gesdd;

        if(svd_rtn == rtn::gesdd) {
            bool is_rank_low   = std::cmp_greater_equal(sizeS,
                                                        rank_lim * 8); // Will keep at least 25% of the singular values
            bool is_rank_lower = std::cmp_greater_equal(sizeS,
                                                        rank_lim * 16); // Will keep at least 10% of the singular values
            if(is_rank_low) { svd_rtn = svd::rtn::gesvdx; }
            if(is_rank_lower) { svd_rtn = svd::rtn::gersvd; }
        }
        //        log->info("sizeS = {} | lim {} | {} {} {} | {} ", sizeS, rank_lim, switchsize_gejsv, switchsize_gesvd, switchsize_gesdd,
        //        enum2sv(svd_rtn));
    }
    if(svd_save != svd::save::NONE) {
        saveMetaData.rank_max         = rank_max;
        saveMetaData.rank_min         = rank_min;
        saveMetaData.truncation_lim   = truncation_lim;
        saveMetaData.switchsize_gejsv = switchsize_gejsv;
        saveMetaData.switchsize_gesvd = switchsize_gesvd;
        saveMetaData.switchsize_gesdd = switchsize_gesdd;
        saveMetaData.svd_lib          = svd_lib;
        saveMetaData.svd_rtn          = svd_rtn;
        saveMetaData.svd_save         = svd_save;
        if(not saveMetaData.at_quick_exit) {
            std::at_quick_exit(save_svd);
            saveMetaData.at_quick_exit = true;
        }
        saveMetaData.svd_is_running = true;
    }

#pragma omp atomic
    count++;
    switch(svd_lib) {
        case svd::lib::lapacke: {
            try {
                if(svd_rtn == svd::rtn::gersvd)
                    return do_svd_rsvd(mat_ptr, rows, cols);
                else
                    return do_svd_lapacke(mat_ptr, rows, cols);
            } catch(const std::exception &ex) {
                try {
                    log->warn("{} {} failed to perform SVD: {} | Trying Lapacke gejsv", enum2sv(svd_lib), enum2sv(svd_rtn), std::string_view(ex.what()));
                    auto svd_rtn_backup = svd_rtn; // Restore after
                    auto svd_log_level  = log->level();
                    log->set_level(spdlog::level::trace);
                    svd_rtn         = rtn::gejsv;
                    auto [U, S, VT] = do_svd_lapacke(mat_ptr, rows, cols);
                    svd_rtn         = svd_rtn_backup;
                    log->set_level(svd_log_level);
                    return {U, S, VT};
                } catch(const std::exception &ex) {
                    log->warn("{} {} failed to perform SVD: {} | Trying Eigen JacobiSVD", enum2sv(svd_lib), enum2sv(svd_rtn), std::string_view(ex.what()));
                    auto svd_rtn_backup = svd_rtn; // Restore after
                    auto svd_log_level  = log->level();
                    log->set_level(spdlog::level::trace);
                    svd_rtn         = rtn::gejsv;
                    auto [U, S, VT] = do_svd_eigen(mat_ptr, rows, cols);
                    svd_rtn         = svd_rtn_backup;
                    log->set_level(svd_log_level);
                    return {U, S, VT};
                }
            }
            break;
        }
        case svd::lib::eigen: {
            try {
                if(svd_rtn == svd::rtn::gersvd)
                    return do_svd_rsvd(mat_ptr, rows, cols);
                else
                    return do_svd_eigen(mat_ptr, rows, cols);
            } catch(const std::exception &ex) {
                log->warn("{} {} failed to perform SVD: {} | Trying Eigen", enum2sv(svd_lib), enum2sv(svd_rtn), std::string_view(ex.what()));
                return do_svd_lapacke(mat_ptr, rows, cols);
            }
            break;
        }
        default: throw std::logic_error("Unrecognized svd library");
    }
    throw std::logic_error("Unrecognized svd library");
}

using real = double;
using cplx = std::complex<double>;

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'double'
template std::tuple<svd::MatrixType<real>, svd::VectorType<real>, svd::MatrixType<real>> svd::solver::do_svd_ptr(const real *, long, long, const svd::config &);

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'std::complex<double>'
template std::tuple<svd::MatrixType<cplx>, svd::VectorType<cplx>, svd::MatrixType<cplx>> svd::solver::do_svd_ptr(const cplx *, long, long, const svd::config &);

template<typename Scalar>
void svd::solver::print_matrix([[maybe_unused]] const Scalar *mat_ptr, [[maybe_unused]] long rows, [[maybe_unused]] long cols,
                               [[maybe_unused]] std::string_view tag, [[maybe_unused]] long dec) const {
#if !defined(NDEBUG)
    auto A = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(mat_ptr, rows, cols);
    log->warn("Matrix [{}] with dimensions {}x{}\n", tag, rows, cols);
    for(long r = 0; r < A.rows(); r++) {
        if constexpr(std::is_same_v<Scalar, std::complex<double>>)
            for(long c = 0; c < A.cols(); c++) fmt::print("({1:.{0}f},{2:+.{0}f}) ", dec, std::real(A(r, c)), std::imag(A(r, c)));
        else
            for(long c = 0; c < A.cols(); c++) fmt::print("{1:.{0}f} ", dec, A(r, c));
        fmt::print("\n");
    }
#endif
}

template<typename Scalar>
void svd::solver::print_vector([[maybe_unused]] const Scalar *vec_ptr, [[maybe_unused]] long size, [[maybe_unused]] std::string_view tag,
                               [[maybe_unused]] long dec) const {
#if !defined(NDEBUG)
    auto V = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(vec_ptr, size);
    log->warn("Vector [{}] with size {}\n", tag, size);
    if constexpr(std::is_same_v<Scalar, std::complex<double>>)
        for(long i = 0; i < V.size(); i++) fmt::print("({1:.{0}f},{2:+.{0}f})\n", dec, std::real(V[i]), std::imag(V[i]));
    else
        for(long i = 0; i < V.size(); i++) fmt::print("{1:.{0}f}\n", dec, V[i]);
#endif
}

template void svd::solver::print_matrix<real>(const real *vec_ptr, long rows, long cols, std::string_view tag, long dec) const;

template void svd::solver::print_matrix<cplx>(const cplx *vec_ptr, long rows, long cols, std::string_view tag, long dec) const;

template void svd::solver::print_vector<real>(const real *vec_ptr, long size, std::string_view tag, long dec) const;

template void svd::solver::print_vector<cplx>(const cplx *vec_ptr, long size, std::string_view tag, long dec) const;


// template<typename Scalar>
std::pair<long, double> svd::solver::get_rank_from_truncation_error(const VectorType<double> &S) const {
    //        assert(std::abs(S.norm() - 1.0) < 1e-10); // make sure this is normalized
    VectorType<double> truncation_errors(S.size() + 1);
    for(long s = 0; s <= S.size(); s++) { truncation_errors[s] = S.bottomRows(S.size() - s).norm(); } // Last one should be zero, i.e. no truncation
    auto rank_    = (truncation_errors.array() >= truncation_lim).count();
    auto rank_lim = S.size();
    if(rank_max > 0) rank_lim = std::min(S.size(), rank_max);
    rank_ = std::min(rank_, rank_lim);
    if(rank_min > 0) rank_ = std::max(rank_, std::min(S.size(),
                                                      rank_min)); // Make sure we don't overtruncate in some cases (e.g. when stashing)

    //    tools::log->info("Size {} | Rank {} | Rank limit {} | truncation error limit {:8.5e} | error {:8.5e}", S.size(), rank_, rank_lim,
    //                     truncation_lim, truncation_errors[rank_]);
    //    tools::log->info("Size {} | Rank {} | Rank limit {} | truncation error limit {:8.5e} | error {:8.5e} truncation errors: {:8.5e}", S.size(), rank_,
    //    rank_lim,
    //                     truncation_lim, truncation_errors[rank_], fmt::join(truncation_errors, ", "));
    if(rank_ <= 0) {
        if(log)
            log->error("Size {} | Rank {} | Rank limit {} | truncation error limit {:8.2e} | error {:8.2e} truncation errors: {:8.2e}", S.size(), rank_,
                       rank_lim, truncation_lim, truncation_errors[rank_], fmt::join(truncation_errors, ", "));
        throw std::logic_error("rank <= 0");
    }
    return {rank_, truncation_errors[rank_]};
}

template<typename Scalar>
auto dropfilter(svd::MatrixType<Scalar> &U, svd::VectorType<Scalar> &S, svd::MatrixType<Scalar> &V, double threshold, double min_log_drop_size) {
    long   drop_rank = 0;
    double drop_size = 0;
    for(long idx = 0; idx < S.size(); ++idx) {
        if(idx + 1 < S.size()) drop_size = std::abs(std::abs(std::log10(std::real(S[idx]))) - std::abs(std::log10(std::real(S[idx + 1]))));
        if(std::real(S[idx]) >= threshold and drop_size <= min_log_drop_size) {
            drop_rank = idx + 1; // Keep 0,1,2...idx+1 singular values
        } else
            break;
    }
    // fmt::print("drop_rank: {} | drop_size {}\n", drop_rank, drop_size);
    if(drop_rank > 0) {
        U = U.leftCols(drop_rank).eval();
        S = S.head(drop_rank).eval();
        V = V.topRows(drop_rank).eval();
    }
    return S.size();
}

template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 4>, Eigen::Tensor<Scalar, 2>> svd::solver::split_mpo_l2r(const Eigen::Tensor<Scalar, 4> &mpo, const svd::config &svd_cfg) {
    /*
     * Compress an MPO left to right using SVD as described in https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.035129
     *
     *          (2) d
     *             |
     *  (0) m ---[mpo]--- (1) m
     *             |
     *         (3) d
     *
     * is shuffled into
     *
     *          (0) d
     *             |
     *  (2) m ---[mpo]--- (3) m
     *             |
     *         (1) d
     *
     * and reshaped like
     *
     * ddm (012) ---[mpo]--- (3) m
     *
     * and subjected to the typical SVD so mpo = USV. This is then reshaped back into
     *
     *            d
     *            |
     *     m ---[mpo]---  m'   m'---[S]---'m m'---[V]---m
     *            |
     *            d
     *
     * where hopefully m' < m and the transfer matrix T = SV is multiplied onto the mpo on the right later.
     *
     * To stablize the compression, it is useful to insert avgS *  1/avgS, where one factor is put into U and the other into S.
     *
     *
     */

    auto                     dim0      = mpo.dimension(2);
    auto                     dim1      = mpo.dimension(3);
    auto                     dim2      = mpo.dimension(0);
    auto                     dim3      = mpo.dimension(1);
    auto                     dim_ddm   = dim0 * dim1 * dim2;
    Eigen::Tensor<Scalar, 2> mpo_rank2 = mpo.shuffle(tenx::array4{2, 3, 0, 1}).reshape(tenx::array2{dim_ddm, dim3});
    auto [U, S, VT]                     = do_svd_ptr(mpo_rank2.data(), mpo_rank2.dimension(0), mpo_rank2.dimension(1), svd_cfg);
    auto Smin = S.real().minCoeff();
    auto Smax = S.real().maxCoeff();
    // Stabilize by inserting avgS *  1/avgS
    auto avgS = num::next_power_of_two<double>(S.head(S.nonZeros()).real().mean()); // Nearest power of two larger than S.mean()
    if(avgS > 1) {
        S /= avgS;
        std::tie(rank, truncation_error) = get_rank_from_truncation_error(S.real().head(S.real().nonZeros()));
        U  = U.leftCols(rank).eval();
        S  = S.head(rank).eval();
        VT = VT.topRows(rank).eval();
        U *= avgS;
    }
    // rank = dropfilter(U, S, V, svd_cfg.truncation_limit.value_or(1e-16), 8);
    VT    = S.asDiagonal() * VT; // Rescaled singular values
    fmt::print("S l2r min {:.5e} | max {:.5e} -->  min {:.5e} | max {:.5e} | size {}\n", Smin, Smax, S.real().minCoeff(), S.real().maxCoeff(), S.size());
    /* clang-format off */
    return std::make_tuple(
            tenx::TensorMap(U).reshape(tenx::array4{dim0, dim1, dim2, rank}).shuffle(
                    tenx::array4{2, 3, 0, 1}).template cast<Scalar>(),
            tenx::TensorMap(VT).template cast<Scalar>());
    /* clang-format on */
}

template std::tuple<Eigen::Tensor<real, 4>, Eigen::Tensor<real, 2>> svd::solver::split_mpo_l2r(const Eigen::Tensor<real, 4> &mpo, const svd::config &svd_cfg);

template std::tuple<Eigen::Tensor<cplx, 4>, Eigen::Tensor<cplx, 2>> svd::solver::split_mpo_l2r(const Eigen::Tensor<cplx, 4> &mpo, const svd::config &svd_cfg);

template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 2>, Eigen::Tensor<Scalar, 4>> svd::solver::split_mpo_r2l(const Eigen::Tensor<Scalar, 4> &mpo, const svd::config &svd_cfg) {
    /*
     * Splits an MPO right to left using SVD as described in https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.035129
     *
     *          (2) d
     *             |
     *  (0) m ---[mpo]--- (1) m
     *            |
     *        (3) d
     *
     * is shuffled into
     *
     *          (1) d
     *             |
     *  (0) m ---[mpo]--- (3) m
     *            |
     *        (2) d
     *
     * and reshaped like
     *
     * d (0) ---[mpo]--- (123) ddm
     *
     * and subjected to the typical SVD so mpo = USV. This is then reshaped back into
     *
     *                                          d
     *                                          |
     *   m---[U]---m'   m'---[S]---'m   m' ---[mpo]---  m
     *                                          |
     *                                          d
     *
     * where hopefully m' < m and the transfer matrix T = US is multiplied onto the mpo on the left later.
     *
     * To stablize the compression, it is useful to insert avgS *  1/avgS, where one factor is put into V and the other into S.
     *
     *
     */

    auto dim0    = mpo.dimension(0);
    auto dim1    = mpo.dimension(2);
    auto dim2    = mpo.dimension(3);
    auto dim3    = mpo.dimension(1);
    auto dim_ddm = dim1 * dim2 * dim3;

    Eigen::Tensor<Scalar, 2> mpo_rank2 = mpo.shuffle(tenx::array4{0, 2, 3, 1}).reshape(tenx::array2{dim0, dim_ddm});
    auto [U, S, VT]                     = do_svd_ptr(mpo_rank2.data(), mpo_rank2.dimension(0), mpo_rank2.dimension(1), svd_cfg);

    auto Smin = S.real().minCoeff();
    auto Smax = S.real().maxCoeff();

    // Stabilize by inserting avgS *  1/avgS
    auto avgS = num::next_power_of_two<double>(S.head(S.nonZeros()).real().mean()); // Nearest power of two larger than S.mean()
    if(avgS > 1) {
        S /= avgS;
        std::tie(rank, truncation_error) = get_rank_from_truncation_error(S.real().head(S.real().nonZeros()));
        U  = U.leftCols(rank).eval();
        S  = S.head(rank).eval();
        VT = VT.topRows(rank).eval();
        VT *= avgS;

    }
    // TRY PRINTING S before rescaling
    fmt::print("S r2l min {:.5e} | max {:.5e} -->  min {:.5e} | max {:.5e} | size {}\n", Smin, Smax, S.real().minCoeff(), S.real().maxCoeff(), S.size());
    // rank = dropfilter(U, S, V, svd_cfg.truncation_limit.value_or(1e-16), 8);
    U    = U * S.asDiagonal();

    /* clang-format off */
    return std::make_tuple(
            tenx::TensorMap(U).template cast<Scalar>(),
            tenx::TensorMap(VT).reshape(tenx::array4{rank, dim1, dim2, dim3}).shuffle(
                    tenx::array4{0, 3, 1, 2}).template cast<Scalar>());
    /* clang-format on */
}

template std::tuple<Eigen::Tensor<real, 2>, Eigen::Tensor<real, 4>> svd::solver::split_mpo_r2l(const Eigen::Tensor<real, 4> &mpo, const svd::config &svd_cfg);

template std::tuple<Eigen::Tensor<cplx, 2>, Eigen::Tensor<cplx, 4>> svd::solver::split_mpo_r2l(const Eigen::Tensor<cplx, 4> &mpo, const svd::config &svd_cfg);
