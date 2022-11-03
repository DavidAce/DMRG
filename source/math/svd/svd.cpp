#include "math/svd.h"
#include "tid/tid.h"
#include <Eigen/QR>

long long svd::solver::count = 0;

svd::solver::solver() { setLogLevel(2); }

svd::config::config(long rank_max_) : rank_max(rank_max_) {}
svd::config::config(double truncation_lim_) : truncation_lim(truncation_lim_) {}
svd::config::config(long rank_max_, double truncation_lim_) : rank_max(rank_max_), truncation_lim(truncation_lim_) {}
svd::config::config(std::optional<long> rank_max_) : rank_max(rank_max_) {}
svd::config::config(std::optional<double> truncation_lim_) : truncation_lim(truncation_lim_) {}
svd::config::config(std::optional<long> rank_max_, std::optional<double> truncation_lim_) : rank_max(rank_max_), truncation_lim(truncation_lim_) {}

void svd::solver::copy_config(const svd::config &svd_cfg) {
    if(svd_cfg.rank_max) rank_max = svd_cfg.rank_max.value();
    if(svd_cfg.rank_min) rank_min = svd_cfg.rank_min.value();
    if(svd_cfg.truncation_lim) truncation_lim = svd_cfg.truncation_lim.value();
    if(svd_cfg.switchsize_bdc) switchsize_bdc = svd_cfg.switchsize_bdc.value();
    if(svd_cfg.loglevel) setLogLevel(svd_cfg.loglevel.value());
    if(svd_cfg.use_bdc) use_bdc = svd_cfg.use_bdc.value();
    if(svd_cfg.svd_lib) svd_lib = svd_cfg.svd_lib.value();
    if(svd_cfg.save_fail) save_fail = svd_cfg.save_fail.value();
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

void svd::solver::setLogLevel(size_t logLevel) {
    if(not svd::log)
        tools::Logger::setLogger(svd::log, "svd", logLevel, true);
    else
        tools::Logger::setLogLevel(svd::log, logLevel);
    svd::log->set_pattern("[%Y-%m-%d %H:%M:%S.%e][%n]%^[%=8l]%$ %v");
}

double    svd::solver::get_truncation_error() const { return truncation_error; }
long      svd::solver::get_rank() const { return rank; }
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
std::tuple<svd::solver::MatrixType<Scalar>, svd::solver::VectorType<Scalar>, svd::solver::MatrixType<Scalar>>
    svd::solver::do_svd_ptr(const Scalar *mat_ptr, long rows, long cols, const svd::config &svd_cfg) {
#if defined(DMRG_ENABLE_RSVD)
    constexpr bool has_rsvd = true;
#else
    constexpr bool has_rsvd = false;
#endif

    auto t_svd = tid::tic_scope("svd", tid::level::highest);

    copy_config(svd_cfg);
    auto sizeS    = std::min(rows, cols);
    long rank_lim = rank_max > 0 ? std::min(sizeS, rank_max) : sizeS;
    bool use_rsvd = has_rsvd and num::cmp_greater(sizeS, rank_lim * 10) and num::cmp_greater(sizeS, switchsize_rnd);
#pragma omp atomic
    count++;
    switch(svd_lib) {
        case svd::lib::lapacke: {
            try {
                if(use_rsvd) // Make sure the problem is large enough so that it pays to use rsvd
                    return do_svd_rsvd(mat_ptr, rows, cols);
                else
                    return do_svd_lapacke(mat_ptr, rows, cols);
            } catch(const std::exception &ex) {
                svd::log->warn("Lapacke failed to perform SVD: {} | Trying Eigen", std::string_view(ex.what()));
                return do_svd_eigen(mat_ptr, rows, cols);
            }
            break;
        }
        case svd::lib::eigen: {
            try {
                if(use_rsvd) // Make sure the problem is large enough so that it pays to use rsvd
                    return do_svd_rsvd(mat_ptr, rows, cols);
                else
                    return do_svd_eigen(mat_ptr, rows, cols);
            } catch(const std::exception &ex) {
                svd::log->warn("Eigen failed to perform SVD: {} | Trying Lapacke", ex.what());
                return do_svd_lapacke(mat_ptr, rows, cols);
            }
            break;
        }
        case svd::lib::rsvd: {
            try {
                return do_svd_rsvd(mat_ptr, rows, cols);
            } catch(const std::exception &ex) {
                svd::log->warn("Rsvd failed to perform SVD: {} | Trying Lapacke", ex.what());
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
template std::tuple<svd::solver::MatrixType<real>, svd::solver::VectorType<real>, svd::solver::MatrixType<real>>
    svd::solver::do_svd_ptr(const real *, long, long, const svd::config &);

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'std::complex<double>'
template std::tuple<svd::solver::MatrixType<cplx>, svd::solver::VectorType<cplx>, svd::solver::MatrixType<cplx>>
    svd::solver::do_svd_ptr(const cplx *, long, long, const svd::config &);

template<typename Scalar>
void svd::solver::print_matrix([[maybe_unused]] const Scalar *mat_ptr, [[maybe_unused]] long rows, [[maybe_unused]] long cols,
                               [[maybe_unused]] long dec) const {
#if !defined(NDEBUG)
    auto A = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(mat_ptr, rows, cols);
    svd::log->warn("Print matrix of dimensions {}x{}\n", rows, cols);
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
void svd::solver::print_vector([[maybe_unused]] const Scalar *vec_ptr, [[maybe_unused]] long size, [[maybe_unused]] long dec) const {
#if !defined(NDEBUG)
    auto V = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(vec_ptr, size);
    svd::log->warn("Print matrix of size {}\n", size);
    if constexpr(std::is_same_v<Scalar, std::complex<double>>)
        for(long i = 0; i < V.size(); i++) fmt::print("({1:.{0}f},{2:+.{0}f})\n", dec, std::real(V[i]), std::imag(V[i]));
    else
        for(long i = 0; i < V.size(); i++) fmt::print("{1:.{0}f}\n", dec, V[i]);
#endif
}

template void svd::solver::print_matrix<real>(const real *vec_ptr, long rows, long cols, long dec) const;
template void svd::solver::print_matrix<cplx>(const cplx *vec_ptr, long rows, long cols, long dec) const;
template void svd::solver::print_vector<real>(const real *vec_ptr, long size, long dec) const;
template void svd::solver::print_vector<cplx>(const cplx *vec_ptr, long size, long dec) const;

template<typename Scalar>
Eigen::Tensor<Scalar, 2> svd::solver::pseudo_inverse(const Eigen::Tensor<Scalar, 2> &tensor) {
    auto t_psinv = tid::tic_token("psinv");
    if(tensor.dimension(0) <= 0) { throw std::runtime_error("pseudo_inverse error: Dimension is zero: tensor.dimension(0)"); }
    if(tensor.dimension(1) <= 0) { throw std::runtime_error("pseudo_inverse error: Dimension is zero: tensor.dimension(1)"); }
    Eigen::Map<const MatrixType<Scalar>> mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    return tenx::TensorCast(mat.completeOrthogonalDecomposition().pseudoInverse());
}

//! \relates svd::class_SVD
//! \brief force instantiation of pseudo_inverse for type 'double'
template Eigen::Tensor<double, 2> svd::solver::pseudo_inverse(const Eigen::Tensor<double, 2> &tensor);
//! \relates svd::class_SVD
//! \brief force instantiation of pseudo_inverse for type 'std::complex<double>'
template Eigen::Tensor<cplx, 2> svd::solver::pseudo_inverse(const Eigen::Tensor<cplx, 2> &tensor);

// template<typename Scalar>
std::pair<long, double> svd::solver::get_rank_from_truncation_error(const VectorType<double> &S) const {
    assert(std::abs(S.norm() - 1.0) < 1e-10); // make sure this is normalized
    VectorType<double> truncation_errors(S.size() + 1);
    for(long s = 0; s <= S.size(); s++) { truncation_errors[s] = S.bottomRows(S.size() - s).norm(); } // Last one should be zero, i.e. no truncation
    auto rank_    = (truncation_errors.array() >= truncation_lim).count();
    auto rank_lim = S.size();
    if(rank_max > 0) rank_lim = std::min(S.size(), rank_max);
    rank_ = std::min(rank_, rank_lim);
    if(rank_min > 0) rank_ = std::max(rank_, std::min(S.size(), rank_min)); // Make sure we don't overtruncate in some cases (e.g. when stashing)

    //    tools::log->info("Size {} | Rank {} | Rank limit {} | truncation error limit {:8.5e} | error {:8.5e}", S.size(), rank_, rank_lim,
    //                     truncation_lim, truncation_errors[rank_]);
    //    tools::log->info("Size {} | Rank {} | Rank limit {} | truncation error limit {:8.5e} | error {:8.5e} truncation errors: {:8.5e}", S.size(), rank_,
    //    rank_lim,
    //                     truncation_lim, truncation_errors[rank_], fmt::join(truncation_errors, ", "));
    if(rank_ <= 0) {
        svd::log->error("Size {} | Rank {} | Rank limit {} | truncation error limit {:8.2e} | error {:8.2e} truncation errors: {:8.2e}", S.size(), rank_,
                        rank_lim, truncation_lim, truncation_errors[rank_], fmt::join(truncation_errors, ", "));
        throw std::logic_error("rank <= 0");
    }
    return {rank_, truncation_errors[rank_]};
}