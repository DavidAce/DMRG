#include <Eigen/QR>
#include <math/svd.h>
#include <tid/tid.h>

std::optional<long long> svd::solver::count = 0;

svd::solver::solver() {
    setLogLevel(2);
    if(not count) count = 0;
}
void svd::solver::copy_settings(const svd::settings &svd_settings) {
    if(svd_settings.threshold) threshold = svd_settings.threshold.value();
    if(svd_settings.switchsize_bdc) switchsize_bdc = svd_settings.switchsize_bdc.value();
    if(svd_settings.loglevel) setLogLevel(svd_settings.loglevel.value());
    if(svd_settings.use_bdc) use_bdc = svd_settings.use_bdc.value();
    if(svd_settings.svd_lib) svd_lib = svd_settings.svd_lib.value();
    if(svd_settings.save_fail) save_fail = svd_settings.save_fail.value();
}

svd::solver::solver(const svd::settings &svd_settings) : solver() { copy_settings(svd_settings); }

svd::solver::solver(std::optional<svd::settings> svd_settings) : solver() {
    if(svd_settings) copy_settings(svd_settings.value());
}

void svd::solver::setLogLevel(size_t logLevel) {
    if(not svd::log)
        tools::Logger::setLogger(svd::log, "svd", logLevel, true);
    else
        tools::Logger::setLogLevel(svd::log, logLevel);
    svd::log->set_pattern("[%Y-%m-%d %H:%M:%S.%e][%n]%^[%=8l]%$ %v");
}

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
    svd::solver::do_svd_ptr(const Scalar *mat_ptr, long rows, long cols, std::optional<long> rank_max) {
    auto t_svd = tid::tic_scope("svd");
    if(not rank_max.has_value())
        rank_max = std::min(rows, cols);
    else
        rank_max = std::min({rank_max.value(), rows, cols});
    auto minrc    = std::min(rows, cols);
    bool use_rsvd = num::cmp_greater(minrc, rank_max.value() * 10) and num::cmp_greater(minrc, switchsize_rnd);

    switch(svd_lib) {
        case SVDLib::lapacke: {
            try {
                if(use_rsvd) // Make sure the problem is large enough so that it pays to use rsvd
                    return do_svd_rsvd(mat_ptr, rows, cols, rank_max);
                else
                    return do_svd_lapacke(mat_ptr, rows, cols, rank_max);
            } catch(const std::exception &ex) {
                svd::log->warn("Lapacke failed to perform SVD: {} | Trying Eigen", std::string_view(ex.what()));
                return do_svd_eigen(mat_ptr, rows, cols, rank_max);
            }
            break;
        }
        case SVDLib::eigen: {
            try {
                if(use_rsvd) // Make sure the problem is large enough so that it pays to use rsvd
                    return do_svd_rsvd(mat_ptr, rows, cols, rank_max);
                else
                    return do_svd_eigen(mat_ptr, rows, cols, rank_max);
            } catch(const std::exception &ex) {
                svd::log->warn("Eigen failed to perform SVD: {} | Trying Lapacke", ex.what());
                return do_svd_lapacke(mat_ptr, rows, cols, rank_max);
            }
            break;
        }
        case SVDLib::rsvd: {
            try {
                return do_svd_rsvd(mat_ptr, rows, cols, rank_max);
            } catch(const std::exception &ex) {
                svd::log->warn("Rsvd failed to perform SVD: {} | Trying Lapacke", ex.what());
                return do_svd_lapacke(mat_ptr, rows, cols, rank_max);
            }
            break;
        }
        default: throw std::logic_error("Unrecognized svd library");
    }
    throw std::logic_error("Unrecognized svd library");
}

template<typename Scalar>
void svd::solver::print_matrix(const Scalar *mat_ptr, long rows, long cols, long dec) {
    auto A = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(mat_ptr, rows, cols);
    svd::log->warn("Print matrix of dimensions {}x{}\n", rows, cols);
    for(long r = 0; r < A.rows(); r++) {
        if constexpr(std::is_same_v<Scalar, std::complex<double>>)
            for(long c = 0; c < A.cols(); c++) fmt::print("({1:.{0}f},{2:+.{0}f}) ", dec, std::real(A(r, c)), std::imag(A(r, c)));
        else
            for(long c = 0; c < A.cols(); c++) fmt::print("{1:.{0}f} ", dec, A(r, c));
        fmt::print("\n");
    }
}
template<typename Scalar>
void svd::solver::print_vector(const Scalar *vec_ptr, long size, long dec) {
    auto V = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(vec_ptr, size);
    svd::log->warn("Print matrix of size {}\n", size);
    if constexpr(std::is_same_v<Scalar, std::complex<double>>)
        for(long i = 0; i < V.size(); i++) fmt::print("({1:.{0}f},{2:+.{0}f})\n", dec, std::real(V[i]), std::imag(V[i]));
    else
        for(long i = 0; i < V.size(); i++) fmt::print("{1:.{0}f}\n", dec, V[i]);
}

using cplx = std::complex<double>;
using real = double;

template void svd::solver::print_matrix<real>(const real *vec_ptr, long rows, long cols, long dec);
template void svd::solver::print_matrix<cplx>(const cplx *vec_ptr, long rows, long cols, long dec);
template void svd::solver::print_vector<real>(const real *vec_ptr, long size, long dec);
template void svd::solver::print_vector<cplx>(const cplx *vec_ptr, long size, long dec);

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'double'
template std::tuple<svd::solver::MatrixType<double>, svd::solver::VectorType<double>, svd::solver::MatrixType<double>, long>
    svd::solver::do_svd_ptr(const double *, long, long, std::optional<long>);

//! \relates svd::class_SVD
//! \brief force instantiation of do_svd for type 'std::complex<double>'
template std::tuple<svd::solver::MatrixType<cplx>, svd::solver::VectorType<cplx>, svd::solver::MatrixType<cplx>, long>
    svd::solver::do_svd_ptr(const cplx *, long, long, std::optional<long>);

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