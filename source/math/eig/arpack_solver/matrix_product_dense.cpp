//
// Created by david on 2018-11-16.
//

#include "matrix_product_dense.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <general/class_tic_toc.h>
#include <memory>
#define profile_matrix_product_dense 1

// Function definitions
template<typename Scalar>
using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
template<typename Scalar>
using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
template<typename Scalar>
using VectorTypeT = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;

namespace dense_lu {
    std::optional<Eigen::PartialPivLU<MatrixType<double>>>               lu_real;
    std::optional<Eigen::PartialPivLU<MatrixType<std::complex<double>>>> lu_cplx;
    void                                                                 reset() {
        lu_real.reset();
        lu_cplx.reset();
    }
    template<typename Scalar>
    void init() {
        if constexpr(std::is_same_v<Scalar, double>) { dense_lu::lu_real = Eigen::PartialPivLU<MatrixType<Scalar>>(); }
        if constexpr(std::is_same_v<Scalar, std::complex<double>>) dense_lu::lu_cplx = Eigen::PartialPivLU<MatrixType<Scalar>>();
    }

}

template<typename Scalar>
MatrixProductDense<Scalar>::~MatrixProductDense() {
    dense_lu::reset();
}

// Pointer to data constructor, copies the matrix into an internal Eigen matrix.
template<typename Scalar>
MatrixProductDense<Scalar>::MatrixProductDense(const Scalar *const A_, const long L_, const bool copy_data, const eig::Form form_, const eig::Side side_)
    : A_ptr(A_), L(L_), form(form_), side(side_) {
    if(copy_data) {
        A_stl.resize(static_cast<size_t>(L * L));
        std::copy(A_ptr, A_ptr + static_cast<size_t>(L * L), A_stl.begin());
        A_ptr = A_stl.data();
    }
    dense_lu::init<Scalar>();
    init_profiling();
}

template<typename Scalar>
void MatrixProductDense<Scalar>::FactorOP()

/*  Partial pivot LU decomposition
 *  Factors P(A-sigma*I) = LU
 */
{
    if(readyFactorOp) return; // happens only once
    if(not readyShift) throw std::runtime_error("Cannot FactorOP: Shift value sigma has not been set.");
    Eigen::Map<const MatrixType<Scalar>> A_matrix(A_ptr, L, L);
    t_factorOP->tic();
    if constexpr(std::is_same_v<Scalar, eig::real>) { dense_lu::lu_real.value().compute(A_matrix); }
    if constexpr(std::is_same_v<Scalar, eig::cplx>) { dense_lu::lu_cplx.value().compute(A_matrix); }

    readyFactorOp = true;
    t_factorOP->toc();
}

template<typename Scalar>
void MatrixProductDense<Scalar>::MultOPv(Scalar *x_in_ptr, Scalar *x_out_ptr) {
    assert(readyFactorOp and "FactorOp() has not been run yet.");
    t_multOPv->tic();
    switch(side) {
        case eig::Side::R: {
            Eigen::Map<VectorType<Scalar>> x_in(x_in_ptr, L);
            Eigen::Map<VectorType<Scalar>> x_out(x_out_ptr, L);
            if constexpr(std::is_same_v<Scalar, double>) x_out.noalias() = dense_lu::lu_real.value().solve(x_in);
            if constexpr(std::is_same_v<Scalar, std::complex<double>>) x_out.noalias() = dense_lu::lu_cplx.value().solve(x_in);
            break;
        }
        case eig::Side::L: {
            Eigen::Map<VectorTypeT<Scalar>> x_in(x_in_ptr, L);
            Eigen::Map<VectorTypeT<Scalar>> x_out(x_out_ptr, L);
            if constexpr(std::is_same_v<Scalar, double>) x_out.noalias() = x_in * dense_lu::lu_real.value().inverse();
            if constexpr(std::is_same_v<Scalar, std::complex<double>>) x_out.noalias() = x_in * dense_lu::lu_cplx.value().inverse();
            break;
        }
        case eig::Side::LR: {
            throw std::runtime_error("eigs cannot handle sides L and R simultaneously");
        }
    }
    t_multOPv->toc();
    counter++;
}

template<typename Scalar>
void MatrixProductDense<Scalar>::MultAx(Scalar *x_in, Scalar *x_out) {
    t_multAx->tic();
    Eigen::Map<const MatrixType<Scalar>> A_matrix(A_ptr, L, L);
    switch(form) {
        case eig::Form::NSYM:
            switch(side) {
                case eig::Side::R: {
                    Eigen::Map<VectorType<Scalar>> x_vec_in(x_in, L);
                    Eigen::Map<VectorType<Scalar>> x_vec_out(x_out, L);
                    x_vec_out.noalias() = A_matrix * x_vec_in;
                    break;
                }
                case eig::Side::L: {
                    Eigen::Map<VectorTypeT<Scalar>> x_vec_in(x_in, L);
                    Eigen::Map<VectorTypeT<Scalar>> x_vec_out(x_out, L);
                    x_vec_out.noalias() = x_vec_in * A_matrix;
                    break;
                }
                case eig::Side::LR: {
                    throw std::runtime_error("eigs cannot handle sides L and R simultaneously");
                }
            }
            break;
        case eig::Form::SYMM: {
            Eigen::Map<VectorType<Scalar>> x_vec_in(x_in, L);
            Eigen::Map<VectorType<Scalar>> x_vec_out(x_out, L);
            x_vec_out.noalias() = A_matrix.template selfadjointView<Eigen::Lower>() * x_vec_in;
            break;
        }
    }
    t_multAx->toc();
    counter++;
}

template<typename Scalar>
void MatrixProductDense<Scalar>::print() const {
    Eigen::Map<const MatrixType<Scalar>> A_matrix(A_ptr, L, L);
}

template<typename Scalar>
void MatrixProductDense<Scalar>::set_shift(std::complex<double> sigma_) {
    if(readyShift) { return; }
    sigma = sigma_;
    if(A_stl.empty()) {
        A_stl.resize(static_cast<size_t>(L * L));
        std::copy(A_ptr, A_ptr + static_cast<size_t>(L * L), A_stl.begin());
        A_ptr = A_stl.data();
    }
    Eigen::Map<MatrixType<Scalar>> A_matrix(A_stl.data(), L, L);
    if constexpr(std::is_same_v<Scalar, eig::real>) A_matrix -= Eigen::MatrixXd::Identity(L, L) * std::real(sigma);
    if constexpr(std::is_same_v<Scalar, eig::cplx>) A_matrix -= Eigen::MatrixXd::Identity(L, L) * sigma;
    readyShift = true;
}

template<typename Scalar>
void MatrixProductDense<Scalar>::set_mode(const eig::Form form_) {
    form = form_;
}
template<typename Scalar>
void MatrixProductDense<Scalar>::set_side(const eig::Side side_) {
    side = side_;
}
template<typename Scalar>
const eig::Form &MatrixProductDense<Scalar>::get_form() const {
    return form;
}
template<typename Scalar>
const eig::Side &MatrixProductDense<Scalar>::get_side() const {
    return side;
}

template<typename Scalar>
void MatrixProductDense<Scalar>::init_profiling() {
    t_factorOP = std::make_unique<class_tic_toc>(profile_matrix_product_dense, 5, "Time FactorOp");
    t_multOPv  = std::make_unique<class_tic_toc>(profile_matrix_product_dense, 5, "Time MultOpv");
    t_multAx   = std::make_unique<class_tic_toc>(profile_matrix_product_dense, 5, "Time MultAx");
}

// Explicit instantiations

template class MatrixProductDense<double>;
template class MatrixProductDense<std::complex<double>>;
