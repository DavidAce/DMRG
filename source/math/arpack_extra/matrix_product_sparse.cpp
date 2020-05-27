
#include "matrix_product_sparse.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SparseLU>

template<typename T>
using SparseMatrixType = Eigen::SparseMatrix<T>;
template<typename T>
using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
template<typename T>
using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;
template<typename T>
using VectorTypeT = Eigen::Matrix<T, 1, Eigen::Dynamic, Eigen::RowMajor>;

namespace sparse {
    std::optional<Eigen::SparseMatrix<double>>                             A_real_sparse  = std::nullopt;
    std::optional<Eigen::SparseMatrix<std::complex<double>>>               A_cplx_sparse  = std::nullopt;
    std::optional<Eigen::SparseLU<SparseMatrixType<double>>>               lu_real_sparse = {};
    std::optional<Eigen::SparseLU<SparseMatrixType<std::complex<double>>>> lu_cplx_sparse = {};
    std::optional<Eigen::PartialPivLU<MatrixType<double>>>                 lu_real_dense  = std::nullopt;
    std::optional<Eigen::PartialPivLU<MatrixType<std::complex<double>>>>   lu_cplx_dense  = std::nullopt;

    void reset() {
        A_real_sparse.reset();
        A_cplx_sparse.reset();
        lu_real_sparse.reset();
        lu_cplx_sparse.reset();
        lu_real_dense.reset();
        lu_cplx_dense.reset();
    }

}

template<typename Scalar, bool sparseLU>
SparseMatrixProduct<Scalar, sparseLU>::~SparseMatrixProduct(){
    sparse::reset();
}



template<typename Scalar, bool sparseLU>
SparseMatrixProduct<Scalar, sparseLU>::SparseMatrixProduct(const Scalar *A_, const int L_, const bool copy_data, const eigutils::eigSetting::Form form_, const eigutils::eigSetting::Side side_)
    : A_ptr(A_), L(L_), form(form_), side(side_) {
    if (copy_data){
        A_stl.resize(static_cast<size_t>(L*L));
        std::copy(A_ptr,A_ptr + static_cast<size_t>(L*L), A_stl.begin());
        A_ptr = A_stl.data();
    }


    if constexpr(sparseLU) {
        Eigen::Map<const MatrixType<Scalar>> A_matrix(A_ptr, L, L);
        if constexpr(std::is_same_v<Scalar, double>) {
            sparse::A_real_sparse = A_matrix.sparseView();
            sparse::A_real_sparse.value().makeCompressed();
        }
        if constexpr(std::is_same_v<Scalar, std::complex<double>>) {
            sparse::A_cplx_sparse = A_matrix.sparseView();
            sparse::A_cplx_sparse.value().makeCompressed();
        }

    }
    init_profiling();
}

// Function definitions

template<typename Scalar, bool sparseLU>
void SparseMatrixProduct<Scalar, sparseLU>::print() const {
    Eigen::Map<const MatrixType<Scalar>> A_matrix(A_ptr, L, L);
    std::cout << "A_matrix: \n" << A_matrix << std::endl;
}

template<typename Scalar, bool sparseLU>
void SparseMatrixProduct<Scalar, sparseLU>::FactorOP()

/*  Sparse decomposition
 *  Factors P(A-sigma*I) = LU
 */

{
    if(readyFactorOp) { return; }
    assert(readyShift and "Shift value sigma has not been set.");
    t_factorOp.tic();
    Eigen::Map<const MatrixType<Scalar>> A_matrix(A_ptr, L, L);

    Scalar sigma;
    if constexpr(std::is_same_v<Scalar, double>) sigma = sigmaR;
    if constexpr(std::is_same_v<Scalar, std::complex<double>>) sigma = std::complex<double>(sigmaR, sigmaI);
    // Real
    if constexpr(std::is_same_v<Scalar, double> and not sparseLU) {
        sparse::lu_real_dense = Eigen::PartialPivLU<MatrixType<Scalar>>();
        sparse::lu_real_dense.value().compute(A_matrix - sigma * Eigen::MatrixXd::Identity(L, L));
    }
    if constexpr(std::is_same_v<Scalar, double> and sparseLU) {
        //        container::lu_real_sparse = Eigen::SparseLU<SparseMatrixType<Scalar>>();
        sparse::lu_real_sparse.value().compute(sparse::A_real_sparse.value() - sigma * Eigen::MatrixXd::Identity(L, L));
    }
    // Complex
    if constexpr(std::is_same_v<Scalar, std::complex<double>> and not sparseLU) {
        sparse::lu_cplx_dense = Eigen::PartialPivLU<MatrixType<Scalar>>();
        sparse::lu_cplx_dense.value().compute(A_matrix - sigma * Eigen::MatrixXd::Identity(L, L));
    }
    if constexpr(std::is_same_v<Scalar, std::complex<double>> and sparseLU) {
        //        container::lu_cplx_sparse = Eigen::SparseLU<SparseMatrixType<Scalar>>();
        sparse::lu_cplx_sparse.value().compute(sparse::A_cplx_sparse.value() - sigma * Eigen::MatrixXcd::Identity(L, L));
    }

    t_factorOp.toc();
    readyFactorOp = true;
    std::cout << "Time Factor Op [ms]: " << std::fixed << std::setprecision(3) << t_factorOp.get_last_time_interval() * 1000 << '\n';
}

template<typename Scalar, bool sparseLU>
void SparseMatrixProduct<Scalar, sparseLU>::MultOPv(Scalar *x_in_ptr, Scalar *x_out_ptr) {
    using namespace eigutils::eigSetting;
    assert(readyFactorOp and "FactorOp() has not been run yet.");
    using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    Eigen::Map<VectorType> x_in(x_in_ptr, L);
    Eigen::Map<VectorType> x_out(x_out_ptr, L);

    switch(side) {
        case Side::R: {
            if constexpr(std::is_same_v<Scalar, double> and not sparseLU)
                x_out.noalias() = sparse::lu_real_dense.value().solve(x_in);
            else if constexpr(std::is_same_v<Scalar, std::complex<double>> and not sparseLU)
                x_out.noalias() = sparse::lu_cplx_dense.value().solve(x_in);
            else if constexpr(std::is_same_v<Scalar, double> and sparseLU)
                x_out.noalias() = sparse::lu_real_sparse.value().solve(x_in);
            else if constexpr(std::is_same_v<Scalar, std::complex<double>> and sparseLU)
                x_out.noalias() = sparse::lu_cplx_sparse.value().solve(x_in);
            break;
        }
        case Side::L: {
            if constexpr(std::is_same_v<Scalar, double> and not sparseLU)
                x_out.noalias() = x_in * sparse::lu_real_dense.value().inverse();
            else if constexpr(std::is_same_v<Scalar, std::complex<double>> and not sparseLU)
                x_out.noalias() = x_in * sparse::lu_cplx_dense.value().inverse();
            else {
                throw std::runtime_error("Left sided sparse shift invert hasn't been implemented yet...");
            }
            break;
        }
    }
    counter++;
}

template<typename Scalar, bool sparseLU>
void SparseMatrixProduct<Scalar, sparseLU>::MultAx(Scalar *x_in, Scalar *x_out) {
    using namespace eigutils::eigSetting;
    Eigen::Map<const MatrixType<Scalar>> A_matrix(A_ptr, L, L);
    switch(form) {
        case Form::NONSYMMETRIC:
            switch(side) {
                case Side::R: {
                    using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
                    Eigen::Map<VectorType> x_vec_in(x_in, L);
                    Eigen::Map<VectorType> x_vec_out(x_out, L);
                    if constexpr(not sparseLU) x_vec_out.noalias() = A_matrix * x_vec_in;
                    else if constexpr(std::is_same_v<Scalar,double> and sparseLU) x_vec_out.noalias() = sparse::A_real_sparse.value() * x_vec_in;
                    else if constexpr(std::is_same_v<Scalar,std::complex<double>> and sparseLU) x_vec_out.noalias() = sparse::A_cplx_sparse.value() * x_vec_in;
                    break;
                }
                case Side::L: {
                    using VectorTypeT = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;
                    Eigen::Map<VectorTypeT> x_vec_in(x_in, L);
                    Eigen::Map<VectorTypeT> x_vec_out(x_out, L);
                    if constexpr(not sparseLU) x_vec_out.noalias() = x_vec_in * A_matrix;
                    else if constexpr(std::is_same_v<Scalar,double> and sparseLU) x_vec_out.noalias() = x_vec_in * sparse::A_real_sparse.value();
                    else if constexpr(std::is_same_v<Scalar,std::complex<double>> and sparseLU) x_vec_out.noalias() = x_vec_in * sparse::A_cplx_sparse.value();
                    break;
                }
            }
            break;
        case Form::SYMMETRIC: {
            using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
            Eigen::Map<VectorType> x_vec_in(x_in, L);
            Eigen::Map<VectorType> x_vec_out(x_out, L);
            if constexpr(not sparseLU) x_vec_out.noalias() = A_matrix.template selfadjointView<Eigen::Upper>() * x_vec_in;
            if constexpr(std::is_same_v<Scalar,double> and sparseLU) x_vec_out.noalias() = sparse::A_real_sparse.value().template selfadjointView<Eigen::Upper>() * x_vec_in;
            if constexpr(std::is_same_v<Scalar,std::complex<double>> and sparseLU) x_vec_out.noalias() = sparse::A_cplx_sparse.value().template selfadjointView<Eigen::Upper>() * x_vec_in;
            break;
        }
    }
    counter++;
}

// Explicit instantiations
template class SparseMatrixProduct<double, true>;
template class SparseMatrixProduct<double, false>;
template class SparseMatrixProduct<std::complex<double>, true>;
template class SparseMatrixProduct<std::complex<double>, false>;
