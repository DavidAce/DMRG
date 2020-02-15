
#include "matrix_product_dense.h"
#include <optional>
#include <Eigen/Core>
#include <Eigen/LU>

template<typename T> using MatrixType  = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
template<typename T> using VectorType  = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;
template<typename T> using VectorTypeT = Eigen::Matrix<T, 1, Eigen::Dynamic, Eigen::RowMajor>;

namespace dense_lu {
    std::optional<Eigen::PartialPivLU<MatrixType<double>>>               lu_real       = std::nullopt;
    std::optional<Eigen::PartialPivLU<MatrixType<std::complex<double>>>> lu_cplx       = std::nullopt;

    void reset(){
        lu_real.reset();
        lu_cplx.reset();
    }
    template<typename Scalar>
    void init(){
        if constexpr (std::is_same_v<Scalar,double>)
            dense_lu::lu_real = Eigen::PartialPivLU<MatrixType<Scalar>>();
        if constexpr (std::is_same_v<Scalar,std::complex<double>>)
            dense_lu::lu_cplx = Eigen::PartialPivLU<MatrixType<Scalar>>();
    }
}


template<typename Scalar>
DenseMatrixProduct<Scalar>::~DenseMatrixProduct(){
    dense_lu::reset();
}



// Pointer to data constructor, copies the matrix into an internal Eigen matrix.
template<typename Scalar>
DenseMatrixProduct<Scalar>::DenseMatrixProduct(const Scalar *const A_, const int L_,const bool copy_data, const eigutils::eigSetting::Form form_, const eigutils::eigSetting::Side side_)
    : A_ptr(A_), L(L_), form(form_), side(side_) {
    if (copy_data){
        A_stl.resize(L*L);
        std::copy(A_ptr,A_ptr + L*L, A_stl.begin());
        A_ptr = A_stl.data();
    }
    dense_lu::init<Scalar>();
    init_profiling();
}


template<typename Scalar>
void DenseMatrixProduct<Scalar>::print() const {
    Eigen::Map<const MatrixType<Scalar>> A_matrix (A_ptr, L, L);
    std::cout << "A_matrix: \n" << A_matrix << std::endl;
}

// Function definitions

template<typename Scalar>
void DenseMatrixProduct<Scalar>::FactorOP()

/*  Partial pivot LU decomposition
 *  Factors P(A-sigma*I) = LU
 */
{
    if(readyFactorOp) { return; }
    std::cout << "Starting LU \n";
    Eigen::Map<const MatrixType<Scalar>> A_matrix (A_ptr, L, L);

    t_factorOp.tic();
    assert(readyShift and "Shift value sigma has not been set.");

    if constexpr(std::is_same_v<Scalar, double>) {
        dense_lu::lu_real = Eigen::PartialPivLU<MatrixType<Scalar>>();
        dense_lu::lu_real.value().compute(A_matrix - sigmaR * Eigen::MatrixXd::Identity(L, L));
    }
    if constexpr(std::is_same_v<Scalar, std::complex<double>>) {
        Scalar sigma = std::complex<double>(sigmaR, sigmaI);
        dense_lu::lu_cplx = Eigen::PartialPivLU<MatrixType<Scalar>>();
        dense_lu::lu_cplx.value().compute(A_matrix - sigma * Eigen::MatrixXd::Identity(L, L));
    }

    readyFactorOp = true;
    t_factorOp.toc();
    std::cout << "Finished LU \n";
    std::cout << "Time LU Op [ms]: " << std::fixed << std::setprecision(3) << t_factorOp.get_last_time_interval() * 1000 << '\n';
}

template<typename Scalar>
void DenseMatrixProduct<Scalar>::MultOPv(Scalar *x_in_ptr, Scalar *x_out_ptr) {
    using namespace eigutils::eigSetting;
    assert(readyFactorOp and "FactorOp() has not been run yet.");
    switch(side) {
        case Side::R: {
            Eigen::Map<VectorType<Scalar>> x_in(x_in_ptr, L);
            Eigen::Map<VectorType<Scalar>> x_out(x_out_ptr, L);
            if constexpr(std::is_same_v<Scalar, double>)
                x_out.noalias() = dense_lu::lu_real.value().solve(x_in);
            if constexpr(std::is_same_v<Scalar, std::complex<double>>)
                x_out.noalias() = dense_lu::lu_cplx.value().solve(x_in);
            break;
        }
        case Side::L: {
            Eigen::Map<VectorTypeT<Scalar>> x_in(x_in_ptr, L);
            Eigen::Map<VectorTypeT<Scalar>> x_out(x_out_ptr, L);
            if constexpr(std::is_same_v<Scalar, double>)
                x_out.noalias() = x_in * dense_lu::lu_real.value().inverse();
            if constexpr(std::is_same_v<Scalar, std::complex<double>>)
                x_out.noalias() = x_in * dense_lu::lu_cplx.value().inverse();
            break;
        }
    }
    counter++;
}

template<typename Scalar>
void DenseMatrixProduct<Scalar>::MultAx(Scalar *x_in, Scalar *x_out) {
    using namespace eigutils::eigSetting;
    Eigen::Map<const MatrixType<Scalar>> A_matrix (A_ptr, L, L);
    switch(form) {
        case Form::NONSYMMETRIC:
            switch(side) {
                case Side::R: {
                    Eigen::Map<VectorType<Scalar>> x_vec_in(x_in, L);
                    Eigen::Map<VectorType<Scalar>> x_vec_out(x_out, L);
                    x_vec_out.noalias() = A_matrix * x_vec_in;
                    break;
                }
                case Side::L: {
                    Eigen::Map<VectorTypeT<Scalar>> x_vec_in(x_in, L);
                    Eigen::Map<VectorTypeT<Scalar>> x_vec_out(x_out, L);
                    x_vec_out.noalias() = x_vec_in * A_matrix;
                    break;
                }
            }
            break;
        case Form::SYMMETRIC: {
            Eigen::Map<VectorType<Scalar>> x_vec_in(x_in, L);
            Eigen::Map<VectorType<Scalar>> x_vec_out(x_out, L);
            x_vec_out.noalias() = A_matrix.template selfadjointView<Eigen::Upper>() * x_vec_in;
            break;
        }
    }
    counter++;
}

// Explicit instantiations
template class DenseMatrixProduct<double>;
template class DenseMatrixProduct<std::complex<double>>;
