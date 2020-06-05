//
// Created by david on 2018-11-16.
//

#include "matrix_product_stl.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <general/class_tic_toc.h>
#define profile_matrix_product_dense 1

// Function definitions
template <typename Scalar>
using MatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;

template <typename Scalar>
using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;

template <typename Scalar>
using VectorTypeT = Eigen::Matrix<Scalar,1,Eigen::Dynamic>;


namespace stl_lu{
    std::optional<Eigen::PartialPivLU<MatrixType<double>>               >lu_real;
    std::optional<Eigen::PartialPivLU<MatrixType<std::complex<double>>> >lu_cplx;
    void reset(){
        lu_real.reset();
        lu_cplx.reset();
    }
    template<typename Scalar>
    void init(){
        if constexpr (std::is_same_v<Scalar,double>){
            stl_lu::lu_real   = Eigen::PartialPivLU<MatrixType<Scalar>>();
        }
        if constexpr (std::is_same_v<Scalar,std::complex<double>>)
            stl_lu::lu_cplx = Eigen::PartialPivLU<MatrixType<Scalar>>();
    }

}

template<typename Scalar>
void StlMatrixProduct<Scalar>::init_profiling(){
    t_factorOP = std::make_unique<class_tic_toc>(profile_matrix_product_dense, 5,"Time FactorOp");
    t_multOPv = std::make_unique<class_tic_toc>(profile_matrix_product_dense, 5,"Time MultOpv");
    t_multAx = std::make_unique<class_tic_toc>(profile_matrix_product_dense, 5,"Time MultAx");
}


template<typename Scalar>
StlMatrixProduct<Scalar>::~StlMatrixProduct(){
    stl_lu::reset();
}

// Pointer to data constructor, copies the matrix into an internal Eigen matrix.
template<typename Scalar>
StlMatrixProduct<Scalar>::StlMatrixProduct(
    const Scalar * const A_,
    const int L_,
    const bool copy_data,
    const eigutils::eigSetting::Form form_,
    const eigutils::eigSetting::Side side_):
    A_ptr(A_) ,L(L_), form(form_), side(side_)
{
    if (copy_data){
        A_stl.resize(static_cast<size_t>(L*L));
        std::copy(A_ptr,A_ptr + static_cast<size_t>(L*L), A_stl.begin());
        A_ptr = A_stl.data();
    }
    stl_lu::init<Scalar>();
    init_profiling();
}

template<typename Scalar>
void StlMatrixProduct<Scalar>::print() const {
    Eigen::Map<const MatrixType<Scalar>> A_matrix (A_ptr,L,L);
    std::cout << "A_matrix: \n" << A_matrix << std::endl;
}


template<typename Scalar>
void StlMatrixProduct<Scalar>::FactorOP()

/*  Partial pivot LU decomposition
 *  Factors P(A-sigma*I) = LU
 */
{
    if(readyFactorOp){return;}
//    lu_real_ptr = std::make_shared<LU_REAL>( LU_REAL() );
    Eigen::Map<const MatrixType<Scalar>> A_matrix (A_ptr,L,L);
    t_factorOP->tic();
    assert(readyShift and "Shift value sigma has not been set.");
    if constexpr(std::is_same_v<Scalar,double>){
        stl_lu::lu_real.value().compute(A_matrix - sigmaR * Eigen::MatrixXd::Identity(L,L));
    }
    if constexpr(std::is_same_v<Scalar,std::complex<double>>){
        Scalar sigma = std::complex<double>(sigmaR,sigmaI);
        stl_lu::lu_cplx.value().compute(A_matrix - sigma * Eigen::MatrixXd::Identity(L,L));
    }

    readyFactorOp = true;
    t_factorOP->toc();
//    std::cout << "Time Factor Op [ms]: " << std::fixed << std::setprecision(3) << t_factorOP.get_last_time_interval() * 1000 << '\n';
}




template<typename Scalar>
void StlMatrixProduct<Scalar>::MultOPv(Scalar* x_in_ptr, Scalar* x_out_ptr) {
    using namespace eigutils::eigSetting;
    assert(readyFactorOp and "FactorOp() has not been run yet.");
    t_multOPv->tic();
    switch (side){
        case Side::R: {
            Eigen::Map<VectorType<Scalar>>       x_in    (x_in_ptr,L);
            Eigen::Map<VectorType<Scalar>>       x_out   (x_out_ptr,L);
            if constexpr(std::is_same_v<Scalar,double>)
                x_out.noalias() = stl_lu::lu_real.value().solve(x_in);
            if constexpr(std::is_same_v<Scalar,std::complex<double>>)
                x_out.noalias() = stl_lu::lu_cplx.value().solve(x_in);
            break;
        }
        case Side::L: {
            Eigen::Map<VectorTypeT<Scalar>>       x_in    (x_in_ptr,L);
            Eigen::Map<VectorTypeT<Scalar>>       x_out   (x_out_ptr,L);
            if constexpr(std::is_same_v<Scalar,double>)
                x_out.noalias() = x_in *stl_lu::lu_real.value().inverse();
            if constexpr(std::is_same_v<Scalar,std::complex<double>>)
                x_out.noalias() = x_in *stl_lu::lu_cplx.value().inverse();
            break;
        }
    }
    t_multOPv->toc();
    counter++;
}




template<typename Scalar>
void StlMatrixProduct<Scalar>::MultAx(Scalar* x_in, Scalar* x_out) {
    using namespace eigutils::eigSetting;
    t_multAx->tic();
    Eigen::Map<const MatrixType<Scalar>> A_matrix (A_ptr,L,L);
    switch (form){
        case Form::NONSYMMETRIC:
            switch (side) {
                case Side::R: {
                    Eigen::Map<VectorType<Scalar>> x_vec_in (x_in,  L);
                    Eigen::Map<VectorType<Scalar>> x_vec_out(x_out, L);
                    x_vec_out.noalias() = A_matrix * x_vec_in ;
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
            x_vec_out.noalias() = A_matrix.template selfadjointView<Eigen::Lower>() * x_vec_in;
            break;
        }
    }
    t_multAx->toc();
    counter++;
}






// Explicit instantiations

template class StlMatrixProduct<double>;
template class StlMatrixProduct<std::complex<double>>;
