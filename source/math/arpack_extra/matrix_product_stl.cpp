//
// Created by david on 2018-11-16.
//


#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <memory>
#include "matrix_product_stl.h"

// Function definitions
template <typename Scalar>
using MatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;

template <typename Scalar>
using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;

template <typename Scalar>
using VectorTypeT = Eigen::Matrix<Scalar,1,Eigen::Dynamic>;

static Eigen::PartialPivLU<MatrixType<double>>               lu_real;
static Eigen::PartialPivLU<MatrixType<std::complex<double>>> lu_cplx;


template<typename Scalar>
StlMatrixProduct<Scalar>::~StlMatrixProduct(){
    lu_real = Eigen::PartialPivLU<MatrixType<double>>              ();
    lu_cplx = Eigen::PartialPivLU<MatrixType<std::complex<double>>>();
//    A_stl.clear();
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
    t_factorOp.tic();
    assert(readyShift and "Shift value sigma has not been set.");
    if constexpr(std::is_same<Scalar,double>::value)
    {
        lu_real.compute(A_matrix - sigmaR * Eigen::MatrixXd::Identity(L,L));
    }
    else
    {
        Scalar sigma = std::complex<double>(sigmaR,sigmaI);
        lu_cplx.compute(A_matrix - sigma * Eigen::MatrixXd::Identity(L,L));
    }

    readyFactorOp = true;
    t_factorOp.toc();
//    std::cout << "Time Factor Op [ms]: " << std::fixed << std::setprecision(3) << t_factorOp.get_last_time_interval() * 1000 << '\n';
}




template<typename Scalar>
void StlMatrixProduct<Scalar>::MultOPv(Scalar* x_in_ptr, Scalar* x_out_ptr) {
    using namespace eigutils::eigSetting;
    assert(readyFactorOp and "FactorOp() has not been run yet.");
    switch (side){
        case Side::R: {
            Eigen::Map<VectorType<Scalar>>       x_in    (x_in_ptr,L);
            Eigen::Map<VectorType<Scalar>>       x_out   (x_out_ptr,L);
            if constexpr(std::is_same <Scalar,double>::value){
                x_out.noalias() = lu_real.solve(x_in);
            }else{
                x_out.noalias() = lu_cplx.solve(x_in);
            }
            break;
        }
        case Side::L: {
            Eigen::Map<VectorTypeT<Scalar>>       x_in    (x_in_ptr,L);
            Eigen::Map<VectorTypeT<Scalar>>       x_out   (x_out_ptr,L);
            if constexpr(std::is_same <Scalar,double>::value){
                x_out.noalias() = x_in *lu_real.inverse();
            }else{
                x_out.noalias() = x_in *lu_cplx.inverse();
            }
            break;
        }
    }
    counter++;
}




template<typename Scalar>
void StlMatrixProduct<Scalar>::MultAx(Scalar* x_in, Scalar* x_out) {
    using namespace eigutils::eigSetting;
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
            x_vec_out.noalias() = A_matrix.template selfadjointView<Eigen::Upper>() * x_vec_in;
            break;
        }
    }
    counter++;
}



// Explicit instantiations

template class StlMatrixProduct<double>;
template class StlMatrixProduct<std::complex<double>>;
