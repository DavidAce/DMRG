//
// Created by david on 2018-05-08.
//

#pragma once

//
#ifdef EIGEN_USE_BLAS
#define EIGEN_USE_BLAS_SUSPEND
#undef EIGEN_USE_BLAS
#endif
#ifdef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL_SUSPEND
#undef EIGEN_USE_MKL_ALL
#endif
#ifdef EIGEN_USE_LAPACKE_STRICT
#define EIGEN_USE_LAPACKE_STRICT_SUSPEND
#undef EIGEN_USE_LAPACKE_STRICT
#endif



#include <math/nmspc_eigutils.h>
#include <general/class_tic_toc.h>
#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include <complex.h>
#undef I
#include <Eigen/LU>
#define profile_matrix_product_dense 1

#if defined(_MKL_LAPACK_H_)
#pragma message("_MKL_LAPACK_H_ IS NOT SUPPOSED TO BE DEFINED HERE")
#endif

#if defined(LAPACK_H)
#pragma message("LAPACK IS NOT SUPPOSED TO BE DEFINED HERE")
#endif



template <typename Scalar_>
class DenseMatrixProduct {
public:
    using Scalar      = Scalar_;
    constexpr static bool  can_shift = true;

private:
    using MatrixType  = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::ColMajor>;
    using VectorType  = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    using VectorTypeT = Eigen::Matrix<Scalar,1,Eigen::Dynamic>;

    const MatrixType A_matrix;           // The actual matrix. Given matrices will be copied into this one.
    const int L;                         // The linear matrix dimension
    eigutils::eigSetting::Form form;     // Chooses SYMMETRIC / NONSYMMETRIC mode
    eigutils::eigSetting::Side side;     // Chooses whether to find (R)ight or (L)eft eigenvectors

    // Shift-invert mode stuff
    Eigen::PartialPivLU<MatrixType> lu;                         // Object for dense LU decomposition used in shift-invert mode
    double sigmaR = 0.0;   // The real part of the shift
    double sigmaI = 0.0;   // The imag part of the shift
    bool readyFactorOp = false;                                 // Flag to make sure LU factorization has occurred
    bool readyShift = false;                                    // Flag to make sure

public:
    // Pointer to data constructor, copies the matrix into an internal Eigen matrix.
    DenseMatrixProduct(
            const Scalar * const A_,
            const int L_,
            const eigutils::eigSetting::Form form_ = eigutils::eigSetting::Form::SYMMETRIC,
            const eigutils::eigSetting::Side side_ = eigutils::eigSetting::Side::R

    ): A_matrix(Eigen::Map<const MatrixType>(A_,L_,L_)),
       L(L_), form(form_), side(side_)
    {
        init_profiling();
    }

    // Eigen type constructor. Pass any copy-assignable eigen type into an internal Eigen matrix.
    template<typename Derived>
    explicit DenseMatrixProduct(
            const Eigen::EigenBase<Derived> &matrix_,
            const eigutils::eigSetting::Form form_ = eigutils::eigSetting::Form::SYMMETRIC,
            const eigutils::eigSetting::Side side_ = eigutils::eigSetting::Side::R)
            : A_matrix(matrix_), L(A_matrix.rows()), form(form_), side(side_)
    {
        init_profiling();
    }

    // Functions used in in Arpack++ solver
    int rows() const {return L;};
    int cols() const {return L;};
    void FactorOP();                                      //  Factors (A-sigma*I) into PLU
    void MultOPv(Scalar* x_in_ptr, Scalar* x_out_ptr);    //   Computes the matrix-vector product x_out <- inv(A-sigma*I)*x_in.
    void MultAx (Scalar* x_in_ptr, Scalar* x_out_ptr);    //   Computes the matrix-vector multiplication x_out <- A*x_in.

    // Various utility functions
    int counter = 0;
    void print()const;
    void set_shift(std::complex<double> sigma_)   {if(readyShift){return;} sigmaR=std::real(sigma_);sigmaI=std::imag(sigma_) ;readyShift = true;}
    void set_shift(double               sigma_)   {if(readyShift){return;} sigmaR=sigma_, sigmaI = 0.0;readyShift = true;}
    void set_shift(double sigmaR_, double sigmaI_){if(readyShift){return;} sigmaR=sigmaR_;sigmaI = sigmaI_ ;readyShift = true;}
    void set_mode(const eigutils::eigSetting::Form form_){form = form_;}
    void set_side(const eigutils::eigSetting::Side side_){side = side_;}
    const MatrixType & get_matrix()const{return A_matrix;}
    const eigutils::eigSetting::Form &get_form()const{return form;}
    const eigutils::eigSetting::Side &get_side()const{return side;}

    // Profiling
    void init_profiling(){
        t_factorOp.set_properties(profile_matrix_product_dense, 5,"Time FactorOp");
        t_multOpv.set_properties(profile_matrix_product_dense, 5,"Time MultOpv");
        t_multax.set_properties(profile_matrix_product_dense, 5,"Time MultAx");
    }
    class_tic_toc t_factorOp;
    class_tic_toc t_multOpv;
    class_tic_toc t_multax;
};




// Function definitions



template<typename Scalar>
void DenseMatrixProduct<Scalar>::print() const {
    std::cout << "A_matrix: \n" << A_matrix << std::endl;
}


template<typename Scalar>
void DenseMatrixProduct<Scalar>::FactorOP()

/*  Partial pivot LU decomposition
 *  Factors P(A-sigma*I) = LU
 */
{
    if(readyFactorOp){return;}
    std::cout << "Starting LU \n";

    t_factorOp.tic();
    assert(readyShift and "Shift value sigma has not been set.");
   
    if constexpr(std::is_same<Scalar,double>::value)
    {
        lu.compute(A_matrix - sigmaR * Eigen::MatrixXd::Identity(L,L));
    }
    else
    {
        Scalar sigma = std::complex<double>(sigmaR,sigmaI);
        lu.compute(A_matrix - sigma * Eigen::MatrixXd::Identity(L,L));
    }

    readyFactorOp = true;
    t_factorOp.toc();
    std::cout << "Finished LU \n";
    std::cout << "Time LU Op [ms]: " << std::fixed << std::setprecision(3) << t_factorOp.get_last_time_interval() * 1000 <<'\n';

}




template<typename Scalar>
void DenseMatrixProduct<Scalar>::MultOPv(Scalar* x_in_ptr, Scalar* x_out_ptr) {
    using namespace eigutils::eigSetting;
    assert(readyFactorOp and "FactorOp() has not been run yet.");
    switch (side){
        case Side::R: {
            Eigen::Map<VectorType>       x_in    (x_in_ptr,L);
            Eigen::Map<VectorType>       x_out   (x_out_ptr,L);
            x_out.noalias() = lu.solve(x_in);
            break;
        }
        case Side::L: {
            Eigen::Map<VectorTypeT>       x_in    (x_in_ptr,L);
            Eigen::Map<VectorTypeT>       x_out   (x_out_ptr,L);
            x_out.noalias() = x_in * lu.inverse();
            break;
        }
    }
    counter++;
}




template<typename Scalar>
void DenseMatrixProduct<Scalar>::MultAx(Scalar* x_in, Scalar* x_out) {
    using namespace eigutils::eigSetting;
    switch (form){
        case Form::NONSYMMETRIC:
            switch (side) {
                case Side::R: {
                    Eigen::Map<VectorType> x_vec_in (x_in,  L);
                    Eigen::Map<VectorType> x_vec_out(x_out, L);
                    x_vec_out.noalias() = A_matrix * x_vec_in ;
                    break;
                }
                case Side::L: {
                    Eigen::Map<VectorTypeT> x_vec_in(x_in, L);
                    Eigen::Map<VectorTypeT> x_vec_out(x_out, L);
                    x_vec_out.noalias() = x_vec_in * A_matrix;
                    break;
                }
            }
            break;
        case Form::SYMMETRIC: {
            Eigen::Map<VectorType> x_vec_in(x_in, L);
            Eigen::Map<VectorType> x_vec_out(x_out, L);
            x_vec_out.noalias() = A_matrix.template selfadjointView<Eigen::Upper>() * x_vec_in;
            break;
        }
    }
    counter++;
}



#ifdef EIGEN_USE_BLAS_SUSPEND
#define EIGEN_USE_BLAS
#undef EIGEN_USE_BLAS_SUSPEND
#endif
#ifdef EIGEN_USE_MKL_ALL_SUSPEND
#define EIGEN_USE_MKL_ALL
#undef EIGEN_USE_MKL_ALL_SUSPEND
#endif

#ifdef EIGEN_USE_LAPACKE_STRICT_SUSPEND
#define EIGEN_USE_LAPACKE_STRICT
#undef EIGEN_USE_LAPACKE_STRICT_SUSPEND
#endif

