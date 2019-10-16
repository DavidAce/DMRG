//
// Created by david on 2018-05-08.
//

#pragma once


#ifdef EIGEN_USE_BLAS
#define EIGEN_USE_BLAS_SUSPEND
#undef EIGEN_USE_BLAS
#endif

#include <math/nmspc_eigutils.h>
#include <general/class_tic_toc.h>
#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
//#include <Eigen/SuperLUSupport>
#define profile_matrix_product_sparse 1


template <typename Scalar_>
class SparseMatrixProduct {
public:
    using Scalar          = Scalar_;
    constexpr static bool  can_shift = true;
private:
    using MatrixType      = Eigen::SparseMatrix<Scalar>;
    using DenseMatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
    using VectorType      = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    using VectorTypeT     = Eigen::Matrix<Scalar,1,Eigen::Dynamic>;

    MatrixType A_matrix;          // The actual matrix. Given matrices will be copied into this one.
    const int L;                        // The linear matrix dimension
    eigutils::eigSetting::Form form;     // Chooses SYMMETRIC / NONSYMMETRIC mode
    eigutils::eigSetting::Side side;     // Chooses whether to find (R)ight or (L)eft eigenvectors

    // Shift-invert mode stuff
//    MatrixType A_shift;
//    Eigen::SuperLU<MatrixType>           lu_dense;
    Eigen::SparseLU<MatrixType>          lu_sparse;              // Object for sparse LU decomposition used in shift-invert mode
    Eigen::PartialPivLU<DenseMatrixType> lu_dense;               // Object for dense LU decomposition used in shift-invert mode
//    Eigen::SparseLU<Eigen::SparseMatrix<Scalar, Eigen::ColMajor, Eigen::Index>, Eigen::COLAMDOrdering<Eigen::Index> > lu_sparse;

    double sigmaR = std::numeric_limits<double>::quiet_NaN();   // The real part of the shift
    double sigmaI = std::numeric_limits<double>::quiet_NaN();   // The imag part of the shift
    bool sparseLU      = false;
    bool readyFactorOp = false;                                 // Flag to make sure LU factorization has occurred
    bool readyShift = false;

public:
    // Pointer to data constructor, copies the matrix into an internal Eigen matrix.
    SparseMatrixProduct(
            const Scalar * A_,
            const int L_,
            const eigutils::eigSetting::Form form_ = eigutils::eigSetting::Form::SYMMETRIC,
            const eigutils::eigSetting::Side side_ = eigutils::eigSetting::Side::R)
            : A_matrix(Eigen::Map<const DenseMatrixType>(A_,L_,L_).sparseView()), L(L_), form(form_), side(side_)
    {
        A_matrix.makeCompressed();
        init_profiling();
    }

    // Eigen type constructor. Pass any copy-assignable eigen type into an internal Eigen matrix.
    template<typename Derived>
    explicit SparseMatrixProduct(
            const Eigen::EigenBase<Derived> &matrix_,
            const eigutils::eigSetting::Form form_ = eigutils::eigSetting::Form::SYMMETRIC,
            const eigutils::eigSetting::Side side_ = eigutils::eigSetting::Side::R)
            : A_matrix(matrix_.derived().sparseView()), L(A_matrix.rows()), form(form_), side(side_)
    {
        A_matrix.makeCompressed();
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
    const MatrixType & get_matrix()const {return A_matrix;}
    const eigutils::eigSetting::Form &get_form()const{return form;}
    const eigutils::eigSetting::Side &get_side()const{return side;}

    // Profiling
    void init_profiling(){
        t_factorOp.set_properties(profile_matrix_product_sparse, 5,"Time FactorOp");
        t_multOpv.set_properties(profile_matrix_product_sparse, 5,"Time MultOpv");
        t_multax.set_properties(profile_matrix_product_sparse, 5,"Time MultAx");
    }
    class_tic_toc t_factorOp;
    class_tic_toc t_multOpv;
    class_tic_toc t_multax;

};




// Function definitions



template<typename Scalar>
void SparseMatrixProduct<Scalar>::print() const {
    std::cout << "A_matrix: \n" << A_matrix << std::endl;
}


template<typename Scalar>
void SparseMatrixProduct<Scalar>::FactorOP()

/*  Sparse decomposition
 *  Factors P(A-sigma*I) = LU
 */

{   if(readyFactorOp){return;}
    assert(readyShift and "Shift value sigma has not been set.");
    Scalar sigma;
    t_factorOp.tic();
    if constexpr(std::is_same<Scalar,double>::value)
    {
        sigma = sigmaR;
        lu_dense.compute(A_matrix - sigmaR * Eigen::MatrixXd::Identity(L,L));
    }
    else
    {
        sigma = std::complex<double>(sigmaR,sigmaI);
        lu_dense.compute(A_matrix - sigma * Eigen::MatrixXd::Identity(L,L));
    }

    sparseLU = false;
    t_factorOp.toc();
    readyFactorOp = true;
    std::cout << "Time Factor Op [ms]: " << std::fixed << std::setprecision(3) << t_factorOp.get_last_time_interval() * 1000 <<'\n';
}

template<typename Scalar>
void SparseMatrixProduct<Scalar>::MultOPv(Scalar* x_in_ptr, Scalar* x_out_ptr) {
    using namespace eigutils::eigSetting;
    assert(readyFactorOp and "FactorOp() has not been run yet.");
    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    Eigen::Map<VectorType>       x_in    (x_in_ptr,L);
    Eigen::Map<VectorType>       x_out   (x_out_ptr,L);
    switch (side){
        case Side::R: {
            if(sparseLU){x_out.noalias() = lu_sparse.solve(x_in);}
            else        {x_out.noalias() = lu_dense.solve(x_in);}
            break;
        }
        case Side::L: {
            std::cerr << "Left sided sparse shift invert hasn't been implemented yet..." << std::endl;
            exit(1);
            break;
        }
    }
    counter++;
}

template<typename Scalar>
void SparseMatrixProduct<Scalar>::MultAx(Scalar* x_in, Scalar* x_out) {
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
            Eigen::Map<VectorType> x_vec_in (x_in,  L);
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


