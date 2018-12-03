//
// Created by david on 2018-11-16.
//

#ifndef MATRIX_PRODUCT_STL_H
#define MATRIX_PRODUCT_STL_H




#include <general/nmspc_eigutils.h>
#include <general/class_tic_toc.h>
#include <iostream>
#include <iomanip>
#include <vector>
#define profile_matrix_product_dense 1


template <typename Scalar_>
class StlMatrixProduct {
public:
    using Scalar      = Scalar_;
private:

    std::vector<Scalar> A_stl;           // The actual matrix. Given matrices will be copied into this one.
    const Scalar *A_ptr;                 // A pointer to the matrix, to allow optional copying of the matrix. Note that PartialPivLU stores LU in A.
    const int L;                         // The linear matrix dimension
    eigutils::eigSetting::Form form;     // Chooses SYMMETRIC / NONSYMMETRIC mode
    eigutils::eigSetting::Side side;     // Chooses whether to find (R)ight or (L)eft eigenvectors

    // Shift-invert mode stuff
    double sigmaR = 0.0;   // The real part of the shift
    double sigmaI = 0.0;   // The imag part of the shift
    bool readyFactorOp = false;                                 // Flag to make sure LU factorization has occurred
    bool readyShift = false;                                    // Flag to make sure

public:
    // Pointer to data constructor, copies the matrix into an internal Eigen matrix.
    StlMatrixProduct(
            const Scalar * const A_,
            const int L_,
            const eigutils::eigSetting::Form form_ = eigutils::eigSetting::Form::SYMMETRIC,
            const eigutils::eigSetting::Side side_ = eigutils::eigSetting::Side::R,
            const bool copy_data = false

    ): A_ptr(A_) ,L(L_), form(form_), side(side_)
    {
        if (copy_data){
            A_stl.resize(L*L);
            std::copy(A_ptr,A_ptr + L*L, A_stl.begin());
            A_ptr = A_stl.data();
        }
        init_profiling();
    }

    ~StlMatrixProduct();


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
    const std::vector<Scalar> get_matrix()const{return std::vector<Scalar>(A_ptr,A_ptr+L*L);}
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







#endif



