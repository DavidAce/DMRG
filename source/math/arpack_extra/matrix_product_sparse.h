//
// Created by david on 2018-05-08.
//

#pragma once
#include "matrix_product_fwd_decl.h"
#include <general/class_tic_toc.h>
#include <iomanip>
#include <iostream>
#include <math/nmspc_eigutils.h>

#define profile_matrix_product_sparse 1

template<typename Scalar_, bool sparseLU = false>
class SparseMatrixProduct {
    public:
    using Scalar                    = Scalar_;
    constexpr static bool can_shift = true;

    private:
    std::vector<Scalar>        A_stl; // The actual matrix. Given matrices will be copied into this one if desired.
    const Scalar *             A_ptr; // A pointer to the matrix, to allow optional copying of the matrix.
    const int                  L;     // The linear matrix dimension
    eigutils::eigSetting::Form form;  // Chooses SYMMETRIC / NONSYMMETRIC mode
    eigutils::eigSetting::Side side;  // Chooses whether to find (R)ight or (L)eft eigenvectors

    double sigmaR        = std::numeric_limits<double>::quiet_NaN(); // The real part of the shift
    double sigmaI        = std::numeric_limits<double>::quiet_NaN(); // The imag part of the shift
    bool   readyFactorOp = false;                                    // Flag to make sure LU factorization has occurred
    bool   readyShift    = false;

    public:
    // Pointer to data constructor, copies the matrix into an internal Eigen matrix.
    SparseMatrixProduct(const Scalar *                   A_,
                        const int                        L_,
                        const bool                       copy_data,
                        const eigutils::eigSetting::Form form_     = eigutils::eigSetting::Form::SYMMETRIC,
                        const eigutils::eigSetting::Side side_     = eigutils::eigSetting::Side::R);
    ~SparseMatrixProduct();
    // Functions used in in Arpack++ solver
    int  rows() const { return L; };
    int  cols() const { return L; };
    void FactorOP();                                   //  Factors (A-sigma*I) into PLU
    void MultOPv(Scalar *x_in_ptr, Scalar *x_out_ptr); //   Computes the matrix-vector product x_out <- inv(A-sigma*I)*x_in.
    void MultAx(Scalar *x_in_ptr, Scalar *x_out_ptr);  //   Computes the matrix-vector multiplication x_out <- A*x_in.

    // Various utility functions
    int  counter = 0;
    void print() const;
    void set_shift(std::complex<double> sigma_) {
        if(readyShift) { return; }
        sigmaR     = std::real(sigma_);
        sigmaI     = std::imag(sigma_);
        readyShift = true;
    }
    void set_shift(double sigma_) {
        if(readyShift) { return; }
        sigmaR = sigma_, sigmaI = 0.0;
        readyShift = true;
    }
    void set_shift(double sigmaR_, double sigmaI_) {
        if(readyShift) { return; }
        sigmaR     = sigmaR_;
        sigmaI     = sigmaI_;
        readyShift = true;
    }
    void                              set_mode(const eigutils::eigSetting::Form form_) { form = form_; }
    void                              set_side(const eigutils::eigSetting::Side side_) { side = side_; }
    const eigutils::eigSetting::Form &get_form() const { return form; }
    const eigutils::eigSetting::Side &get_side() const { return side; }

    // Profiling
    void init_profiling() {
        t_factorOp.set_properties(profile_matrix_product_sparse, 5, "Time FactorOp");
        t_multOpv.set_properties(profile_matrix_product_sparse, 5, "Time MultOpv");
        t_multax.set_properties(profile_matrix_product_sparse, 5, "Time MultAx");
    }
    class_tic_toc t_factorOp;
    class_tic_toc t_multOpv;
    class_tic_toc t_multax;
};
