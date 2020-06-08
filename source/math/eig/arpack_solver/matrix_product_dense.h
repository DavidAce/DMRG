//
// Created by david on 2018-11-16.
//

#pragma once

#include "math/eig/enums.h"
#include <iomanip>
#include <iostream>
#include <vector>
class class_tic_toc;

template<typename Scalar_>
class MatrixProductDense {
    public:
    using Scalar                                   = Scalar_;
    constexpr static bool         can_shift_invert = true;
    constexpr static bool         can_shift        = true;
    constexpr static eig::Storage storage          = eig::Storage::DENSE;

    private:
    std::vector<Scalar> A_stl; // The actual matrix. Given matrices will be copied into this one if desired
    const Scalar *      A_ptr; // A pointer to the matrix, to allow optional copying of the matrix. Note that PartialPivLU stores LU in A.
    const long          L;     // The linear matrix dimension
    eig::Form           form;  // Chooses SYMMETRIC / NONSYMMETRIC mode
    eig::Side           side;  // Chooses whether to find (R)ight or (L)eft eigenvectors
    // Shift and shift-invert mode stuff
    std::complex<double> sigma         = 0.0;   // A possibly complex-valued shift
    bool                 readyFactorOp = false; // Flag to make sure LU factorization has occurred
    bool                 readyShift    = false; // Flag to make sure the shift has occurred

    public:
    // Pointer to data constructor, copies the matrix into an internal Eigen matrix.
    MatrixProductDense(const Scalar *const A_, const long L_, const bool copy_data, const eig::Form form_ = eig::Form::SYMM,
                       const eig::Side side_ = eig::Side::R);

    ~MatrixProductDense();

    // Functions used in in Arpack++ solver
    [[nodiscard]] int rows() const { return static_cast<int>(L); };
    [[nodiscard]] int cols() const { return static_cast<int>(L); };
    void              FactorOP();                                   //  Factors (A-sigma*I) into PLU
    void              MultOPv(Scalar *x_in_ptr, Scalar *x_out_ptr); //   Computes the matrix-vector product x_out <- inv(A-sigma*I)*x_in.
    void              MultAx(Scalar *x_in_ptr, Scalar *x_out_ptr);  //   Computes the matrix-vector multiplication x_out <- A*x_in.

    // Various utility functions
    int              counter = 0;
    void             print() const;
    void             set_shift(std::complex<double> sigma_);
    void             set_mode(const eig::Form form_);
    void             set_side(const eig::Side side_);
    const eig::Form &get_form() const;
    const eig::Side &get_side() const;
    // Profiling
    void                           init_profiling();
    std::unique_ptr<class_tic_toc> t_factorOP;
    std::unique_ptr<class_tic_toc> t_multOPv;
    std::unique_ptr<class_tic_toc> t_multAx;
};
