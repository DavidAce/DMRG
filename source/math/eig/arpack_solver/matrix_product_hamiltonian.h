//
// Created by david on 2018-10-30.
//

#pragma once
//#ifdef EIGEN_USE_BLAS
//#define EIGEN_USE_BLAS_SUSPEND
//#undef EIGEN_USE_BLAS
//#endif
//#ifdef EIGEN_USE_MKL_ALL
//#define EIGEN_USE_MKL_ALL_SUSPEND
//#undef EIGEN_USE_MKL_ALL
//#endif
//#ifdef EIGEN_USE_LAPACKE_STRICT
//#define EIGEN_USE_LAPACKE_STRICT_SUSPEND
//#undef EIGEN_USE_LAPACKE_STRICT
//#endif
#include "math/eig/enums.h"
#include <array>
#include <general/class_tic_toc.h>
#include <vector>

#define profile_matrix_product_hamiltonian 0

class OMP;

template<class Scalar_>
class MatrixProductHamiltonian {
    public:
    using Scalar                                   = Scalar_;
    constexpr static bool         can_shift_invert = false;
    constexpr static bool         can_shift        = true;
    constexpr static eig::Storage storage          = eig::Storage::TENSOR;

    private:
    const Scalar_ *      Lblock;
    const Scalar_ *      Rblock;
    const Scalar_ *      mpo;
    std::array<long, 3>  shape_mps;
    std::array<long, 4>  shape_mpo;
    std::vector<Scalar>  shift_mpo;
    long                 mps_size;
    eig::Form            form = eig::Form::SYMM;
    eig::Side            side = eig::Side::R;
    std::shared_ptr<OMP> omp;
    std::optional<int>   num_threads = std::nullopt; /*!< Number of threads */

    // Shift and shift-invert mode stuff
    std::complex<double> sigma      = 0.0;   // The real part of the shift
    bool                 readyShift = false; // Flag to make sure the shift has occurred

    public:
    MatrixProductHamiltonian(const Scalar_ *     Lblock_,                    /*!< The left block tensor.  */
                             const Scalar_ *     Rblock_,                    /*!< The right block tensor.  */
                             const Scalar_ *     mpo_,                       /*!< The Hamiltonian MPO's  */
                             std::array<long, 3> shape_mps_,                 /*!< An array containing the shapes of theta  */
                             std::array<long, 4> shape_mpo_,                 /*!< An array containing the shapes of the MPO  */
                             std::optional<int>  num_threads_ = std::nullopt /*!< Number of threads */
    );

    // Functions used in in Arpack++ solver
    [[nodiscard]] int rows() const { return static_cast<int>(mps_size); }; /*!< Linear size\f$d^2 \times \chi_L \times \chi_R \f$  */
    [[nodiscard]] int cols() const { return static_cast<int>(mps_size); }; /*!< Linear size\f$d^2 \times \chi_L \times \chi_R \f$  */

    void MultAx(Scalar_ *theta_in_, Scalar_ *theta_out_); //   Computes the matrix-vector multiplication x_out <- A*x_in.

    // Various utility functions
    int                            counter = 0;
    void                           print() const;
    void                           set_shift(std::complex<double> sigma_);
    void                           set_mode(eig::Form form_);
    void                           set_side(eig::Side side_);
    [[nodiscard]] const eig::Form &get_form() const;
    [[nodiscard]] const eig::Side &get_side() const;
    // Profiling
    void                           init_profiling();
    std::unique_ptr<class_tic_toc> t_multAx;
};
