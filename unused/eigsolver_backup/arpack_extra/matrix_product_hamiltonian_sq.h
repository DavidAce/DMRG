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
#include "math/nmspc_eigutils.h"
#include <array>
#include <general/class_tic_toc.h>
#include <vector>

#define profile_matrix_product_hamiltonian_sq 0

class OMP;

template<class Scalar_>
class DenseHamiltonianSqProduct {
    public:
    using Scalar                    = Scalar_;
    constexpr static bool can_shift = false;

    private:
    const Scalar_ *            env2L_ptr;
    const Scalar_ *            env2R_ptr;
    const Scalar_ *            mpo_ptr;
    std::array<long, 3>        shape_theta;
    std::array<long, 4>        shape_mpo;
    long theta_size;
    eigutils::eigSetting::Form form = eigutils::eigSetting::Form::SYMMETRIC;
    eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R;
    std::shared_ptr<OMP> omp;
    int num_threads = 1; /*!< Number of threads */
    public:
    DenseHamiltonianSqProduct(const Scalar_ *           env2L_,          /*!< The left block tensor.  */
                              const Scalar_ *           env2R_,          /*!< The right block tensor.  */
                              const Scalar_ *           mpo_,            /*!< The left Hamiltonian MPO's  */
                              const std::array<long, 3> shape_theta_,    /*!< An array containing the dimensions of the multisite theta  */
                              const std::array<long, 4> shape_mpo_,      /*!< An array containing the dimensions of the multisite mpo  */
                              const int                 num_threads_ = 1 /*!< Number of threads */
    ) :                           env2L_ptr(env2L_),
                                  env2R_ptr(env2R_),
                                  mpo_ptr(mpo_),
                                  shape_theta(shape_theta_),
                                  shape_mpo(shape_mpo_),
                                  num_threads(num_threads_){
        t_mul.set_properties(profile_matrix_product_hamiltonian_sq, 10, "Time multiplying");
        if(env2L_ptr == nullptr) throw std::runtime_error("env2L is a nullptr!");
        if(env2R_ptr == nullptr) throw std::runtime_error("env2R is a nullptr!");
        if(mpo_ptr   == nullptr) throw std::runtime_error("mpo is a nullptr!");
        theta_size = shape_theta[0]*shape_theta[1]*shape_theta[2];
    }

    // Functions used in in Arpack++ solver
    int rows() const { return (int) theta_size; }; /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */
    int cols() const { return (int) theta_size; }; /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */

    void MultAx(Scalar_ *theta_in_, Scalar_ *theta_out_); //   Computes the matrix-vector multiplication x_out <- A*x_in.

    // Various utility functions
    int                               counter = 0;
    void                              print() const {};
    void                              set_mode(const eigutils::eigSetting::Form form_) { form = form_; }
    void                              set_side(const eigutils::eigSetting::Side side_) { side = side_; }
    const eigutils::eigSetting::Form &get_form() const { return form; }
    const eigutils::eigSetting::Side &get_side() const { return side; }

    // Profiling
    class_tic_toc t_mul;
};

