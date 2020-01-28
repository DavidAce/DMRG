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

#define profile_matrix_product_hamiltonian 0

template<class Scalar_>
class DenseHamiltonianProduct {
    public:
    using Scalar                    = Scalar_;
    constexpr static bool can_shift = false;

    private:
    const Scalar_ *            Lblock;
    const Scalar_ *            Rblock;
    const Scalar_ *            HA;
    const Scalar_ *            HB;
    std::array<long, 4>        shape_theta4;
    std::array<long, 2>        shape_theta2;
    std::array<long, 1>        shape_theta1;
    std::array<long, 4>        shape_mpo4;
    eigutils::eigSetting::Form form = eigutils::eigSetting::Form::SYMMETRIC;
    eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R;
    //    std::shared_ptr<OMP> omp;
    public:
    DenseHamiltonianProduct(const Scalar_ *           Lblock_,        /*!< The left block tensor.  */
                            const Scalar_ *           Rblock_,        /*!< The right block tensor.  */
                            const Scalar_ *           HA_,            /*!< The left Hamiltonian MPO's  */
                            const Scalar_ *           HB_,            /*!< The right Hamiltonian MPO's */
                            const std::array<long, 4> shape_theta4_,  /*!< An array containing the shapes of theta  */
                            const std::array<long, 4> shape_mpo4_,    /*!< An array containing the shapes of the MPO  */
                            const size_t              num_threads = 1 /*!< Number of threads */
    );

    // Functions used in in Arpack++ solver
    int rows() const { return (int) shape_theta1[0]; }; /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */
    int cols() const { return (int) shape_theta1[0]; }; /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */

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

template<typename T>
DenseHamiltonianProduct<T>::DenseHamiltonianProduct(const Scalar *            Lblock_,       /*!< The left block tensor.  */
                                                    const Scalar *            Rblock_,       /*!< The right block tensor.  */
                                                    const Scalar *            HA_,           /*!< The left Hamiltonian MPO's  */
                                                    const Scalar *            HB_,           /*!< The right Hamiltonian MPO's */
                                                    const std::array<long, 4> shape_theta4_, /*!< An array containing the shapes of theta  */
                                                    const std::array<long, 4> shape_mpo4_,   /*!< An array containing the shapes of the MPO  */
                                                    const size_t              num_threads)
    : Lblock(Lblock_), Rblock(Rblock_), HA(HA_), HB(HB_), shape_theta4(shape_theta4_),
      shape_theta2({shape_theta4[0] * shape_theta4[1], shape_theta4[2] * shape_theta4[3]}),
      shape_theta1({shape_theta4[0] * shape_theta4[1] * shape_theta4[2] * shape_theta4[3]}), shape_mpo4(shape_mpo4_)
//    omp(std::make_shared<OMP>(num_threads));

{
    t_mul.set_properties(profile_matrix_product_hamiltonian, 10, "Time multiplying");
    if(Lblock == nullptr) throw std::runtime_error("Lblock is a nullptr!");
    if(Rblock == nullptr) throw std::runtime_error("Rblock is a nullptr!");
    if(HA == nullptr) throw std::runtime_error("HA is a nullptr!");
    if(HB == nullptr) throw std::runtime_error("HB is a nullptr!");
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
