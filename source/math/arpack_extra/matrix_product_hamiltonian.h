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

class OMP;

template<class Scalar_>
class DenseHamiltonianProduct {
    public:
    using Scalar                    = Scalar_;
    constexpr static bool can_shift = false;

    private:
    const Scalar_ *            Lblock;
    const Scalar_ *            Rblock;
    const Scalar_ *            mpo;
    std::array<long, 3>        shape_mps;
    std::array<long, 4>        shape_mpo;
    int                        mps_size;
    eigutils::eigSetting::Form form = eigutils::eigSetting::Form::SYMMETRIC;
    eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R;
    std::shared_ptr<OMP>       omp;
    std::optional<int>         num_threads = std::nullopt; /*!< Number of threads */
    public:
    DenseHamiltonianProduct(const Scalar_ *           Lblock_,                    /*!< The left block tensor.  */
                            const Scalar_ *           Rblock_,                    /*!< The right block tensor.  */
                            const Scalar_ *           mpo_,                       /*!< The Hamiltonian MPO's  */
                            const std::array<long, 3> shape_mps_,                 /*!< An array containing the shapes of theta  */
                            const std::array<long, 4> shape_mpo_,                 /*!< An array containing the shapes of the MPO  */
                            std::optional<int>        num_threads_ = std::nullopt /*!< Number of threads */
                            )
        : Lblock(Lblock_), Rblock(Rblock_), mpo(mpo_), shape_mps(shape_mps_), shape_mpo(shape_mpo_), num_threads(num_threads_) {
        t_mul.set_properties(profile_matrix_product_hamiltonian, 10, "Time multiplying");
        if(Lblock == nullptr) throw std::runtime_error("Lblock is a nullptr!");
        if(Rblock == nullptr) throw std::runtime_error("Rblock is a nullptr!");
        if(mpo == nullptr) throw std::runtime_error("mpo is a nullptr!");
        mps_size = static_cast<int>(shape_mps[0] * shape_mps[1] * shape_mps[2]);
    }

    // Functions used in in Arpack++ solver
    [[nodiscard]] int rows() const { return mps_size; }; /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */
    [[nodiscard]] int cols() const { return mps_size; }; /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */

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
