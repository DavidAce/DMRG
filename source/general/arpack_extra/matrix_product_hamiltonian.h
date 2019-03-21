//
// Created by david on 2018-10-30.
//

#ifndef MATRIX_PRODUCT_HAMILTONIAN_H
#define MATRIX_PRODUCT_HAMILTONIAN_H
#include <general/class_tic_toc.h>
#include <array>
#include <vector>
#include "general/nmspc_eigutils.h"
#include <general/nmspc_tensor_extra.h>

#define profile_matrix_product_hamiltonian 0



template<class Scalar_>
class DenseHamiltonianProduct {
public:
    using Scalar          = Scalar_;
    constexpr static bool  can_shift = false;
private:
    const Scalar_ *Lblock;
    const Scalar_ *Rblock;
    const Scalar_ *HA;
    const Scalar_ *HB;
    std::array<long,4> shape_theta4;
    std::array<long,2> shape_theta2;
    std::array<long,1> shape_theta1;
    std::array<long,4> shape_mpo4;
    eigutils::eigSetting::Form form = eigutils::eigSetting::Form::SYMMETRIC;
    eigutils::eigSetting::Side side = eigutils::eigSetting::Side::R;

public:

    DenseHamiltonianProduct(
            const Scalar_ *Lblock_,                            /*!< The left block tensor.  */
            const Scalar_ *Rblock_,                            /*!< The right block tensor.  */
            const Scalar_ *HA_,                                /*!< The left Hamiltonian MPO's  */
            const Scalar_ *HB_,                                /*!< The right Hamiltonian MPO's */
            const std::array<long,4> shape_theta4_,      /*!< An array containing the shapes of theta  */
            const std::array<long,4> shape_mpo4_         /*!< An array containing the shapes of the MPO  */
    ):                                                   /*!< Initializes the custom contraction. */
            Lblock(Lblock_),
            Rblock(Rblock_),
            HA(HA_),
            HB(HB_),
            shape_theta4(shape_theta4_),
            shape_theta2({shape_theta4[0] * shape_theta4[1] , shape_theta4[2] * shape_theta4[3]}),
            shape_theta1({shape_theta4[0] * shape_theta4[1] * shape_theta4[2] * shape_theta4[3]}),
            shape_mpo4(shape_mpo4_)
    {
        t_mul.set_properties(profile_matrix_product_hamiltonian, 10,"Time multiplying");

    }

    // Functions used in in Arpack++ solver
    int rows()const {return (int)shape_theta1[0];};               /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */
    int cols()const {return (int)shape_theta1[0];};               /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */

    void MultAx(Scalar_* theta_in_, Scalar_* theta_out_); //   Computes the matrix-vector multiplication x_out <- A*x_in.

    // Various utility functions
    int counter = 0;
    void print()const{};
    void set_mode(const eigutils::eigSetting::Form form_){form = form_;}
    void set_side(const eigutils::eigSetting::Side side_){side = side_;}
    const eigutils::eigSetting::Form &get_form()const{return form;}
    const eigutils::eigSetting::Side &get_side()const{return side;}

    //Profiling
    class_tic_toc t_mul;
};


template<class T>
void DenseHamiltonianProduct<T>::MultAx(T* theta_in_, T* theta_out_) {
    t_mul.tic();
    Eigen::TensorMap<Eigen::Tensor<const T, 3>>       Lblock_map(Lblock,shape_theta4[1], shape_theta4[1], shape_mpo4[0] );
    Eigen::TensorMap<Eigen::Tensor<const T, 3>>       Rblock_map(Rblock,shape_theta4[3], shape_theta4[3], shape_mpo4[1] );
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       HA_map    (HA, shape_mpo4);
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       HB_map    (HB, shape_mpo4);
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       theta_in  (theta_in_, shape_theta4);
    Eigen::TensorMap<Eigen::Tensor<T, 4>>             theta_out (theta_out_, shape_theta4);



    //Best yet! The sparcity of the effective hamiltonian (Lblock HA HB Rblock) is about 58% nonzeros.
    //L have shown this to be the fastest contraction ordering
    theta_out = Lblock_map
            .contract(theta_in,    Textra::idx({0},{1}))
            .contract(HA_map ,     Textra::idx({1,2},{0,2}))//  idx({1,2,3},{0,4,5}))
            .contract(HB_map ,     Textra::idx({3,1},{0,2}))//  idx({1,2,3},{0,4,5}))
            .contract(Rblock_map,  Textra::idx({1,3},{0,2}))
            .shuffle(Textra::array4{1,0,2,3});
    counter++;
    t_mul.toc();
}










//template<class Scalar_>
//class SparseHamiltonianProduct {
//public:
//    using Scalar          = Scalar_;
//    using MatrixType      = Eigen::SparseMatrix<Scalar>;
//    using DenseMatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
//    using VectorType      = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
//    using VectorTypeT     = Eigen::Matrix<Scalar,1,Eigen::Dynamic>;
//private:
//    const Scalar_ *Lblock;
//    const Scalar_ *Rblock;
//    const Scalar_ *HA;
//    const Scalar_ *HB;
//    std::array<long,4> shape_theta4;
//    std::array<long,2> shape_theta2;
//    std::array<long,1> shape_theta1;
//    std::array<long,4> shape_mpo4;
//public:
//    int rows()const {return (int)shape_theta1[0];};               /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */
//    int cols()const {return (int)shape_theta1[0];};               /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */
//
//    void MultMv(Scalar_* theta_in_, Scalar_* theta_out_);               /*!< The function that contracts.  */
//    int counter = 0;
//    SparseHamiltonianProduct(
//            const Scalar_ *Lblock_,                            /*!< The left block tensor.  */
//            const Scalar_ *Rblock_,                            /*!< The right block tensor.  */
//            const Scalar_ *HA_,                                /*!< The left Hamiltonian MPO's  */
//            const Scalar_ *HB_,                                /*!< The right Hamiltonian MPO's */
//            const std::array<long,4> shape_theta4_,      /*!< An array containing the shapes of theta  */
//            const std::array<long,4> shape_mpo4_         /*!< An array containing the shapes of the MPO  */
//    ):                                                   /*!< Initializes the custom contraction. */
//            Lblock(Lblock_),
//            Rblock(Rblock_),
//            HA(HA_),
//            HB(HB_),
//            shape_theta4(shape_theta4_),
//            shape_theta2({shape_theta4[0] * shape_theta4[1] , shape_theta4[2] * shape_theta4[3]}),
//            shape_theta1({shape_theta4[0] * shape_theta4[1] * shape_theta4[2] * shape_theta4[3]}),
//            shape_mpo4(shape_mpo4_)
//    {
//        t_mul.set_properties(profile_matrix_product_hamiltonian, 10,"Time multiplying");
//
//    }
//
//
//    //Profiling
//    class_tic_toc t_mul;
//};


#endif //MATRIX_PRODUCT_HAMILTONIAN_H
