//
// Created by david on 2018-05-08.
//

#ifndef DMRG_CLASS_ARPACKPP_CUSTOM_MULTIPLICATION_H
#define DMRG_CLASS_ARPACKPP_CUSTOM_MULTIPLICATION_H

#include "class_arpack_eigsolver.h"

template <typename T, Form form = Form::GENERAL>
class DenseMatrixProduct {
private:
    T *A; //A pointer to the matrix "A"
    int n;    // The dimension of the problem.
    Side side;
public:
    DenseMatrixProduct(int n_, T *A_, Side side_):A(A_), n(n_), side(side_){}

    int counter = 0;
    int rows() const {return n;};
    int cols() const {return n;};
    void MultMv(T* x_in, T* x_out);
};



template<class T>
class DenseHamiltonianProduct {

private:
    const T *Lblock;
    const T *Rblock;
    const T *HA;
    const T *HB;
    std::array<long,4> shape_theta4;
    std::array<long,2> shape_theta2;
    std::array<long,1> shape_theta1;
    std::array<long,4> shape_mpo4;
public:
    int rows()const {return (int)shape_theta1[0];};               /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */
    int cols()const {return (int)shape_theta1[0];};               /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */

    void MultMv(T* theta_in_, T* theta_out_);               /*!< The function that contracts.  */
    int counter = 0;
    DenseHamiltonianProduct(
            const T *Lblock_,                            /*!< The left block tensor.  */
            const T *Rblock_,                            /*!< The right block tensor.  */
            const T *HA_,                                /*!< The left Hamiltonian MPO's  */
            const T *HB_,                                /*!< The right Hamiltonian MPO's */
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
    }
};





#endif //DMRG_CLASS_ARPACKPP_CUSTOM_MULTIPLICATION_H
