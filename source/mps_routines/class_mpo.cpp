//
// Created by david on 2019-01-17.
//

#include "class_mpo.h"


#include <general/nmspc_quantum_mechanics.h>

Eigen::Tensor<class_mpo::Scalar,4> class_mpo::parity(const Eigen::Matrix2cd & paulimatrix)
/*! Builds the MPO for measuring parity on spin 1/2 systems.
 *      P = Π  s_{i}
 * where Π is the product over all sites, and s_{i} is the given pauli matrix for site i.
 *
 *       | I   0   0 |
 * P_i = | s   0   0 |
 *       | 0   s   I |
 *
 *        2
 *        |
 *    0---P---1
 *        |
 *        3
 *
 */
{
    int spin_dim = 2;
    Eigen::array<long, 4> extent4 = {1, 1, spin_dim, spin_dim};                    /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2> extent2 = {spin_dim, spin_dim};                          /*!< Extent of pauli matrices in a rank-2 tensor */

    Eigen::Tensor<Scalar,4> MPO;
    MPO.resize(3, 3, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(qm::spinOneHalf::I);
    MPO.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(paulimatrix);
    MPO.slice(Eigen::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(paulimatrix);
    MPO.slice(Eigen::array<long, 4>{2, 2, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(qm::spinOneHalf::I);
    return MPO;
}


Eigen::Tensor<class_mpo::Scalar,4> class_mpo::parity(const Eigen::Matrix3cd & paulimatrix)
/*! Builds the MPO for measuring parity on spin 1 systems.
 *      P = Π  s_{i}
 * where Π is the product over all sites, and s_{i} is the given pauli matrix for site i.
 *
 *       | I   0   0 |
 * P_i = | s   0   0 |
 *       | 0   s   I |
 *
 *        2
 *        |
 *    0---P---1
 *        |
 *        3
 *
 */
{
    int spin_dim = 3;
    Eigen::array<long, 4> extent4 = {1, 1, spin_dim, spin_dim};                    /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2> extent2 = {spin_dim, spin_dim};                          /*!< Extent of pauli matrices in a rank-2 tensor */

    Eigen::Tensor<Scalar,4> MPO;
    MPO.resize(3, 3, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(qm::spinOne::I);
    MPO.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(paulimatrix);
    MPO.slice(Eigen::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(paulimatrix);
    MPO.slice(Eigen::array<long, 4>{2, 2, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(qm::spinOne::I);
    return MPO;
}