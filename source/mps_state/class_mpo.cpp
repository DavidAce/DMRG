//
// Created by david on 2019-01-17.
//

#include "class_mpo.h"



std::tuple<
        Eigen::Tensor<class_mpo::Scalar,4>,
        Eigen::Tensor<class_mpo::Scalar,3>,
        Eigen::Tensor<class_mpo::Scalar,3>>
class_mpo::pauli_mpo(const Eigen::MatrixXcd paulimatrix)
/*! Builds the MPO for measuring parity on spin systems.
 *      P = Π  s_{i}
 * where Π is the product over all sites, and s_{i} is the given pauli matrix for site i.
 *
 * MPO = | s |
 *
 *        2
 *        |
 *    0---s---1
 *        |
 *        3
 *
 */
{
    long spin_dim = paulimatrix.rows();
    Eigen::array<long, 4> extent4 = {1, 1, spin_dim, spin_dim};                    /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2> extent2 = {spin_dim, spin_dim};                          /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar,4> MPO(1, 1, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(paulimatrix);

    //Create compatible edges
    Eigen::Tensor<Scalar,3> Ledge(1,1,1); // The left  edge
    Eigen::Tensor<Scalar,3> Redge(1,1,1); // The right edge
    Ledge(0,0,0) = 1;
    Redge(0,0,0) = 1;
    return std::make_tuple(MPO,Ledge,Redge);
}




std::tuple<
        Eigen::Tensor<class_mpo::Scalar,4>,
        Eigen::Tensor<class_mpo::Scalar,3>,
        Eigen::Tensor<class_mpo::Scalar,3>>
class_mpo::parity_selector_mpo(const Eigen::MatrixXcd paulimatrix, const int sector)
/*! Builds the MPO that projects out the MPS component in a parity sector.
 *      |psi+->  = O |psi>=  (1 +- P) |psi>
 *  Note that |psi+-> aren't normalized after applying this MPO!
 *
 *       | I   0    |
 * O   = | 0   +-s  |
 *
 *
 *        2
 *        |
 *    0---O---1
 *        |
 *        3
 *
 */
{
    long spin_dim = paulimatrix.rows();
    auto I = Eigen::MatrixXcd::Identity(spin_dim,spin_dim);
    Eigen::array<long, 4> extent4 = {1, 1, spin_dim, spin_dim};                    /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2> extent2 = {spin_dim, spin_dim};                          /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar,4> MPO(2, 2, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(I);
    MPO.slice(Eigen::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(sector * paulimatrix);

    //Create compatible edges
    Eigen::Tensor<Scalar,3> Ledge(1,1,2); // The left  edge
    Eigen::Tensor<Scalar,3> Redge(1,1,2); // The right edge
    Ledge(0,0,0) = 1;
    Ledge(0,0,1) = 1;
    Redge(0,0,0) = 1;
    Redge(0,0,1) = 1;
    return std::make_tuple(MPO,Ledge,Redge);
}



std::tuple<
        std::list<Eigen::Tensor<class_mpo::Scalar,4>>,
        Eigen::Tensor<class_mpo::Scalar,3>,
        Eigen::Tensor<class_mpo::Scalar,3>>
class_mpo::parity_projector_mpos(const Eigen::MatrixXcd paulimatrix, const size_t sites, const int sector)/*! Builds the MPO that projects out the MPS component in a parity sector.
 *      |psi+->  = O |psi>=  (1/2) (1 +- P) |psi>
 *      Here 1 = outer product of 2x2 identity matrices, "sites" times.
 *      Also P = outer product of 2x2 pauli matrices, "sites" times.
 *      The sign in "sector" and the factor 1/2 is put into the left edge at the end.
 *
 *            | I   0  |
 * O   =      | 0   s  |
 *
 *
 *        2
 *        |
 *    0---O---1
 *        |
 *        3
 *
 */
{
    long spin_dim = paulimatrix.rows();
    auto I = Eigen::MatrixXcd::Identity(spin_dim,spin_dim);
    Eigen::array<long, 4> extent4 = {1, 1, spin_dim, spin_dim};                    /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2> extent2 = {spin_dim, spin_dim};                          /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar,4> MPO(2, 2, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(I);
    MPO.slice(Eigen::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(paulimatrix);

    std::list<Eigen::Tensor<class_mpo::Scalar,4>> mpos(sites,MPO);
    //Create compatible edges
    Eigen::Tensor<Scalar,3> Ledge(1,1,2); // The left  edge
    Eigen::Tensor<Scalar,3> Redge(1,1,2); // The right edge
    Ledge(0,0,0) = 1.0/2.0;
    Ledge(0,0,1) = 1.0/2.0 * sector;
    Redge(0,0,0) = 1;
    Redge(0,0,1) = 1;

    return std::make_tuple(mpos,Ledge,Redge);
}