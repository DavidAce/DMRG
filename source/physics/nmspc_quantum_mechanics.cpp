//
// Created by david on 2019-03-20.
//
//#include <complex.h>
//#undef I

#include <Eigen/Core>
#include <vector>

#include "general/nmspc_tensor_extra.h"
#include "nmspc_quantum_mechanics.h"
#include <math/num.h>
#include <math/rnd.h>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

using Scalar = std::complex<double>;
using namespace Eigen;



std::vector<Eigen::MatrixXcd> qm::gen_manybody_spin(const Eigen::MatrixXcd &s, int sites) {
    std::vector<MatrixXcd> S;
    MatrixXcd              id = MatrixXcd::Identity(s.rows(), s.cols());
    MatrixXcd              tmp;
    for(int idx = 0; idx < sites; idx++) {
        tmp = idx == 0 ? s : id;
        for(int j = 1; j < sites; j++) {
            tmp = kroneckerProduct(tmp, idx == j ? s : id).eval();
        }
        S.emplace_back(tmp);
    }
    return S;
}

namespace qm::spinOneHalf {

    /* clang-format off */
    Matrix2cd sx = (Matrix2cd() <<
            0.0, 1.0,
            1.0, 0.0).finished();
    Matrix2cd sy = (Matrix2cd() <<
            0.0, imn,
            imp, 0.0).finished();
    Matrix2cd sz = (Matrix2cd() <<
            1.0, 0.0,
            0.0, -1.0).finished();
    Matrix2cd sp = (Matrix2cd() <<
            0.0, 2.0,
            0.0, 0.0).finished();
    Matrix2cd sm = (Matrix2cd() <<
            0.0, 0.0,
            2.0, 0.0).finished();
    Matrix2cd id  = (Matrix2cd() << 1.0, 0.0,
            0.0, 1.0).finished();

    std::array<Vector2cd,2> sx_spinors{(Vector2cd() << 1.0, 1.0).finished()/std::sqrt(2),
                                       (Vector2cd() << 1.0,-1.0).finished()/std::sqrt(2)};

    std::array<Vector2cd,2> sy_spinors{(Vector2cd() << 1.0, imp).finished()/std::sqrt(2),
                                       (Vector2cd() << 1.0, imn).finished()/std::sqrt(2)};

    std::array<Vector2cd,2> sz_spinors{(Vector2cd() << 1.0, 0.0).finished()/std::sqrt(2),
                                       (Vector2cd() << 0.0, 1.0).finished()/std::sqrt(2)};

    std::vector<Eigen::MatrixXcd> SX;
    std::vector<Eigen::MatrixXcd> SY;
    std::vector<Eigen::MatrixXcd> SZ;
    std::vector<Eigen::MatrixXcd> II;
    /* clang-format on */

    std::vector<Eigen::Matrix4cd> gen_twobody_spin(const Eigen::Matrix2cd &s) {
        std::vector<Matrix4cd> S;
        MatrixXcd tmp;
        for(int idx = 0; idx < 2; idx++) {
            tmp = idx == 0 ? s : id;
            for(int j = 1; j < 2; j++) {
                tmp = kroneckerProduct(tmp, idx == j ? s : id).eval();
            }
            S.emplace_back(tmp);
        }
        return S;
    }


}

namespace qm::SpinOne {
    /* clang-format off */

    Matrix3cd sx = (Matrix3cd() <<  0.0, 1.0, 0.0,
                                    1.0, 0.0, 1.0,
                                    0.0, 1.0, 0.0).finished();
    Matrix3cd sy = (Matrix3cd() <<  0.0 , imn, 0.0,
                                    imp,  0.0 ,imn,
                                    0.0 ,  imp, 0.0).finished();
    Matrix3cd sz = (Matrix3cd() <<  1.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0,
                                    0.0, 0.0,-1.0).finished();
    Matrix3cd id = (Matrix3cd()  << 1.0, 0.0, 0.0,
                                    0.0, 1.0, 0.0,
                                    0.0, 0.0, 1.0).finished();

    std::vector<Eigen::MatrixXcd> SX;
    std::vector<Eigen::MatrixXcd> SY;
    std::vector<Eigen::MatrixXcd> SZ;
    std::vector<Eigen::MatrixXcd> II;
    /* clang-format on */

    std::vector<Eigen::MatrixXcd> gen_twobody_spin(const Eigen::Matrix3cd &s) {
        std::vector<MatrixXcd> S;
        MatrixXcd tmp;
        for(int idx = 0; idx < 2; idx++) {
            tmp = idx == 0 ? s : id;
            for(int j = 1; j < 2; j++) {
                tmp = kroneckerProduct(tmp, idx == j ? s : id).eval();
            }
            S.emplace_back(tmp);
        }
        return S;
    }

}

namespace qm::timeEvolution {

    std::vector<Eigen::MatrixXcd> Suzuki_Trotter_1st_order(std::complex<double> t, const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd) {
        return {(t * h_evn).exp(), (t * h_odd).exp()};
    }

    std::vector<Eigen::MatrixXcd> Suzuki_Trotter_2nd_order(std::complex<double> t, const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd) {
        return {(t * h_evn / 2.0).exp(), (t * h_odd).exp(), (t * h_evn / 2.0).exp()};
    }

    std::vector<Eigen::MatrixXcd> Suzuki_Trotter_4th_order(std::complex<double> t, const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd)
    /*!
     * Implementation based on
     * Janke, W., & Sauer, T. (1992).
     * Properties of higher-order Trotter formulas.
     * Physics Letters A, 165(3), 199–205.
     * https://doi.org/10.1016/0375-9601(92)90035-K
     *
     */
    {
        double cbrt2 = pow(2.0, 1.0 / 3.0);
        double beta1 = 1.0 / (2.0 - cbrt2);
        double beta2 = -cbrt2 * beta1;
        double alph1 = 0.5 * beta1;
        double alph2 = (1.0 - cbrt2) / 2.0 * beta1;

        std::vector<Eigen::MatrixXcd> temp;

        temp.emplace_back((alph1 * t * h_evn).exp());
        temp.emplace_back((beta1 * t * h_odd).exp());
        temp.emplace_back((alph2 * t * h_evn).exp());
        temp.emplace_back((beta2 * t * h_odd).exp());
        temp.emplace_back((alph2 * t * h_evn).exp());
        temp.emplace_back((beta1 * t * h_odd).exp());
        temp.emplace_back((alph1 * t * h_evn).exp());
        return temp;
    }

    std::vector<Eigen::Tensor<std::complex<double>, 2>> get_2site_evolution_gates(std::complex<double> t, size_t susuki_trotter_order,
                                                                                  const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd)
    /*! Returns a set of 2-site unitary gates, using Suzuki Trotter decomposition to order 1, 2 or 3.
     * These gates need to be applied to the MPS one at a time with a swap in between.
     */
    {
        long                          spin_dim_evn = h_evn.rows();
        long                          spin_dim_odd = h_odd.rows();
        std::vector<Eigen::MatrixXcd> matrix_vec;
        switch(susuki_trotter_order) {
            case 1: matrix_vec = Suzuki_Trotter_1st_order(t, h_evn, h_odd); break;
            case 2: matrix_vec = Suzuki_Trotter_2nd_order(t, h_evn, h_odd); break;
            case 4: matrix_vec = Suzuki_Trotter_4th_order(t, h_evn, h_odd); break;
            default: matrix_vec = Suzuki_Trotter_2nd_order(t, h_evn, h_odd); break;
        }
        std::vector<Eigen::Tensor<std::complex<double>, 2>> tensor_vec;
        for(auto &m : matrix_vec) {
            tensor_vec.emplace_back(Textra::MatrixTensorMap(m, spin_dim_evn, spin_dim_odd)); // spin_dim^2  *
        }
        return tensor_vec;
    }



    std::vector<Eigen::Tensor<std::complex<double>, 2>> compute_G(const std::complex<double> a, size_t susuki_trotter_order, const Eigen::MatrixXcd &h_evn,
                                                                  const Eigen::MatrixXcd &h_odd)
    /*! Returns the moment generating function, or characteristic function (if a is imaginary) for the Hamiltonian as a rank 2 tensor.
     *  The legs contain two physical spin indices each
    *   G := exp(iaM) or exp(aM), where a is a small parameter and M is an MPO.
    *   Note that G(-a) = G(a)* if  exp(iaM) !
    *
    @verbatim
                     0
                     |
                [ exp(aH) ]
                     |
                     1
    @endverbatim
    */
    {
        return get_2site_evolution_gates(a, susuki_trotter_order, h_evn, h_odd);
    }

}

std::vector<Eigen::Tensor<qm::cplx,2>> qm::lbit::get_unitary_twosite_operators(size_t sites, double fmix){
    /*! Returns a set of unitary two site operators used to transform between physical and l-bit representations
     *
     *
     @verbatim
                   0      1                0
                   |      |                |
                 [ exp(-ifH) ]  ---> [ exp(-ifH) ]
                   |      |               |
                   2      3               1
     @endverbatim
    */

    tools::log->trace("Generating twosite unitaries");
    std::vector<Eigen::Tensor<qm::cplx,2>> unitaries(sites-1);
    auto SZ  = qm::spinOneHalf::gen_twobody_spin(qm::spinOneHalf::sz);
    auto SP  = qm::spinOneHalf::gen_twobody_spin(qm::spinOneHalf::sp);
    auto SM  = qm::spinOneHalf::gen_twobody_spin(qm::spinOneHalf::sm);
    auto ID  = qm::spinOneHalf::gen_twobody_spin(qm::spinOneHalf::id);
    auto N   = std::vector<Eigen::Matrix4cd> {0.5 * (ID[0] + SZ[0]), 0.5 * (ID[1] + SZ[1])};

    for (size_t idx = 0; idx < sites-1; idx++){
        double th0 = rnd::uniform_double_box(-1,1);
        double th1 = rnd::uniform_double_box(-1,1);
        double th2 = rnd::uniform_double_box(-1,1);
        double th3 = rnd::uniform_double_box(-1,1);
        std::complex<double> t(rnd::uniform_double_box(-1,1), rnd::uniform_double_box(-1,1));

        Eigen::Matrix4cd H =
            th3*N[0]*N[1] +
            th2*N[1]*(ID[0]-N[0]) +
            th1*N[0]*(ID[1]-N[1]) +
            th0*(ID[0]-N[0])*(ID[1]-N[1]) +
            SP[0]* SM[1]*t +
            SP[1]* SM[0]*std::conj(t);

        unitaries[idx] = Textra::MatrixToTensor((imn * fmix * H).exp(),4,4);
    }
    return unitaries;
}


std::tuple<Eigen::Tensor<Scalar, 4>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>> qm::mpo::pauli_mpo(const Eigen::MatrixXcd &paulimatrix)
/*! Builds the MPO string for measuring  spin on many-body systems.
*      P = Π  s_{i}
* where Π is the product over all sites, and s_{i} is the given pauli matrix for site i.
*
* MPO = | s | (a 1 by 1 matrix with a single pauli matrix element)
*
*        2
*        |
*    0---s---1
*        |
*        3
*
*/
{
    long                     spin_dim = paulimatrix.rows();
    Eigen::array<long, 4>    extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2>    extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO(1, 1, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix);

    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, 1); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, 1); // The right edge
    Ledge(0, 0, 0) = 1;
    Redge(0, 0, 0) = 1;
    return std::make_tuple(MPO, Ledge, Redge);
}

std::tuple<std::list<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
    qm::mpo::parity_projector_mpos(const Eigen::MatrixXcd &paulimatrix, size_t sites, int sign)
/*! Builds the MPO that projects out the MPS component in a parity sector.
 *      |psi+->  = O |psi>=  (1/2) (1 +- P) |psi>
 *      Here 1 = outer product of 2x2 identity matrices, "sites" times, i.e. Kron_(i=0)^(sites-1) I_(2x2)
 *      Also P = outer product of 2x2 pauli matrices, "sites" times, i.e. Kron_(i=0)^(sites-1) s_(2x2)
 *      The sign and the factor 1/2 is put into the left edge at the end.
 *
 *                  | I   0  |
 * O   =   (1/2) *  | 0   s  |
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
    long                     spin_dim = paulimatrix.rows();
    auto                     I        = Eigen::MatrixXcd::Identity(spin_dim, spin_dim);
    Eigen::array<long, 4>    extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2>    extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO(2, 2, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(I);
    MPO.slice(Eigen::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix);

    std::list<Eigen::Tensor<Scalar, 4>> mpos(sites, MPO);
    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, 2); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, 2); // The right edge
    Ledge(0, 0, 0) = 0.5;
    Ledge(0, 0, 1) = 0.5 * sign;
    Redge(0, 0, 0) = 1;
    Redge(0, 0, 1) = 1;

    return std::make_tuple(mpos, Ledge, Redge);
}

std::tuple<std::list<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
    qm::mpo::random_pauli_mpos(const Eigen::MatrixXcd &paulimatrix, size_t sites)
/*! Builds a string of random pauli matrix MPO's
 *      P = Π  O_i
 * where Π is the product over all sites, and O_i is one of {S, I} on site i, where S and I is a pauli matrix or an identity matrix respectively
 *
 * MPO = | s |
 *
 *        2
 *        |
 *    0---O---1
 *        |
 *        3
 *
 */
{
    long                     spin_dim = paulimatrix.rows();
    Eigen::array<long, 4>    extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2>    extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO_I(1, 1, spin_dim, spin_dim);
    Eigen::Tensor<Scalar, 4> MPO_S(1, 1, spin_dim, spin_dim);
    MPO_I.setZero();
    MPO_S.setZero();
    MPO_I.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix);
    MPO_S.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(Eigen::MatrixXcd::Identity(spin_dim, spin_dim));

    // We have to push in an even number of pauli matrices to retain the parity sector.
    // Choosing randomli
    std::vector<int> binary(sites, -1);
    int              sum = 0;
    while(true) {
        binary[rnd::uniform_integer_box<size_t>(0, sites)] *= -1;
        sum = std::accumulate(binary.begin(), binary.end(), 0);
        if((num::mod<size_t>(sites, 2) == 0 and sum == 0) or (num::mod<size_t>(sites, 2) == 1 and sum == 1)) break;
    }

    std::list<Eigen::Tensor<Scalar, 4>> mpos;
    for(auto &val : binary) {
        if(val < 0)
            mpos.push_back(MPO_S);
        else
            mpos.push_back(MPO_I);
    }

    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, 1); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, 1); // The right edge
    Ledge(0, 0, 0) = 1;
    Redge(0, 0, 0) = 1;
    return std::make_tuple(mpos, Ledge, Redge);
}

std::tuple<std::list<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
    qm::mpo::random_pauli_mpos_x2(const Eigen::MatrixXcd &paulimatrix1, const Eigen::MatrixXcd &paulimatrix2, const size_t sites)
/*! Builds a string of random pauli matrix MPO's
 *      P = Π  O_i
 * where Π is the product over all sites, and O_i is one of {S, I} on site i.
 * S is the sum of pauli matrices s1 and s2, and where I is an identity matrix of the same size
 *            | s1  0  |
 * S   =      | 0   s2 |
 *
 *            | id  0  |
 * I   =      | 0   id |
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
    if(paulimatrix1.rows() != paulimatrix2.rows()) throw std::logic_error("Pauli matrices must be of equal size");
    long                     spin_dim = paulimatrix1.rows();
    auto                     I        = Eigen::MatrixXcd::Identity(spin_dim, spin_dim);
    Eigen::array<long, 4>    extent4  = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2>    extent2  = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO_S(2, 2, spin_dim, spin_dim);
    Eigen::Tensor<Scalar, 4> MPO_I(2, 2, spin_dim, spin_dim);
    Eigen::Tensor<Scalar, 4> MPO_P(2, 2, spin_dim, spin_dim);
    MPO_S.setZero();
    MPO_I.setZero();
    MPO_P.setZero();

    MPO_S.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix1);
    MPO_S.slice(Eigen::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix2);
    MPO_I.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(I);
    MPO_I.slice(Eigen::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(I);
    MPO_P.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(I);
    MPO_P.slice(Eigen::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix1);

    // Push in an even number of operators
    std::vector<int> binary(sites, -1);
    int              sum = 0;
    while(true) {
        binary[rnd::uniform_integer_box<size_t>(0, sites)] *= -1;
        sum = std::accumulate(binary.begin(), binary.end(), 0);
        if((num::mod<size_t>(sites, 2) == 0 and sum == 0) or (num::mod<size_t>(sites, 2) == 1 and sum == 1)) break;
    }
    if(binary.size() != sites) throw std::logic_error("Size mismatch");
    // Generate the list
    std::list<Eigen::Tensor<Scalar, 4>> mpos;
    for(auto &val : binary) {
        if(val < 0)
            mpos.push_back(MPO_S);
        else
            mpos.push_back(MPO_I);
    }

    //    std::list<Eigen::Tensor<Scalar,4>> mpos(sites,MPO_P);
    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, 2); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, 2); // The right edge
    Ledge(0, 0, 0) = 1.0 / 2.0;
    Ledge(0, 0, 1) = 1.0 / 2.0;
    Redge(0, 0, 0) = 1;
    Redge(0, 0, 1) = 1;
    return std::make_tuple(mpos, Ledge, Redge);
}

std::tuple<std::list<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
    qm::mpo::random_pauli_mpos(const std::vector<Eigen::Matrix2cd> &paulimatrices, size_t sites)
/*! Builds a string of random pauli matrix MPO's
 *      P = Π  O_i
 * where Π is the product over all sites, and O_i is one of {S, I} on site i.
 * S is the sum of pauli matrices s1,s2,s3... , and where I is an identity matrix of the same size
 *
 *            | s1  0   0  .  |
 * S   =      | 0   s2  0  .  |
 *            | 0   0  s3  .  |
 *            | .   .   . ... |
 *
 *            | id  0   0  .  |
 * I   =      | 0   id  0  .  |
 *            | 0   0  id  .  |
 *            | .   .   . ... |
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
    if(paulimatrices.empty()) throw std::runtime_error("List of pauli matrices is empty");
    long                     num_paulis = static_cast<long>(paulimatrices.size());
    long                     spin_dim   = 2;
    auto                     I          = Eigen::MatrixXcd::Identity(spin_dim, spin_dim);
    Eigen::array<long, 4>    extent4    = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2>    extent2    = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<Scalar, 4> MPO_S(num_paulis, num_paulis, spin_dim, spin_dim);
    Eigen::Tensor<Scalar, 4> MPO_I(num_paulis, num_paulis, spin_dim, spin_dim);
    MPO_S.setZero();
    MPO_I.setZero();
    for(long diag_pos = 0; diag_pos < num_paulis; diag_pos++) {
        MPO_S.slice(Eigen::array<long, 4>{diag_pos, diag_pos, 0, 0}, extent4).reshape(extent2) =
            Textra::MatrixTensorMap(paulimatrices[static_cast<size_t>(diag_pos)]);
        MPO_I.slice(Eigen::array<long, 4>{diag_pos, diag_pos, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(I);
    }

    // Push in an even number of operators
    // This is so that we get a 50% chance of applying a gate.
    std::vector<int> binary(sites, -1);
    int              sum = 0;
    while(true) {
        binary[rnd::uniform_integer_box<size_t>(0, sites-1)] *= -1;
        sum = std::accumulate(binary.begin(), binary.end(), 0);
        if((num::mod<size_t>(sites, 2) == 0 and sum == 0) or (num::mod<size_t>(sites, 2) == 1 and sum == 1)) break;
    }
    if(binary.size() != sites) throw std::logic_error("Size mismatch");
    // Generate the list
    std::list<Eigen::Tensor<Scalar, 4>> mpos;
    std::vector<std::string> mpos_str;
    for(auto &val : binary) {
        if(val < 0){
            mpos.push_back(MPO_S);
            mpos_str.emplace_back("S");
        }else {
            mpos.push_back(MPO_I);
            mpos_str.emplace_back("I");
        }
    }
    tools::log->warn("Generated random pauli MPO string: {}",mpos_str);
    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, num_paulis); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, num_paulis); // The right edge
    Ledge(0, 0, 0) = 1.0 / static_cast<double>(num_paulis);
    Ledge(0, 0, 1) = 1.0 / static_cast<double>(num_paulis);
    Redge(0, 0, 0) = 1;
    Redge(0, 0, 1) = 1;
    return std::make_tuple(mpos, Ledge, Redge);
}




std::tuple<std::list<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
qm::mpo::random_pauli_mpos(const std::vector<Eigen::Matrix2cd> &paulimatrices, const std::vector<double> & uniform_dist_widths, size_t sites)
/*! Builds a set of MPO's used for randomizing a state  pauli matrix MPO's with random weights picked from a uniform distribution
 *      P = Π  O_i
 * where Π is the product over all sites, and O_i is the MPO sum of pauli matrices with random weights.
 *
 *            | c0*s0   0       0     .   |
 * O_i =      | 0       c1*s1   0     .   |
 *            | 0       0       c2*s2 .   |
 *            | .       .       .     ... |
 *  Here s_i are 2x2 pauli matrices (including identity) and
 *  the weight coefficients c_i are random real numbers drawn from a uniform distribution U(-w,w).
 *
 *        2
 *        |
 *    0---O---1
 *        |
 *        3
 *
 */
{
    if(paulimatrices.empty()) throw std::runtime_error("List of pauli matrices is empty");
    if (paulimatrices.size() != uniform_dist_widths.size()) throw std::runtime_error("List size mismatch: paulimatrices and uniform_dist_widths");
    long                     num_paulis = static_cast<long>(paulimatrices.size());
    long                     spin_dim   = 2;
    auto                     I          = Eigen::MatrixXcd::Identity(spin_dim, spin_dim);
    Eigen::array<long, 4>    extent4    = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2>    extent2    = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */

    std::list<Eigen::Tensor<Scalar, 4>> mpos;
    for(size_t site = 0; site < sites; site++){
        Eigen::Tensor<Scalar, 4> MPO_S(num_paulis, num_paulis, spin_dim, spin_dim);
        MPO_S.setZero();
        for(long idx = 0; idx < num_paulis; idx++) {
            auto uidx = static_cast<size_t>(idx);
            auto coeff   = 1 + rnd::uniform_double_box(uniform_dist_widths[uidx]);
            auto offset4 = Eigen::array<long, 4>{idx, idx, 0, 0};
            const auto & pauli =  paulimatrices[uidx];
            MPO_S.slice(offset4, extent4).reshape(extent2) = Textra::MatrixToTensor(coeff * pauli);
        }
        mpos.emplace_back(MPO_S);
    }

    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, num_paulis); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, num_paulis); // The right edge
    Ledge(0, 0, 0) = 1.0 / static_cast<double>(num_paulis);
    Ledge(0, 0, 1) = 1.0 / static_cast<double>(num_paulis);
    Redge(0, 0, 0) = 1;
    Redge(0, 0, 1) = 1;
    return std::make_tuple(mpos, Ledge, Redge);
}
