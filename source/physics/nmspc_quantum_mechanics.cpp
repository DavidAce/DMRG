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


Eigen::MatrixXcd qm::gen_embedded_spin_operator(const Eigen::MatrixXcd &s, size_t at, size_t sites, bool swap)
/*
 * Returns a spin operator embedded in a larger Hilbert space. For instance, if at == 1 and sites == 4:
 *
 *   σ¹ = i ⊗ σ ⊗ i ⊗ i
 *
 * where each element is a dxd matrix, resulting in a d^4 * d^4 matrix.

 * Note that if this matrix is converted to a rank-8 tensor, the indexing goes like:
 *
 @verbatim
        3 2 1 0
        | | | |
       [  σ¹  ]
       | | | |
       7 6 5 4
 @endverbatim

 * whereas you would normally want left-to-right indexing:
 *
 @verbatim
        0 1 2 3
        | | | |
       [  σ¹  ]
       | | | |
       4 5 6 7
 @endverbatim

 * So don't forget to set "swap = true" if you intend to use the result as a tensor.
 */

{
    if(at >= sites) throw std::logic_error("Expected at < sites. Got [at = " + std::to_string(at) + "] [sites = " + std::to_string(sites) + "]");
    MatrixXcd id = MatrixXcd::Identity(s.rows(),s.cols());
    MatrixXcd result = at == 0 ? s : id;
    if(swap)
        for(size_t site = 1; site < sites; site++) result = kroneckerProduct(site == at ? s : id, result).eval(); // .eval() is required to avoid aliasing!!
    else
        for(size_t site = 1; site < sites; site++) result = kroneckerProduct(result, site == at ? s : id).eval(); // .eval() is required to avoid aliasing!!
    return result;
}


std::vector<Eigen::MatrixXcd> qm::gen_manybody_spins(const Eigen::MatrixXcd &s, int sites, bool swap) {
    std::vector<MatrixXcd> S;
    for(int site = 0; site < sites; site++) S.emplace_back(qm::gen_embedded_spin_operator(s,site,sites,swap));
    return S;
}

namespace qm::spinHalf {

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


    Eigen::MatrixXcd gen_embedded_spin_operator(const Eigen::Matrix2cd &s, size_t at, size_t sites,bool swap){
        // Don't forget to set "swap = true" if you intend to use the result as a tensor.
        return qm::gen_embedded_spin_operator(s,at,sites,swap);
    }


    std::vector<Eigen::Matrix4cd> gen_twobody_spins(const Eigen::Matrix2cd &s, bool swap)
    // Returns a pair of two-body 4x4 spin operators for embedded in a two-site Hilbert space:
    //        (σ ⊗ i, i ⊗ σ)
    // where σ is a 2x2 (pauli) matrix and i is the 2x2 identity matrix.
    // Don't forget to set "swap = true" if you intend to use the result as a tensor.
    {
        std::vector<Matrix4cd> S;
        for (size_t site = 0; site < 2; site++)
            S.emplace_back(qm::gen_embedded_spin_operator(s,site,2,swap));
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

    std::vector<Eigen::MatrixXcd> gen_twobody_spins(const Eigen::Matrix3cd &s, bool swap)
    // Returns a pair of two-body 4x4 spin operators for embedded in a two-site Hilbert space:
    //        (σ ⊗ i, i ⊗ σ)
    // where σ is a 3x3 (pauli) matrix and i is the 3x3 identity matrix.
    // So don't forget to set "swap = true" if you intend to use the result as a tensor.
{
        return qm::gen_manybody_spins(s,2,swap);
    }
}

namespace qm::timeEvolution {

    std::vector<Eigen::Tensor<Scalar,2>> Suzuki_Trotter_1st_order(cplx t, const Eigen::Tensor<Scalar,2> &h_evn, const Eigen::Tensor<Scalar,2> &h_odd) {
        auto h_evn_matrix = Textra::TensorMatrixMap(h_evn);
        auto h_odd_matrix = Textra::TensorMatrixMap(h_odd);
        return
            {
                Textra::MatrixToTensor((t * h_evn_matrix).exp()),
                Textra::MatrixToTensor((t * h_odd_matrix).exp())
            };
    }

    std::vector<Eigen::Tensor<Scalar,2>> Suzuki_Trotter_2nd_order(cplx t, const Eigen::Tensor<Scalar,2> &h_evn, const Eigen::Tensor<Scalar,2> &h_odd) {
        auto h_evn_matrix = Textra::TensorMatrixMap(h_evn);
        auto h_odd_matrix = Textra::TensorMatrixMap(h_odd);
        return
        {
            Textra::MatrixToTensor((t * h_evn_matrix / 2.0).exp()),
            Textra::MatrixToTensor((t * h_odd_matrix).exp()),
            Textra::MatrixToTensor((t * h_evn_matrix / 2.0).exp())
        };
    }

    std::vector<Eigen::Tensor<Scalar,2>> Suzuki_Trotter_4th_order(cplx t, const Eigen::Tensor<Scalar,2> &h_evn, const Eigen::Tensor<Scalar,2> &h_odd)
    /*!
     * Implementation based on
     * Janke, W., & Sauer, T. (1992).
     * Properties of higher-order Trotter formulas.
     * Physics Letters A, 165(3), 199–205.
     * https://doi.org/10.1016/0375-9601(92)90035-K
     *
     */
    {
        auto h_evn_matrix = Textra::TensorMatrixMap(h_evn);
        auto h_odd_matrix = Textra::TensorMatrixMap(h_odd);
        double cbrt2 = pow(2.0, 1.0 / 3.0);
        double beta1 = 1.0 / (2.0 - cbrt2);
        double beta2 = -cbrt2 * beta1;
        double alph1 = 0.5 * beta1;
        double alph2 = (1.0 - cbrt2) / 2.0 * beta1;

        std::vector<Eigen::Tensor<Scalar,2>> temp;
        temp.emplace_back(Textra::MatrixToTensor((alph1 * t * h_evn_matrix).exp()));
        temp.emplace_back(Textra::MatrixToTensor((beta1 * t * h_odd_matrix).exp()));
        temp.emplace_back(Textra::MatrixToTensor((alph2 * t * h_evn_matrix).exp()));
        temp.emplace_back(Textra::MatrixToTensor((beta2 * t * h_odd_matrix).exp()));
        temp.emplace_back(Textra::MatrixToTensor((alph2 * t * h_evn_matrix).exp()));
        temp.emplace_back(Textra::MatrixToTensor((beta1 * t * h_odd_matrix).exp()));
        temp.emplace_back(Textra::MatrixToTensor((alph1 * t * h_evn_matrix).exp()));
        return temp;
    }

    std::vector<Eigen::Tensor<Scalar, 2>> get_twosite_time_evolution_operators(std::complex<double> t, size_t susuki_trotter_order,
                                                                                  const Eigen::Tensor<Scalar,2> &h_evn, const Eigen::Tensor<Scalar,2> &h_odd)
    /*! Returns a set of 2-site unitary gates, using Suzuki Trotter decomposition to order 1, 2 or 3.
     * These gates need to be applied to the MPS one at a time with a swap in between.
     */
    {
        switch(susuki_trotter_order) {
            case 1:  return Suzuki_Trotter_1st_order(t, h_evn, h_odd);
            case 2:  return Suzuki_Trotter_2nd_order(t, h_evn, h_odd);
            case 4:  return Suzuki_Trotter_4th_order(t, h_evn, h_odd);
            default: return Suzuki_Trotter_2nd_order(t, h_evn, h_odd);
        }
    }



    std::vector<Eigen::Tensor<Scalar, 2>> compute_G(const cplx a, size_t susuki_trotter_order, const Eigen::Tensor<Scalar,2> &h_evn, const Eigen::Tensor<Scalar,2> &h_odd)
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
        return get_twosite_time_evolution_operators(a, susuki_trotter_order, h_evn, h_odd);
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
    constexpr bool kroneckerSwap = false;
    auto SZ  = qm::spinHalf::gen_twobody_spins(qm::spinHalf::sz,kroneckerSwap); // We use these as matrices
    auto SP  = qm::spinHalf::gen_twobody_spins(qm::spinHalf::sp,kroneckerSwap); // We use these as matrices
    auto SM  = qm::spinHalf::gen_twobody_spins(qm::spinHalf::sm,kroneckerSwap); // We use these as matrices
    auto ID  = qm::spinHalf::gen_twobody_spins(qm::spinHalf::id,kroneckerSwap); // We use these as matrices
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


        if constexpr (kroneckerSwap){
            // Here the kronecker already has index pattern left-to-right and there is no need to shuffle

            //         0               0      1
            //         |               |      |
            //   [ exp(-ifH) ]  ==  [ exp(-ifH) ]
            //        |               |      |
            //        1               2      3

            unitaries[idx] = Textra::MatrixToTensor2((imn * fmix * H).exp());
        }else{
            // Here we shuffle to get the correct underlying index pattern: Sites are contracted left-to right, but
            // the kronecker product that generated two-site gates above has indexed right-to-left
            //         0                   1      0              0      1               0
            //         |                   |      |              |      |               |
            //   [ exp(-ifH) ]  --->    [ exp(-ifH) ]   --->  [ exp(-ifH) ]  --->  [ exp(-ifH) ]
            //        |                   |      |              |      |                |
            //        1                   3      2              2      3                1
            Eigen::Tensor<Scalar,2> H_shuffled = Textra::MatrixTensorMap(H,2,2,2,2).shuffle(Textra::array4{1,0,3,2}).reshape(Textra::array2{4,4});
            Eigen::MatrixXcd expifH = (imn * fmix * Textra::TensorMatrixMap(H_shuffled)).exp();
            unitaries[idx] = Textra::MatrixTensorMap(expifH);
        }
    }
    // Sanity check
    for(auto && u : unitaries) if(not Textra::TensorMatrixMap(u).isUnitary()) throw std::logic_error("u is not unitary!");
    return unitaries;
}

std::vector<Eigen::Tensor<Scalar,2>> qm::lbit::get_twosite_time_evolution_operators(size_t sites, cplx delta_t, const std::vector<Eigen::Tensor<Scalar,2>> &twosite_hams){
    // In l-bit systems we are aldready in a diagonal basis, so h_{j,j+1} and h_{j+1,j+2} commute. Therefore we can immediately use the relation
    //      exp(-i*dt *[h_{j,j+1} + h_{j+1,j+2} + ... + h_{L-2, L-1}]) =  exp(-i*dt [h_{j,j+1}]) * exp(-i*dt*[h_{j+1,j+2}]) * ... * exp(-i*dt*[h_{L-2, L-1}])
    // without passing through the Suzuki-Trotter decomposition.

    // Here we expect "twosite_hams" to contain operators  h_{j,j+1} + h_{j+1,j+2} + ... + h_{L-2, L-1}.

    if(twosite_hams.size() != sites-1)
        throw std::logic_error(fmt::format("Wrong number of twosite hamiltonians: {}. Expected {}", twosite_hams.size(),sites-1));

    std::vector<Eigen::Tensor<Scalar,2>> time_evolution_operators;
    time_evolution_operators.reserve(sites-1);
    for(auto && h : twosite_hams)
       time_evolution_operators.emplace_back(Textra::MatrixToTensor((delta_t * Textra::TensorMatrixMap(h)).exp()));
    return time_evolution_operators;
}

std::vector<Eigen::Tensor<Scalar,2>> qm::lbit::get_3site_time_evolution_operators(size_t sites, cplx delta_t, const std::vector<Eigen::Tensor<Scalar,2>> &hams_3site){
    // In l-bit systems we are aldready in a diagonal basis, so h_{i,j,k} and h_{l,m,n} commute. Therefore we can immediately use the relation
    // exp(A + B) = exp(A)exp(B)
    // without passing through the Suzuki-Trotter decomposition.

    if(hams_3site.size() != sites-2)
        throw std::logic_error(fmt::format("Wrong number of three-site hamiltonians: {}. Expected {}", hams_3site.size(),sites-2));

    std::vector<Eigen::Tensor<Scalar,2>> time_evolution_operators;
    time_evolution_operators.reserve(sites-1);
    for(auto && h : hams_3site)
        time_evolution_operators.emplace_back(Textra::MatrixToTensor((delta_t * Textra::TensorMatrixMap(h)).exp()));
    return time_evolution_operators;
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
 *      |psi+->  = P |psi>=  1/sqrt(2) (1 +- S) |psi>
 *      Here 1 = outer product of L=sites 2x2 identity matrices, i.e. Kron_(i=0)^(L-1) I_(2x2)
 *      Also S = outer product of L=sites 2x2 pauli matrices, i.e. Kron_(i=0)^(L-1) s_(2x2)
 *      The sign and the factor 1/2 is put into the left edge at the end.
 *
 *                     | I   0  |
 * S   =   1/sqrt(2) * | 0   s  |
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
    Ledge(0, 0, 0) = 1.0/std::sqrt(2);//0.5;
    Ledge(0, 0, 1) = 1.0/std::sqrt(2) * sign;
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
    Ledge(0, 0, 0) = 1.0 / std::sqrt(2);
    Ledge(0, 0, 1) = 1.0 / std::sqrt(2);
    Redge(0, 0, 0) = 1;
    Redge(0, 0, 1) = 1;
    return std::make_tuple(mpos, Ledge, Redge);
}

std::tuple<std::list<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
    qm::mpo::random_pauli_mpos(const std::vector<Eigen::Matrix2cd> &paulimatrices, size_t sites)
/*! Builds a string of random pauli matrix MPO's
 *      P = Π  O_i
 * where Π is the product over all sites, and O_i is one of {S, I} on site i.
 * S is the sum of pauli matrices s0,s1,s2... , and where I is an identity matrix of the same size
 *
 *            | s0  0   0  .  |
 * S   =      | 0   s1  0  .  |
 *            | 0   0  s2  .  |
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
    Ledge(0, 0, 0) = 1.0 / std::sqrt(num_paulis);
    Ledge(0, 0, 1) = 1.0 / std::sqrt(num_paulis);
    Redge(0, 0, 0) = 1;
    Redge(0, 0, 1) = 1;
    return std::make_tuple(mpos, Ledge, Redge);
}



std::tuple<std::list<Eigen::Tensor<Scalar, 4>>, Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 3>>
qm::mpo::sum_of_pauli_mpo(const std::vector<Eigen::Matrix2cd> &paulimatrices, size_t sites, bool shuffle)
/*! Builds a string of MPO's
 *      P = Π  O_i
 * where Π is the product over all sites, and O_i are MPOs with 2x2 (pauli) matrices on the diagonal
 *
 *            | s0  0   0  .  |
 * O_i =      | 0   s1  0  .  |
 *            | 0   0  s2  .  |
 *            | .   .   . ... |
 *
 * The matrices s0, s1, s2 are shuffled if the argument shuffle == true
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
    Eigen::array<long, 4>    extent4    = {1, 1, spin_dim, spin_dim}; /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2>    extent2    = {spin_dim, spin_dim};       /*!< Extent of pauli matrices in a rank-2 tensor */
    std::list<Eigen::Tensor<Scalar, 4>> mpos;
    std::vector<size_t> pauli_idx = num::range<size_t,size_t>(0,num_paulis-1,1);
    for(size_t site = 0; site < sites; site++){
        Eigen::Tensor<Scalar, 4> MPO_S(num_paulis, num_paulis, spin_dim, spin_dim);
        MPO_S.setZero();
        if(shuffle)
            std::shuffle(pauli_idx.begin(),pauli_idx.end(),rnd::internal::rng);
        for(long idx = 0; idx < num_paulis; idx++) {
            auto uidx = static_cast<size_t>(idx);
            const auto & pauli =  paulimatrices[pauli_idx[uidx]];
            auto offset4 = Eigen::array<long, 4>{idx, idx, 0, 0};
            MPO_S.slice(offset4, extent4).reshape(extent2) = Textra::MatrixToTensor(pauli);
        }
        mpos.emplace_back(MPO_S);
    }

    // Create compatible edges
    Eigen::Tensor<Scalar, 3> Ledge(1, 1, num_paulis); // The left  edge
    Eigen::Tensor<Scalar, 3> Redge(1, 1, num_paulis); // The right edge
    Ledge(0, 0, 0) = 1.0 / std::sqrt(num_paulis);
    Ledge(0, 0, 1) = 1.0 / std::sqrt(num_paulis);
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
    Ledge(0, 0, 0) = 1.0 / std::sqrt(num_paulis);
    Ledge(0, 0, 1) = 1.0 / std::sqrt(num_paulis);
    Redge(0, 0, 0) = 1;
    Redge(0, 0, 1) = 1;
    return std::make_tuple(mpos, Ledge, Redge);
}
