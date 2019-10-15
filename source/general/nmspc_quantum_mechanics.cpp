//
// Created by david on 2019-03-20.
//
#include <Eigen/Core>
#include <vector>
#include <complex.h>
#undef I
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
#include "nmspc_quantum_mechanics.h"
#include "nmspc_tensor_extra.h"



using namespace Eigen;

std::vector<Eigen::MatrixXcd> qm::gen_manybody_spin(const Eigen::MatrixXcd &s, int sites) {
    std::vector<MatrixXcd> S;
    MatrixXcd Id = MatrixXcd::Identity(s.rows(),s.cols());
    MatrixXcd tmp;
    for (int idx = 0; idx < sites; idx++) {
        tmp = idx == 0 ? s : Id;
        for (int j = 1; j < sites; j++) {
            tmp = kroneckerProduct(tmp, idx == j ? s : Id).eval();
        }
        S.emplace_back(tmp);
    }
    return S;
}



namespace qm::spinOneHalf {
    auto imp = std::complex<double>(0.0,1.0);
    auto imn = std::complex<double>(0.0,-1.0);
    Matrix2cd sx = (Matrix2cd() <<
            0.0, 1.0,
            1.0, 0.0).finished();
    Matrix2cd sy = (Matrix2cd() <<
            0.0, imn,
            imp, 0.0).finished();
    Matrix2cd sz = (Matrix2cd() <<
            1.0, 0.0,
            0.0, -1.0).finished();
    Matrix2cd Id  = (Matrix2cd() << 1.0, 0.0,
            0.0, 1.0).finished();

    std::array<Vector2cd,2> sx_eigvecs {(Vector2cd() << 1.0, 1.0).finished()/std::sqrt(2),
                                        (Vector2cd() << 1.0,-1.0).finished()/std::sqrt(2)};

    std::array<Vector2cd,2> sy_eigvecs {(Vector2cd() << 1.0, imp).finished()/std::sqrt(2),
                                        (Vector2cd() << 1.0, imn).finished()/std::sqrt(2)};

    std::array<Vector2cd,2> sz_eigvecs {(Vector2cd() << 1.0, 0.0).finished()/std::sqrt(2),
                                        (Vector2cd() << 0.0, 1.0).finished()/std::sqrt(2)};

    std::vector<Eigen::MatrixXcd> SX;
    std::vector<Eigen::MatrixXcd> SY;
    std::vector<Eigen::MatrixXcd> SZ;
    std::vector<Eigen::MatrixXcd> II;

}


namespace qm::SpinOne{
    auto imp = std::complex<double>(0.0,1.0);
    auto imn = std::complex<double>(0.0,-1.0);
    Matrix3cd sx = (Matrix3cd() <<  0.0, 1.0, 0.0,
                                    1.0, 0.0, 1.0,
                                    0.0, 1.0, 0.0).finished();
    Matrix3cd sy = (Matrix3cd() <<  0.0 , imn, 0.0,
                                    imp,  0.0 ,imn,
                                    0.0 ,  imp, 0.0).finished();
    Matrix3cd sz = (Matrix3cd() <<  1.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0,
                                    0.0, 0.0,-1.0).finished();
    Matrix3cd Id = (Matrix3cd()  << 1.0, 0.0, 0.0,
                                    0.0, 1.0, 0.0,
                                    0.0, 0.0, 1.0).finished();

    std::vector<Eigen::MatrixXcd> SX;
    std::vector<Eigen::MatrixXcd> SY;
    std::vector<Eigen::MatrixXcd> SZ;
    std::vector<Eigen::MatrixXcd> II;
}




namespace qm::timeEvolution{

    std::vector<Eigen::MatrixXcd> Suzuki_Trotter_1st_order(const std::complex<double> t, const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd){
        return {(t*h_evn).exp(),
                (t*h_odd).exp() };
    }

    std::vector<Eigen::MatrixXcd> Suzuki_Trotter_2nd_order(const std::complex<double> t, const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd){
        return {(t*h_evn/2.0).exp(),
                (t*h_odd).exp(),
                (t*h_evn/2.0).exp()};
    }


    std::vector<Eigen::MatrixXcd> Suzuki_Trotter_4th_order(const std::complex<double> t, const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd)
    /*!
     * Implementation based on
     * Janke, W., & Sauer, T. (1992).
     * Properties of higher-order Trotter formulas.
     * Physics Letters A, 165(3), 199–205.
     * https://doi.org/10.1016/0375-9601(92)90035-K
     *
     */
    {
        double cbrt2 = pow(2.0,1.0/3.0);
        double beta1 = 1.0/(2.0 - cbrt2);
        double beta2 = - cbrt2 *beta1;
        double alph1 = 0.5*beta1;
        double alph2 = (1.0 - cbrt2)/2.0 * beta1;

        std::vector<Eigen::MatrixXcd> temp;

        temp.emplace_back( (alph1 *  t*h_evn).exp() );
        temp.emplace_back( (beta1 *  t*h_odd).exp() );
        temp.emplace_back( (alph2 *  t*h_evn).exp() );
        temp.emplace_back( (beta2 *  t*h_odd).exp() );
        temp.emplace_back( (alph2 *  t*h_evn).exp() );
        temp.emplace_back( (beta1 *  t*h_odd).exp() );
        temp.emplace_back( (alph1 *  t*h_evn).exp() );
        return temp;
    }


    std::vector<Eigen::Tensor<std::complex<double>,4>> get_2site_evolution_gates(const std::complex<double> t,const int susuki_trotter_order, const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd)
    /*! Returns a set of 2-site unitary gates, using Suzuki Trotter decomposition to order 1, 2 or 3.
     * These gates need to be applied to the MPS one at a time with a swap in between.
     */
    {
        std::vector<Eigen::MatrixXcd> matrix_vec;
        switch (susuki_trotter_order) {
            case 1:  matrix_vec = Suzuki_Trotter_1st_order(t, h_evn, h_odd);break;
            case 2:  matrix_vec = Suzuki_Trotter_2nd_order(t, h_evn, h_odd);break;
            case 4:  matrix_vec = Suzuki_Trotter_4th_order(t, h_evn, h_odd);break;
            default: matrix_vec = Suzuki_Trotter_2nd_order(t, h_evn, h_odd);break;
        }
        std::vector<Eigen::Tensor<std::complex<double> ,4>> tensor_vec;
        for(auto &m : matrix_vec){
            tensor_vec.emplace_back(Textra::MatrixTensorMap(m, 2,2,2,2));
        }
        return tensor_vec;
    }


//        inline void update_evolution_step_size(const std::complex<double> dt, const int susuki_trotter_order, Eigen::MatrixXcd &h_evn, Eigen::MatrixXcd &h_odd){
//            /*! Returns a set of 2-site unitary gates for the time evolution operator. */
//            step_size = std::abs(dt);
//            U = get_2site_evolution_gates(dt, susuki_trotter_order, h_evn,h_odd);
//        }


    std::vector<Eigen::Tensor<std::complex<double>,4>> compute_G(const std::complex<double> a, const int susuki_trotter_order, const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd)
    /*! Returns the moment generating function, or characteristic function (if a is imaginary) for the Hamiltonian as a rank 4 tensor.
    *   G := exp(iaM) or exp(aM), where a is a small parameter and M is an MPO.
    *   Note that G(-a) = G(a)* if  exp(iaM) !
    *
    @verbatim
                0         1
                |         |
                [ exp(aH) ]
                |         |
                2         3
    @endverbatim
    */
    {
        return get_2site_evolution_gates(a, susuki_trotter_order, h_evn,h_odd);
    }

}





std::tuple<
        Eigen::Tensor<qm::mpo::Scalar,4>,
        Eigen::Tensor<qm::mpo::Scalar,3>,
        Eigen::Tensor<qm::mpo::Scalar,3>>
qm::mpo::pauli_mpo(const Eigen::MatrixXcd paulimatrix)
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
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix);

    //Create compatible edges
    Eigen::Tensor<Scalar,3> Ledge(1,1,1); // The left  edge
    Eigen::Tensor<Scalar,3> Redge(1,1,1); // The right edge
    Ledge(0,0,0) = 1;
    Redge(0,0,0) = 1;
    return std::make_tuple(MPO,Ledge,Redge);
}






std::tuple<
        Eigen::Tensor<qm::mpo::Scalar,4>,
        Eigen::Tensor<qm::mpo::Scalar,3>,
        Eigen::Tensor<qm::mpo::Scalar,3>>
qm::mpo::parity_selector_mpo(const Eigen::MatrixXcd paulimatrix, const int sector)
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
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(I);
    MPO.slice(Eigen::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sector * paulimatrix);

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
        std::list<
        Eigen::Tensor<qm::mpo::Scalar,4>>,
        Eigen::Tensor<qm::mpo::Scalar,3>,
        Eigen::Tensor<qm::mpo::Scalar,3>>
qm::mpo::parity_projector_mpos(const Eigen::MatrixXcd paulimatrix, const size_t sites, const int sector)/*! Builds the MPO that projects out the MPS component in a parity sector.
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
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(I);
    MPO.slice(Eigen::array<long, 4>{1, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(paulimatrix);

    std::list<Eigen::Tensor<Scalar,4>> mpos(sites,MPO);
    //Create compatible edges
    Eigen::Tensor<Scalar,3> Ledge(1,1,2); // The left  edge
    Eigen::Tensor<Scalar,3> Redge(1,1,2); // The right edge
    Ledge(0,0,0) = 1.0/2.0;
    Ledge(0,0,1) = 1.0/2.0 * sector;
    Redge(0,0,0) = 1;
    Redge(0,0,1) = 1;

    return std::make_tuple(mpos,Ledge,Redge);
}