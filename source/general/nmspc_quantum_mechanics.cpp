//
// Created by david on 2019-03-20.
//
#include <Eigen/Core>
#include <vector>
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
    Matrix2cd sx = (Matrix2cd() << 0.0, 1.0,
            1.0, 0.0).finished();
    Matrix2cd sy = (Matrix2cd() << 0.0, imn,
            imp, 0.0).finished();
    Matrix2cd sz = (Matrix2cd() << 1.0, 0.0,
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




namespace qm{
    namespace timeEvolution{

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
                tensor_vec.emplace_back(Textra::Matrix_to_Tensor(m, 2,2,2,2));
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





}