//
// Created by david on 2018-04-17.
//

#ifndef NMSPC_QUANTUM_MECHANICS_H
#define NMSPC_QUANTUM_MECHANICS_H
#include <Eigen/Core>
#include <complex>
#include <general/nmspc_tensor_extra.h>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std::complex_literals;
namespace qm{


    using namespace Eigen;
    using Scalar = std::complex<double>;
    inline std::vector<MatrixXcd> gen_manybody_spin(const MatrixXcd &s, int sites) {
        std::vector<MatrixXcd> S;
        MatrixXcd I = MatrixXcd::Identity(s.rows(),s.cols());
        MatrixXcd tmp;
        for (int i = 0; i < sites; i++) {
            tmp = i == 0 ? s : I;
            for (int j = 1; j < sites; j++) {
                tmp = kroneckerProduct(tmp, i == j ? s : I).eval();
            }
            S.emplace_back(tmp);
        }
        return S;
    }


    namespace spinOneHalf {
        inline Matrix2cd sx = (Matrix2cd() << 0.0, 1.0,
                                              1.0, 0.0).finished();
        inline Matrix2cd sy = (Matrix2cd() << 0.0, -1.0i,
                                              1.0i, 0.0).finished();
        inline Matrix2cd sz = (Matrix2cd() << 1.0, 0.0,
                                              0.0, -1.0).finished();
        inline Matrix2cd I  = (Matrix2cd() << 1.0, 0.0,
                                              0.0, 1.0).finished();
        inline std::vector<Eigen::MatrixXcd> SX;
        inline std::vector<Eigen::MatrixXcd> SY;
        inline std::vector<Eigen::MatrixXcd> SZ;
        inline std::vector<Eigen::MatrixXcd> II;

    }

    namespace spinOne{
        inline Matrix3cd sx = (Matrix3cd() << 0.0, 1.0, 0.0,
                                              1.0, 0.0, 1.0,
                                              0.0, 1.0, 0.0).finished();
        inline Matrix3cd sy = (Matrix3cd() << 0.0 , -1.0i, 0.0,
                                              1.0i,  0.0 ,-1.0i,
                                              0.0 ,  1.0i, 0.0).finished();
        inline Matrix3cd sz = (Matrix3cd() << 1.0, 0.0, 0.0,
                                              0.0, 0.0, 0.0,
                                              0.0, 0.0,-1.0).finished();
        inline Matrix3cd I = (Matrix3cd()  << 1.0, 0.0, 0.0,
                                              0.0, 1.0, 0.0,
                                              0.0, 0.0, 1.0).finished();

        inline std::vector<Eigen::MatrixXcd> SX;
        inline std::vector<Eigen::MatrixXcd> SY;
        inline std::vector<Eigen::MatrixXcd> SZ;
        inline std::vector<Eigen::MatrixXcd> II;
    }

    namespace timeEvolution{

        inline std::vector<Textra::MatrixType<Scalar>> Suzuki_Trotter_1st_order(Scalar t, Eigen::MatrixXcd &h_evn, Eigen::MatrixXcd &h_odd){
            return {(t*h_evn).exp(),
                    (t*h_odd).exp() };
        }

        inline std::vector<Textra::MatrixType<Scalar>> Suzuki_Trotter_2nd_order(Scalar t, Eigen::MatrixXcd &h_evn, Eigen::MatrixXcd &h_odd){
            return {(t*h_evn/2.0).exp(),
                    (t*h_odd).exp(),
                    (t*h_evn/2.0).exp()};
        }


        inline std::vector<Textra::MatrixType<Scalar>> Suzuki_Trotter_4th_order(Scalar t, Eigen::MatrixXcd &h_evn, Eigen::MatrixXcd &h_odd)
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

            std::vector<Textra::MatrixType<Scalar>> temp;

            temp.emplace_back( (alph1 *  t*h_evn).exp() );
            temp.emplace_back( (beta1 *  t*h_odd).exp() );
            temp.emplace_back( (alph2 *  t*h_evn).exp() );
            temp.emplace_back( (beta2 *  t*h_odd).exp() );
            temp.emplace_back( (alph2 *  t*h_evn).exp() );
            temp.emplace_back( (beta1 *  t*h_odd).exp() );
            temp.emplace_back( (alph1 *  t*h_evn).exp() );
            return temp;
        }


        inline std::vector<Eigen::Tensor<Scalar,4>> get_2site_evolution_gates(const Scalar t,int susuki_trotter_order, Eigen::MatrixXcd &h_evn, Eigen::MatrixXcd &h_odd)
        /*! Returns a set of 2-site unitary gates, using Suzuki Trotter decomposition to order 1, 2 or 3.
         * These gates need to be applied to the MPS one at a time with a swap in between.
         */
        {
            std::vector<Textra::MatrixType<Scalar>> matrix_vec;
            switch (susuki_trotter_order) {
                case 1:  matrix_vec = Suzuki_Trotter_1st_order(t, h_evn, h_odd);break;
                case 2:  matrix_vec = Suzuki_Trotter_2nd_order(t, h_evn, h_odd);break;
                case 4:  matrix_vec = Suzuki_Trotter_4th_order(t, h_evn, h_odd);break;
                default: matrix_vec = Suzuki_Trotter_2nd_order(t, h_evn, h_odd);break;
            }
            std::vector<Eigen::Tensor<Scalar ,4>> tensor_vec;
            for(auto &m : matrix_vec){
                tensor_vec.emplace_back(Textra::Matrix_to_Tensor(m, 2,2,2,2));
            }
            return tensor_vec;
        }


//        inline void update_evolution_step_size(const Scalar dt, const int susuki_trotter_order, Eigen::MatrixXcd &h_evn, Eigen::MatrixXcd &h_odd){
//            /*! Returns a set of 2-site unitary gates for the time evolution operator. */
//            step_size = std::abs(dt);
//            U = get_2site_evolution_gates(dt, susuki_trotter_order, h_evn,h_odd);
//        }


        inline std::vector<Eigen::Tensor<Scalar,4>> compute_G(Scalar a, int susuki_trotter_order, Eigen::MatrixXcd &h_evn, Eigen::MatrixXcd &h_odd)
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
#endif //NMSPC_QUANTUM_MECHANICS_H

