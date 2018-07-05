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

    namespace SpinOneHalf {
        inline Matrix2cd sx = (Matrix2cd() << 0.0, 1.0,
                                              1.0, 0.0).finished();
        inline Matrix2cd sy = (Matrix2cd() << 0.0, -1.0i,
                                              1.0i, 0.0).finished();
        inline Matrix2cd sz = (Matrix2cd() << 1.0, 0.0,
                                              0.0, -1.0).finished();
        inline Matrix2cd I  = (Matrix2cd() << 1.0, 0.0,
                                              0.0, 1.0).finished();
    }

    namespace SpinOne{
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


    }

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

}
#endif //NMSPC_QUANTUM_MECHANICS_H

