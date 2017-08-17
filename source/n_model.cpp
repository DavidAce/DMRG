//
// Created by david on 4/25/17.
//

#include "n_model.h"
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;

namespace Model{

    Matrix2cd gen_sx(){
        Matrix2cd sx;
        sx << 0.0, 1.0,
             1.0, 0.0  ;
        return sx;
    }

    Matrix2cd gen_sy(){
        Matrix2cd sy;
        sy <<  0.0, -1.0i,
                1.0i, 0.0;
        return sy;
    }

    Matrix2cd gen_sz(){
        Matrix2cd sz;
        sz << 1.0, 0.0,
                0.0, -1.0;
        return sz;
    }
    Matrix2cd gen_I(){
        return Matrix2cd::Identity();
    }

    std::vector<MatrixXcd>  generate_manybody_spin(const int sites, Matrix2cd &s){
        std::vector<MatrixXcd> S;
        MatrixXcd tmp;
        S.clear();
        for (int i = 0; i < sites; i++) {
            tmp = i == 0 ? s : I;
            for (int j = 1; j < sites; j++) {
                tmp = kroneckerProduct(tmp, i == j ? s : I).eval();
            }
            S.emplace_back(tmp);
        }
        return S;
    }

//    void generate_spins(const int sites) {
////        sx << 0.0, 1.0,
////              1.0, 0.0;
////        sy << 0.0, -1.0i,
////              1.0i, 0.0;
////        sz << 1.0, 0.0,
////              0.0, -1.0;
//
//        I.setIdentity();
//        MatrixXcd X, Y, Z;
//        SX.clear();
//        SY.clear();
//        SZ.clear();
//
//        for (int i = 0; i < sites; i++) {
//            X = i == 0 ? sx : I;
//            Y = i == 0 ? sy : I;
//            Z = i == 0 ? sz : I;
//            for (int j = 1; j < sites; j++) {
//                X = kroneckerProduct(X, i == j ? sx : I).eval();
//                Y = kroneckerProduct(Y, i == j ? sy : I).eval();
//                Z = kroneckerProduct(Z, i == j ? sz : I).eval();
//            }
//            SX.emplace_back(X);
//            SY.emplace_back(Y);
//            SZ.emplace_back(Z);
//        }
//    }

    double get_exact_energy(){
        return Math::compute_integral([](double x){return -2.0 * sqrt(1.0 + g*g - 2.0*g*cos(x))/M_PI/2.0 ;},{0, M_PI});
    }


    MatrixType Hamiltonian(const int sites) {
//        generate_spins(sites);
        MatrixType H = MatrixType::Zero(Math::ipow(2, sites), Math::ipow(2, sites));

        for (int i = 0; i < sites; i++) {
            H += 0.5 * (J * SZ[i] * SZ[Math::mod(i + 1, sites)] + g * SX[i]).real();
        }
        return H;
    }


    MatrixType Hamiltonian2(const int sites) {
        MatrixType H = MatrixType::Zero(Math::ipow(2, sites), Math::ipow(2, sites));
        H += 0.5 * (J * SZ[0] * SZ[1] + g * SX[0]).real() / 2.0;
        H += 0.5 * (J * SZ[1] * SZ[0] + g * SX[1]).real() ;
        H += 0.5 * (J * SZ[0] * SZ[1] + g * SX[0]).real() / 2.0;
        return H;
    }

    Tensor4 TimeEvolution_4th_order(const int sites, const double delta) {
        H1 = MatrixType::Zero(Math::ipow(2, sites), Math::ipow(2, sites));
        H2 = MatrixType::Zero(Math::ipow(2, sites), Math::ipow(2, sites));
        H1 += 0.5 * (J * SZ[0] * SZ[1] + g * SX[0]).real();
        H2 += 0.5 * (J * SZ[1] * SZ[0] + g * SX[1]).real();
        double delta1 = delta /(4.0-pow(4,1.0/3.0));
        double delta2 = delta1;
        double delta3 = delta - 2.0*delta1 - 2.0*delta2;
        double factor = -2.0*(delta1 + delta2+ 0.5*delta3);
        return matrix_to_tensor<4>((factor * (H1/2.0 + H2 + H1/2.0)).exp().eval(), array4{2,2,2,2});

    }

//    Tensor4 TimeEvolution_4th_order(const int sites, const double delta) {
//        generate_spins(sites);
//
//        MatrixType H1 = MatrixType::Zero(Math::ipow(2, sites), Math::ipow(2, sites));
//        MatrixType H2 = MatrixType::Zero(Math::ipow(2, sites), Math::ipow(2, sites));
//        H1 += 0.5 * (J * SZ[0] * SZ[1] + g * SX[0]).real();
//        H2 += 0.5 * (J * SZ[1] * SZ[0] + g * SX[1]).real();
//        double delta1 = delta /(4.0-pow(4,1.0/3.0));
//        double delta2 = delta1;
//        double delta3 = delta - 2.0*delta1 - 2.0*delta2;
//
//        return matrix_to_tensor<4>(U(H1,H2,delta1)*U(H1,H2,delta2)*U(H1,H2,delta3)*U(H1,H2,delta2)*U(H1,H2,delta1), array4{2,2,2,2});
//
//    }



    MatrixType U(const double delta){
        return (-delta*H1/2.0).exp() * (-delta*H2).exp() * (-delta * H1/2.0).exp();
    }

    Tensor4 W(const int sites) {
        //Returns the Hamiltonian as A rank 4 MPO. Notation following Schollwöck (2010)
//        generate_spins(sites);
        MatrixXcd W(6,6);
        W.setZero();
        W.block(0,0,2,2) = I;
        W.block(2,0,2,2) = sz;
        W.block(4,0,2,2) = g*sx;
        W.block(4,2,2,2) = J*sz;
        W.block(4,4,2,2) = I;
        return matrix_to_tensor2(W.real()).reshape(array4{2,3,2,3}).shuffle(array4{3,1,2,0});

    }
}