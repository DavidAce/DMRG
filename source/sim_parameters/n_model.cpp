//
// Created by david on 4/25/17.
//

#include "n_model.h"
#include <iomanip>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;

namespace Model{

    const Matrix2cd sx(){
        Matrix2cd sx;
        sx << 0.0, 1.0,
              1.0, 0.0  ;
        return sx;
    }

    const Matrix2cd sy(){
        Matrix2cd sy;
        sy <<  0.0, -1.0i,
               1.0i, 0.0;
        return sy;
    }


    const Matrix2cd sz() {
        Matrix2cd sz;
        sz << 1.0, 0.0,
              0.0,-1.0;
        return sz;
    }

    const Matrix2cd I(){
        return Matrix2cd::Identity();
    }

    std::vector<MatrixXcd>  gen_manybody_spin(const Matrix2cd &s, int sites){
        std::vector<MatrixXcd> S;
        MatrixXcd tmp;
        S.clear();
        for (int i = 0; i < sites; i++) {
            tmp = i == 0 ? s : I();
            for (int j = 1; j < sites; j++) {
                tmp = kroneckerProduct(tmp, i == j ? s : I()).eval();
            }
            S.emplace_back(tmp);
        }
        return S;
    }


    /*! From Zapp, K. Matrix Product States for Lattice Gauge Theories. (2015).
     * and
     * Müller-hermes, B. V. A. Tensor-Network-Methods for Simulating Infinite 1-dimensional Quantum-Many-Body Systems. 1–68 (2010).
     */
    double get_exact_energy(){
        return (-1.0/M_PI/2.0) * Math::compute_integral([](double x){return sqrt(1.0 + g*g - 2.0*g*cos(x));},{-M_PI, M_PI});
    }


    MatrixXcd h(int sites, int position){
        int i = Math::mod(position   ,sites);
        int j = Math::mod(position +1,sites);
        if(spins_must_be_generated){
            SX = gen_manybody_spin(sx(),sites);
            SY = gen_manybody_spin(sy(),sites);
            SZ = gen_manybody_spin(sz(),sites);
            spins_must_be_generated = false;
        }
        return 0.5*(-J * SZ[i] * SZ[j] - 0.5*g*(SX[i] + SX[j]));
    }

    Matrix<Scalar,Dynamic,Dynamic> H(int sites) {
        MatrixXcd hi = MatrixXcd::Zero((long)pow(2,sites), (long)pow(2,sites));
        for (int position = 0; position < sites; position++) {
            hi  += h(sites,position);
        }
        return hi.real();
    }



    Textra::Tensor<4,Scalar> TimeEvolution_4th_order(const double delta, int sites) {
        double delta1 = delta /(4.0-pow(4.0,1.0/3.0));
        double delta2 = delta1;
        double delta3 = delta - 2.0*delta1 - 2.0*delta2;
        return Matrix_to_Tensor<4,Scalar>(U(delta1,sites)*U(delta2,sites)*U(delta3,sites)*U(delta2,sites)*U(delta1,sites), array4{2,2,2,2});

    }
    Textra::Tensor<4,Scalar> TimeEvolution_2nd_order(const double delta, int sites) {
        return Matrix_to_Tensor<4,Scalar>( ((-delta*h(sites,0)).exp() * (-delta*h(sites,1)).exp()).real().eval(), array4{2,2,2,2});
    }
    Textra::Tensor<4,Scalar> TimeEvolution_1st_order(const double delta, int sites) {
        return Matrix_to_Tensor<4,Scalar>(U(delta,sites), array4{2,2,2,2});
    }
    Matrix<Scalar,Dynamic,Dynamic> U(double delta, int sites){
        return ((-delta*h(sites,0)/2.0).exp() * (-delta*h(sites,1)).exp() * (-delta * h(sites,0)/2.0).exp()).real();
    }

    Textra::Tensor<4,Scalar> M() {
        //Returns the H_two as A rank 4 MPO. Notation following Schollwöck (2010)
        MatrixXcd W(6,6);
        W.setZero();
        W.block(0,0,2,2) = I();
        W.block(2,0,2,2) = sz();
        W.block(4,0,2,2) = -g*sx();
        W.block(4,2,2,2) = -J*sz();
        W.block(4,4,2,2) = I();
        return (Matrix_to_Tensor<2,Scalar> (W.real(), {6,6})).reshape(array4{2,3,2,3}).shuffle(array4{3,1,2,0});

    }

    Textra::Tensor<6,Scalar> MM() {
        //Returns a 2-site Hamitlonian MPO of rank 6. Notation following Schollwöck (2010)
        return M().contract(M(), idx<1>({1},{0})).shuffle(Textra::array<6>{0,3,1,4,2,5});

    }
}