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
                0.0, -1.0;
        return sz;
    }

    const Matrix2cd I(){
        return Matrix2cd::Identity();
    }

    std::vector<MatrixXcd>  gen_twosite_spin(Matrix2cd s){
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


    MatrixXd h(int pos){
        int i = Math::mod(pos  ,sites);
        int j = Math::mod(pos+1,sites);
//        return 0.5*(-J * SZ[i] * SZ[j] - g * SX[i]).real();
//        return (-J * SZ[i] * SZ[j] - g * SX[i]).real();
        return 0.5*(-J * SZ[i] * SZ[j] - 0.5*g*(SX[i] + SX[j])).real();
    }

    MatrixXd H() {
//        std::cout << h(0)+ h(1)<<endl << endl;
        return h(0) + h(1);
//        return h(0);
    }

    Tensor4d TimeEvolution_4th_order(const double delta) {
        double delta1 = delta /(4.0-pow(4.0,1.0/3.0));
        double delta2 = delta1;
        double delta3 = delta - 2.0*delta1 - 2.0*delta2;
        return Matrix_to_Tensor<4,double>(U(delta1)*U(delta2)*U(delta3)*U(delta2)*U(delta1), array4{2,2,2,2});

    }
    Tensor4d TimeEvolution_2nd_order(const double delta) {
//        cout <<setprecision(8) << fixed << (-delta*h(0)).exp() << endl;
//        return matrix_to_tensor<4>( (-delta*h(0)).exp().eval(), array4{2,2,2,2});
        return Matrix_to_Tensor<4,double>( ((-delta*h(0)).exp() * (-delta*h(1)).exp()).eval(), array4{2,2,2,2});
    }
    Tensor4d TimeEvolution_1st_order(const double delta) {
        return Matrix_to_Tensor<4,double>(U(delta), array4{2,2,2,2});
    }
    MatrixXd U(double delta){
        return (-delta*h(0)/2.0).exp() * (-delta*h(1)).exp() * (-delta * h(0)/2.0).exp();
    }

    Tensor4d W() {
        //Returns the H_two as A rank 4 MPO. Notation following Schollwöck (2010)
        MatrixXcd W(6,6);
        W.setZero();
        W.block(0,0,2,2) = I();
        W.block(2,0,2,2) = sz();
        W.block(4,0,2,2) = -g*sx();
        W.block(4,2,2,2) = -J*sz();
        W.block(4,4,2,2) = I();

        return Matrix_to_Tensor2<double>(W.real()).reshape(array4{2,3,2,3}).shuffle(array4{3,1,2,0});

    }
}