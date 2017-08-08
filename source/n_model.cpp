//
// Created by david on 4/25/17.
//

#include "n_model.h"
#include <unsupported/Eigen/KroneckerProduct>

using namespace std;

namespace Model{
    double J = 1.0;
    double g = 0.5;
    double energy_exact = -1.063544409973372; //Exact ground state energy.

    long local_dimension = 2;
    Matrix2cd sx;
    Matrix2cd sy;
    Matrix2cd sz;
    Matrix2cd I;
    std::vector<MatrixXcd> SY;
    std::vector<MatrixXcd> SZ;
    std::vector<MatrixXcd> SX;

    void generate_spins(const int sites) {
        sx << 0.0, 1.0,
              1.0, 0.0;
        sy << 0.0, -1.0i,
              1.0i, 0.0;
        sz << 1.0, 0.0,
              0.0, -1.0;

        I.setIdentity();
        MatrixXcd X, Y, Z;
        SX.clear();
        SY.clear();
        SZ.clear();

        for (int i = 0; i < sites; i++) {
            X = i == 0 ? sx : I;
            Y = i == 0 ? sy : I;
            Z = i == 0 ? sz : I;
            for (int j = 1; j < sites; j++) {
                X = kroneckerProduct(X, i == j ? sx : I).eval();
                Y = kroneckerProduct(Y, i == j ? sy : I).eval();
                Z = kroneckerProduct(Z, i == j ? sz : I).eval();
            }
            SX.emplace_back(X);
            SY.emplace_back(Y);
            SZ.emplace_back(Z);
        }
    }

    void get_exact_energy(){


//        f = lambda k, g: -2 * np.sqrt(1 + g ** 2 - 2 * g * np.cos(k)) / np.pi / 2.
//        E0_exact = integrate.quad(f, 0, np.pi, args=(g,))[0]
//        print("E_exact =", E0_exact)
    }


    MatrixType Hamiltonian(const int sites) {
        generate_spins(sites);
        MatrixType H = MatrixType::Zero(ipow(2, sites), ipow(2, sites));

        for (int i = 0; i < sites; i++) {
            H += -0.5 * (J * SZ[i] * SZ[mod(i + 1, sites)] + g * SX[i]).real();
        }
        return H;
    }


    Tensor4 W(const int sites) {
        //Returns the Hamiltonian as A rank 4 MPO. Notation following Schollwöck (2010)
        generate_spins(sites);
        MatrixXcd W(6,6);
        W.block(0,0,2,2) = I;
        W.block(2,0,2,2) = sz;
        W.block(4,0,2,2) = -g*sx;
        W.block(4,2,2,2) = -J*sz;
        W.block(4,4,2,2) = I;
        return matrix_to_tensor2(W.real()).reshape(array4{2,3,2,3}).shuffle(array4{3,1,2,0});

    }
}