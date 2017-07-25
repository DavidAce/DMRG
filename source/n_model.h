//
// Created by david on 4/25/17.
//

#ifndef MODEL_H
#define MODEL_H


#include <vector>
#include <array>
#include <n_tensor_extra.h>
#include "funcs.h"

using namespace Textra;
using namespace Eigen;

namespace Model {

    //XXZ model:

    extern double J;
    extern double g;
    extern long local_dimension;

    //Pauli matrices
    extern Matrix2cd sx;
    extern Matrix2cd sy;
    extern Matrix2cd sz;
    extern Matrix2cd I;

    //Spin variables in L-dimensional hilbert space.
    extern std::vector<MatrixXcd> SX;
    extern std::vector<MatrixXcd> SY;
    extern std::vector<MatrixXcd> SZ;

    void generate_spins(const int sites);
    MatrixType Hamiltonian(const int sites);
    Tensor4 W(const int sites);

};


#endif //MPS_EIGEN_N_MODEL_H
