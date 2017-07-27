//
// Created by david on 7/17/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_HAMILTONIAN_H
#define FINITE_DMRG_EIGEN_CLASS_HAMILTONIAN_H

#include <n_tensor_extra.h>
#include <class_environment.h>
#include <SymEigsSolver.h>

using namespace std;
using namespace Textra;

class class_TwoSiteHamiltonian{
public:
    class_TwoSiteHamiltonian();

    unsigned int sites = 2;     // Two sites
    long local_dimension;       // "Spin" dimension
    string picture;             // Graphical representation
    Tensor4 W;                  // MPO representation of Hamiltonian
    MatrixType asMatrix;        // Matrix   representation of full two-site Hamiltonian
    Tensor4 asTensor4;          // Tensor4 representation of full two-site Hamiltonian.

    double delta = 0.005;
    Tensor4 asTimeEvolution;    // Tensor4 unitary time evolution operator for iTEBD.
};


#endif //FINITE_DMRG_EIGEN_CLASS_HAMILTONIAN_H
