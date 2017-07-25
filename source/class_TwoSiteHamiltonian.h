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
    unsigned int sites = 2;     // Two sites
    long local_dimension;       // "Spin" dimension
    string pic;                 // Graphical representation
    Tensor4 W;                  // MPO representation of Hamiltonian
    MatrixType asMatrix;        // Matrix   representation of full two-site Hamiltonian
    Tensor4 asTensor4;          // Tenesor4 representation of full two-site Hamiltonian.

    class_TwoSiteHamiltonian();
};


#endif //FINITE_DMRG_EIGEN_CLASS_HAMILTONIAN_H
