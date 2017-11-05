//
// Created by david on 7/17/17.
//

#ifndef DMRG_CLASS_HAMILTONIAN_H
#define DMRG_CLASS_HAMILTONIAN_H

#include "general/n_tensor_extra.h"

using namespace std;
using namespace Textra;
using namespace Eigen;
class class_TwoSiteHamiltonian{
public:
    class_TwoSiteHamiltonian();

    unsigned int sites = 2;     // Two sites
    long local_dimension;       // "Spin" dimension
    string picture;             // Graphical representation
    Tensor4d W;                  // MPO representation of H_two
    MatrixType<double> asMatrix2;       // Matrix  representation of full two-site H_two
    MatrixType<double> asMatrix;        // Matrix  representation of full two-site H_two
    Tensor4d asTensor4;          // Tensor4d representation of full two-site H_two.
    Tensor4d asTimeEvolution;    // Tensor4d unitary time evolution operator for iTEBD.
//    Tensor4d asTimeEvolution2;    // Tensor4d unitary time evolution operator for iTEBD.
};


#endif //DMRG_CLASS_HAMILTONIAN_H
