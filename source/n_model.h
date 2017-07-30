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

/*! \brief Setup the Hamiltonian and corresponding values.*/

/*!
 *  \namespace Model
 *  This namespace contains model parameters, such as the coupling (or nearest neighbor exchange energy)  \f$J\f$, and
 *  the transverse field strength \f$g\f$. By default the model describes the quantum Ising model with a transverse field:
 *
 *  \f[
 *  H = \frac{1}{2}\sum_{n=1}^L J S^z_n \otimes S^z_{n+1} + gS^x_{n}
 *  \f]
 *
 */

namespace Model {

    //Transverse field Ising model:

    extern double J;
    extern double g;
    extern long local_dimension;
    extern double energy_exact;
    //Pauli matrices
    extern Matrix2cd sx;
    extern Matrix2cd sy;
    extern Matrix2cd sz;
    extern Matrix2cd I;

    //Spin variables in L-dimensional Hilbert space. L being the particle chain length.
    extern std::vector<MatrixXcd> SX;
    extern std::vector<MatrixXcd> SY;
    extern std::vector<MatrixXcd> SZ;

    void generate_spins(const int sites);
    MatrixType Hamiltonian(const int sites);
    Tensor4 W(const int sites);

};


#endif //MPS_EIGEN_N_MODEL_H
