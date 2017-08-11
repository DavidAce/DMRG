//
// Created by david on 4/25/17.
//

#ifndef MODEL_H
#define MODEL_H


#include <vector>
#include <array>
#include <n_tensor_extra.h>
#include "n_math.h"

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
    extern double get_exact_energy();
    extern void generate_spins(const int sites);
    extern MatrixType Hamiltonian(const int sites);
    extern Tensor4 W(const int sites);


    inline double J = -1.0;
    inline double g = 0.5;
    inline long local_dimension = 2;
    inline double energy_exact  = get_exact_energy();  // = -1.063544409973372 if g = 0.5
    //Pauli matrices
    inline Matrix2cd sx;
    inline Matrix2cd sy;
    inline Matrix2cd sz;
    inline Matrix2cd I;

    //Spin variables in L-dimensional Hilbert space. L being the particle chain length.
    inline std::vector<MatrixXcd> SX;
    inline std::vector<MatrixXcd> SY;
    inline std::vector<MatrixXcd> SZ;



};


#endif //MPS_EIGEN_N_MODEL_H
