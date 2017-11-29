//
// Created by david on 4/25/17.
//

#ifndef MODEL_H
#define MODEL_H


#include <vector>
#include <array>
#include "general/n_tensor_extra.h"
#include "general/n_math.h"

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

//    using Scalar = double; //Type for the groundstate wavefunction. Typically just double if the Hamiltonian is Real and Symmetric or Hermitian.

    inline double J = 1.0;
    inline double g = 1.0;
    inline long local_dimension = 2;

    //Transverse field Ising model:
    extern double get_exact_energy();
    inline double energy_exact  = get_exact_energy();  // = -1.063544409973372 if g = 0.5

    //Pauli matrices
    extern const Matrix2cd sx();
    extern const Matrix2cd sy();
    extern const Matrix2cd sz();
    extern const Matrix2cd I();
    //Pauli spin variables in N-dimensional manybody Hilbert space.
    extern std::vector<MatrixXcd>      gen_manybody_spin(Matrix2cd s);
    inline bool spins_must_be_generated = true;
    inline std::vector<MatrixXcd> SX;
    inline std::vector<MatrixXcd> SY;
    inline std::vector<MatrixXcd> SZ;
    extern MatrixXcd h(int sites, int position);
    extern MatrixXcd H(int sites);
    extern MatrixXcd MPO_asMatrix();



};



#endif //MPS_EIGEN_N_MODEL_H
