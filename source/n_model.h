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
    extern std::vector<MatrixXcd> generate_manybody_spin(const int sites, Matrix2cd &s);

    extern MatrixType Hamiltonian(const int sites);
    extern MatrixType Hamiltonian2(const int sites);
    extern Tensor4 TimeEvolution_4th_order(const int sites, const double delta);
    extern MatrixType U(const double delta);
    extern Tensor4 W(const int sites);

    extern Matrix2cd gen_sx();
    extern Matrix2cd gen_sy();
    extern Matrix2cd gen_sz();

    extern std::vector<MatrixXcd> gen_SX(const int sites);
    extern std::vector<MatrixXcd> gen_SY(const int sites);
    extern std::vector<MatrixXcd> gen_SZ(const int sites);

    inline MatrixType H1;
    inline MatrixType H2;

    inline double J = -1.0;
    inline double g =  1.0;
    inline long local_dimension = 2;
    inline double energy_exact  = get_exact_energy();  // = -1.063544409973372 if g = 0.5

    //Pauli matrices
    inline Matrix2cd sx = gen_sx();
    inline Matrix2cd sy = gen_sy();
    inline Matrix2cd sz = gen_sz();
    inline Matrix2cd I = Matrix2cd::Identity();

    //Spin variables in L-dimensional Hilbert space. L being the particle chain length.
    inline std::vector<MatrixXcd> SX = generate_manybody_spin(2, sx);
    inline std::vector<MatrixXcd> SY = generate_manybody_spin(2, sy);
    inline std::vector<MatrixXcd> SZ = generate_manybody_spin(2, sz);



};


#endif //MPS_EIGEN_N_MODEL_H
