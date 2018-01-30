//
// Created by david on 4/25/17.
//

#ifndef MODEL_H
#define MODEL_H


#include <vector>
#include <array>
#include <Eigen/Core>


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
    extern const Eigen::Matrix2cd sx();
    extern const Eigen::Matrix2cd sy();
    extern const Eigen::Matrix2cd sz();
    extern const Eigen::Matrix2cd I();
    //Pauli spin variables in N-dimensional manybody Hilbert space.
    extern std::vector<Eigen::MatrixXcd>      gen_manybody_spin(Eigen::Matrix2cd s);
    inline bool spins_must_be_generated = true;
    inline std::vector<Eigen::MatrixXcd> SX;
    inline std::vector<Eigen::MatrixXcd> SY;
    inline std::vector<Eigen::MatrixXcd> SZ;
    extern Eigen::MatrixXcd h(int sites, int position);
    extern Eigen::MatrixXcd H(int sites);
    extern Eigen::MatrixXcd MPO_asMatrix();



};



#endif //MPS_EIGEN_N_MODEL_H
