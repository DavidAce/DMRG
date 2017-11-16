//
// Created by david on 7/17/17.
//

#ifndef DMRG_CLASS_HAMILTONIAN_H
#define DMRG_CLASS_HAMILTONIAN_H

#include "general/n_tensor_extra.h"
#include "sim_parameters/n_model.h"
#include "sim_parameters/n_sim_settings.h"
using namespace std;
using namespace Textra;
using namespace Eigen;
/*!
 # Superblock Class
   This class contains the Hamiltonian MPO, current wave function MPS, left and right environment blocks and routines to contract, diagonalize, truncate
   and update them.
*/


class class_hamiltonian{
public:
    class_hamiltonian(){};
    using                       Scalar          = Model::Scalar;
    static constexpr int        mps_sites       = 2;                        // Two site MPS
    long                        local_dimension = Model::local_dimension;   // "Spin" dimension
    Textra::MatrixType<Scalar>  asMatrix        = Model::H(mps_sites);               // Matrix representation of full 2-site Hamiltonian
    Textra::Tensor<4,Scalar>    M               = Model::M();               // MPO representation of 1-site Hamiltonian
    Textra::Tensor<6,Scalar>    MM              = Model::MM();              // MPO representation of 2-site Hamiltonian
    Textra::Tensor<2*mps_sites,Scalar>asTensor  = Model::H_tensor<mps_sites>();           // Rank-4 representation 2-site Hamiltonian (non MPO).
    Textra::Tensor<4,Scalar>    Udt             = Model::TimeEvolution_1st_order(0.01, mps_sites);               // Rank-4 of 2-site unitary time evolution operator for iTEBD (non MPO).

    void update_timestep(double delta_t){
        Udt = Model::TimeEvolution_1st_order(delta_t, mps_sites);
    }
};


#endif //DMRG_CLASS_HAMILTONIAN_H
