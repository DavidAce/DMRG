//
// Created by david on 7/17/17.
//

#ifndef DMRG_CLASS_HAMILTONIAN_H
#define DMRG_CLASS_HAMILTONIAN_H
#include <unsupported/Eigen/MatrixFunctions>
#include "mps_routines/class_mps.h"
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


class class_mpo{
public:
    class_mpo(){};
    using                       Scalar                = double;
    static constexpr int        mps_sites             = 2;                                                                              // Two site MPS
    long                        local_dimension       = Model::local_dimension;                                                         // "Spin" dimension
    MatrixType<Scalar>  H_asMatrix      = Model::H(mps_sites).real();                                                      // Matrix representation of full 2-site Hamiltonian
    Tensor<Scalar,4>    H_asTensor      = Matrix_to_Tensor<std::complex<double>,4>(Model::H(mps_sites), {2,2,2,2}).real(); // Rank-4 representation 2-site Hamiltonian (non MPO).
    MatrixType<Scalar>  H_MPO_asMatrix  = Model::MPO_asMatrix().real();                                                    // Matrix representation of full 2-site Hamiltonian
    Tensor<double,2>    H_MPO_asTensor2 = Matrix_to_Tensor2(Model::MPO_asMatrix()).real();                                                    // Matrix representation of full 2-site Hamiltonian
private:


    std::array<MatrixType<std::complex<double>>,2> h  = {Model::h(mps_sites,0), Model::h(mps_sites,1)};
    Tensor<Scalar,4>               compute_M();
    Tensor<Scalar,6>               compute_MM();

    template<typename T>
    MatrixType<std::complex<double>> Suzuki_Trotter_1st_order(T t){
        return (t*h[0]).exp() * (t*h[1]).exp();
    }

    template<typename T>
    MatrixType<std::complex<double>> Suzuki_Trotter_2nd_order(T t){
        return (t*h[0]/2.0).exp() * (t*h[1]).exp() * (t * h[0]/2.0).exp();
    }



    Tensor<std::complex<double>,4> TimeEvolution_1st_order(const double delta_t);
    Tensor<std::complex<double>,4> TimeEvolution_2nd_order(const double delta_t);
    Tensor<std::complex<double>,4> TimeEvolution_4th_order(const double delta_t);
    Tensor<Scalar,4>               compute_Udt(const double delta_t, const int order);

public:

    Tensor<std::complex<double>,4> compute_F(double a);
    Tensor<std::complex<double>,4> compute_G(double a);
    Tensor<std::complex<double>,4> compute_logG(double a);
    Tensor<Scalar,4>               Udt  = compute_Udt(0.01, 1);      // Rank-4 of 2-site unitary time evolution operator for iTEBD (non MPO).
    Tensor<Scalar,4>               M    = compute_M();               // MPO representation of 1-site Hamiltonian
    Tensor<Scalar,6>               MM   = compute_MM();              // MPO representation of 2-site Hamiltonian
    Tensor<std::complex<double>,4> F    = compute_F(0.0001);           // MPO representation of 1-site moment generating function
    Tensor<std::complex<double>,4> G    = compute_G(0.0001);           // MPO representation of 1-site characteristic function of the Hamiltonian
    void update_timestep(const double delta_t, const int order);
};


#endif //DMRG_CLASS_HAMILTONIAN_H
