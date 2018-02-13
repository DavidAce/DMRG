//
// Created by david on 7/17/17.
//

#ifndef DMRG_CLASS_HAMILTONIAN_H
#define DMRG_CLASS_HAMILTONIAN_H
#include <unsupported/Eigen/MatrixFunctions>
#include "mps_routines/class_mps.h"
#include "general/nmspc_tensor_extra.h"
#include "sim_parameters/nmspc_model.h"
#include "sim_parameters/nmspc_sim_settings.h"
//
//using namespace std;
//using namespace Textra;
//using namespace Eigen;


/*!
 * \class class_mpo
 * \brief MPO representations of the Hamiltonian.
 * Class for the Hamiltonian MPO, either as a rank-4 MPO used in DMRG, or as a time evolution operator \f$ \exp(-iH\delta t) \f$ as used in TEBD.
*/


class class_mpo{
public:
    class_mpo(){};
    using                       Scalar                = double;
    static constexpr int        mps_sites             = 2;                                                                              // Two site MPS
    long                        local_dimension       = Model::local_dimension;    // "Spin" dimension
    double                      timestep              ;
    Textra::MatrixType<Scalar>  H_asMatrix      = Model::H(mps_sites).real();                                                      // Matrix representation of full 2-site Hamiltonian
    Textra::Tensor<Scalar,4>    H_asTensor      = Textra::Matrix_to_Tensor<std::complex<double>,4>(Model::H(mps_sites), {2,2,2,2}).real(); // Rank-4 representation 2-site Hamiltonian (non MPO).
    Textra::MatrixType<Scalar>  H_MPO_asMatrix  = Model::MPO_asMatrix().real();                                                    // Matrix representation of full 2-site Hamiltonian
    Textra::Tensor<double,2>    H_MPO_asTensor2 = Textra::Matrix_to_Tensor2(Model::MPO_asMatrix()).real();                                                    // Matrix representation of full 2-site Hamiltonian
private:


    std::array<Textra::MatrixType<std::complex<double>>,2> h  = {Model::h(mps_sites,0), Model::h(mps_sites,1)};
    Textra::Tensor<Scalar,4>               compute_M();
    Textra::Tensor<Scalar,6>               compute_MM();

    template<typename T>
    Textra::MatrixType<std::complex<double>> Suzuki_Trotter_1st_order(T t){
        return (t*h[0]).exp() * (t*h[1]).exp();
    }

    template<typename T>
    Textra::MatrixType<std::complex<double>> Suzuki_Trotter_2nd_order(T t){
        return (t*h[0]/2.0).exp() * (t*h[1]).exp() * (t * h[0]/2.0).exp();
    }



    Textra::Tensor<std::complex<double>,4> TimeEvolution_1st_order(const double delta_t);
    Textra::Tensor<std::complex<double>,4> TimeEvolution_2nd_order(const double delta_t);
    Textra::Tensor<std::complex<double>,4> TimeEvolution_4th_order(const double delta_t);
    Textra::Tensor<Scalar,4>               compute_Udt(const double delta_t, const int order);

public:

    Textra::Tensor<std::complex<double>,4> compute_F(double a);
    Textra::Tensor<std::complex<double>,4> compute_G(double a);
    Textra::Tensor<std::complex<double>,4> compute_logG(double a);
    Textra::Tensor<Scalar,4>               Udt  = compute_Udt(0.01, 1);      // Rank-4 of 2-site unitary time evolution operator for iTEBD (non MPO).
    Textra::Tensor<Scalar,4>               M    = compute_M();               // MPO representation of 1-site Hamiltonian
    Textra::Tensor<Scalar,6>               MM   = compute_MM();              // MPO representation of 2-site Hamiltonian
    Textra::Tensor<std::complex<double>,4> F    = compute_F(0.0001);           // MPO representation of 1-site moment generating function
    Textra::Tensor<std::complex<double>,4> G    = compute_G(0.0001);           // MPO representation of 1-site characteristic function of the Hamiltonian
    void update_timestep(const double delta_t, const int susuki_trotter_order);
};


#endif //DMRG_CLASS_HAMILTONIAN_H
