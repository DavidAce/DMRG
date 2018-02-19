//
// Created by david on 7/17/17.
//

#ifndef DMRG_CLASS_HAMILTONIAN_H
#define DMRG_CLASS_HAMILTONIAN_H


#include <unsupported/Eigen/MatrixFunctions>
#include "mps_routines/class_mps.h"
#include <general/nmspc_tensor_extra.h>
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

template <typename Scalar>
class class_mpo{
public:
    class_mpo();
    static constexpr int        mps_sites = 2  ;    /*!< Two site MPS */
    long                        local_dimension;    /*!< Local "spin" dimension */
    double                      timestep;
    Textra::MatrixType<Scalar>  H_asMatrix;         /*!< Matrix representation of full 2-site Hamiltonian */
    Textra::Tensor<Scalar,4>    H_asTensor;         /*!< Rank-4 representation 2-site Hamiltonian (non MPO). */
    Textra::MatrixType<Scalar>  H_MPO_asMatrix;     /*!< Matrix representation of full 2-site Hamiltonian */
    Textra::Tensor<Scalar,2>    H_MPO_asTensor2;    /*!< Matrix representation of full 2-site Hamiltonian */

    Textra::Tensor<std::complex<double>,4> compute_F(double a);
    Textra::Tensor<std::complex<double>,4> compute_G(double a);
    Textra::Tensor<std::complex<double>,4> compute_logG(double a);
    Textra::Tensor<Scalar,4>               Udt;                 /*!< Rank-4 of 2-site unitary time evolution operator for iTEBD (non MPO). */
    Textra::Tensor<Scalar,4>               M  ;                 /*!< MPO representation of 1-site Hamiltonian */
    Textra::Tensor<Scalar,6>               MM ;                 /*!< MPO representation of 2-site Hamiltonian */
    Textra::Tensor<Scalar,8>               MMMM ;               /*!< MPO representation of 2-site Hamiltonian */
    Textra::Tensor<std::complex<double>,4> F  ;                 /*!< MPO representation of 1-site moment generating function */
    Textra::Tensor<std::complex<double>,4> G  ;                 /*!< MPO representation of 1-site characteristic function of the Hamiltonian */
    void update_timestep(const double delta_t, const int susuki_trotter_order);

private:


    std::array<Textra::MatrixType<std::complex<double>>,2> h;  /*!< Local even and odd hamiltonians */

    Textra::Tensor<Scalar,4>  compute_M();      /*!< Returns an MPO representation of 1-site Hamiltonian */
    Textra::Tensor<Scalar,6>  compute_MM();     /*!< Returns an MPO representation of 2-site Hamiltonian */
    Textra::Tensor<Scalar,8>  compute_MMMM();   /*!< Returns an MPO representation of 2-site double Hamiltonian */

    template<typename T>
    Textra::MatrixType<std::complex<double>> Suzuki_Trotter_1st_order(T t){
        return (t*h[0]).exp() * (t*h[1]).exp();
    }

    template<typename T>
    Textra::MatrixType<std::complex<double>> Suzuki_Trotter_2nd_order(T t){
        return (t*h[0]/2.0).exp() * (t*h[1]).exp() * (t * h[0]/2.0).exp();
    }


    template<typename T>
    Textra::MatrixType<std::complex<double>> Suzuki_Trotter_4th_order(T t){
        T t1 = t /(4.0-pow(4.0,1.0/3.0));
        T t2 = t1;
        T t3 = t - 2.0*t1 - 2.0*t2;
        return Suzuki_Trotter_2nd_order(t1)
               *Suzuki_Trotter_2nd_order(t2)
               *Suzuki_Trotter_2nd_order(t3)
               *Suzuki_Trotter_2nd_order(t2)
               *Suzuki_Trotter_2nd_order(t1);
    }

    Textra::Tensor<std::complex<double>,4> TimeEvolution_1st_order(const double delta_t);
    Textra::Tensor<std::complex<double>,4> TimeEvolution_2nd_order(const double delta_t);
    Textra::Tensor<std::complex<double>,4> TimeEvolution_4th_order(const double delta_t);
    Textra::Tensor<Scalar,4>               compute_Udt(const double delta_t, const int order);

};


#endif //DMRG_CLASS_HAMILTONIAN_H
