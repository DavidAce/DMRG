//
// Created by david on 7/17/17.
//

#ifndef DMRG_CLASS_HAMILTONIAN_H
#define DMRG_CLASS_HAMILTONIAN_H


#include <general/nmspc_tensor_extra.h>

using namespace std::complex_literals;


/*!
 * \class class_mpo
 * \brief MPO representations of the Hamiltonian.
 * Class for the Hamiltonian MPO, either as a rank-4 MPO used in DMRG, or as a time evolution operator \f$ \exp(-iH\delta t) \f$ as used in TEBD.
*/

class class_mpo{

public:
    using Scalar = std::complex<double>;
    class_mpo();
    static constexpr int        mps_sites = 2  ;   /*!< Two site MPS */
    long                        local_dimension;   /*!< Local "spin" dimension */
    double                      step_size;
    Textra::MatrixType<Scalar>  H_asMatrix;        /*!< Matrix representation of full 2-site Hamiltonian */
    Textra::Tensor<Scalar,4>    H_asTensor;        /*!< Rank-4 representation 2-site Hamiltonian (non MPO). */
    Textra::Tensor<Scalar,4>    H_asTensor_sq;     /*!< Rank-4 representation 2-site Hamiltonian squared (non MPO). */
    Textra::MatrixType<Scalar>  H_MPO_asMatrix;    /*!< Matrix representation of full 2-site Hamiltonian */
    Textra::Tensor<Scalar,2>    H_MPO_asTensor2;   /*!< Matrix representation of full 2-site Hamiltonian */
    std::array<Textra::MatrixType<Scalar>,2> h;    /*!< Local even and odd hamiltonians */

    std::vector<Textra::Tensor<Scalar,4>>  U  ;    /*!< Set of rank-4 of 2-site unitary evolution operators (non-MPO). */
    Textra::Tensor<Scalar,4>               M  ;    /*!< MPO representation of 1-site Hamiltonian */
//    Textra::Tensor<Scalar,6>               MM ;    /*!< MPO representation of 2-site Hamiltonian */
//    Textra::Tensor<Scalar,8>               MMMM ;  /*!< MPO representation of 2-site Hamiltonian */
    std::vector<Textra::Tensor<Scalar,4>>  G0;     /*!< MPO representation of 1-site moment generating function */
    std::vector<Textra::Tensor<Scalar,4>>  G1;     /*!< MPO representation of 1-site characteristic function of the Hamiltonian */

    std::vector<Textra::Tensor<Scalar,4>> compute_G(Scalar a, int susuki_trotter_order);
    std::vector<Textra::Tensor<Scalar,4>> get_2site_evolution_gates(const Scalar t,int susuki_trotter_order);
    void update_evolution_step_size(const Scalar dt, const int susuki_trotter_order);

    Textra::Tensor<Scalar,4>  compute_M   (Scalar k = 0.0);         /*!< Returns an MPO representation of 1-site Hamiltonian */
//    Textra::Tensor<Scalar,6>  compute_MM  (Scalar k = 0.0);        /*!< Returns an MPO representation of 2-site Hamiltonian */
//    Textra::Tensor<Scalar,8>  compute_MMMM(Scalar k = 0.0);      /*!< Returns an MPO representation of 2-site double Hamiltonian */
private:

    std::vector<Textra::MatrixType<Scalar>> Suzuki_Trotter_1st_order(Scalar t);
    std::vector<Textra::MatrixType<Scalar>> Suzuki_Trotter_2nd_order(Scalar t);
    std::vector<Textra::MatrixType<Scalar>> Suzuki_Trotter_4th_order(Scalar t);


};


#endif //DMRG_CLASS_HAMILTONIAN_H
