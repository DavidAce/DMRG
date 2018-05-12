//
// Created by david on 2018-04-17.
//

#include "class_hamiltonian.h"
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
using namespace qm;
using Scalar = std::complex<double>;

void class_hamiltonian::build_mpo()
/*! Builds the MPO hamiltonian as a rank 4 tensor. Notation following Schollwöck (2010)
 *
 *          2
 *          |
 *      0---H---1
 *          |
 *          3
 */
{
    /*! Returns the MPO as a tensor. Notation following Schollwöck (2010) */
    MPO.resize(3, 3, 2, 2);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(I);
    MPO.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(sx);
    MPO.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-g * sz - e * I);
    MPO.slice(Eigen::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-J * sx);
    MPO.slice(Eigen::array<long, 4>{2, 2, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(I);
}

Eigen::Tensor<Scalar,4> class_hamiltonian::MPO_reduced() {
    if (e == 0){return MPO;}
    Eigen::Tensor<Scalar,4> temp  = MPO;
    temp.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-g * sz - e * I);
    return temp;
}

Eigen::Tensor<Scalar,4> class_hamiltonian::MPO_reduced(double energy_e_temporary) {
    if (energy_e_temporary == 0){return MPO;}
    Eigen::Tensor<Scalar,4> temp  = MPO;
    temp.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-g * sz - energy_e_temporary * I);
    return temp;
}

void class_hamiltonian::set_parameters(double coupling_J, double field_g, double energy_e) {
    J = coupling_J;
    e = energy_e;
    g = field_g;
    build_mpo();
}

void class_hamiltonian::update_site_coupling(double coupling_J) {
    J = coupling_J;
    MPO        .slice(Eigen::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-J * sx);

}
void class_hamiltonian::update_site_energy(double energy_e) {
    e = energy_e;
}

void class_hamiltonian::update_site_field(double field_g) {
    g = field_g;
    MPO        .slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-g * sz);
}


double class_hamiltonian::get_site_coupling()const {return J;}
double class_hamiltonian::get_site_energy()const {return e;}
double class_hamiltonian::get_site_field()const {return g;}


