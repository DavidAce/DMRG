//
// Created by david on 2018-07-04.
//

#include "class_tf_ising.h"
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_random_numbers.h>
#include <iomanip>

using namespace qm::SpinOneHalf;
using Scalar = std::complex<double>;

void class_tf_ising::build_mpo()
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
    MPO.resize(3, 3, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(I);
    MPO.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(sx);
    MPO.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-(g_field + random_field) * sz);
    MPO.slice(Eigen::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-J_coupling * sx);
    MPO.slice(Eigen::array<long, 4>{2, 2, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(I);
}

void class_tf_ising::randomize_field(){
    random_field = rn::uniform_double(-randomness_strength,randomness_strength);
    MPO.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-(g_field+random_field) * sz);
}

Eigen::Tensor<Scalar,4> class_tf_ising::MPO_reduced_view() const {
    if (energy_reduced == 0){return MPO;}
    return MPO_reduced_view(energy_reduced);
}

Eigen::Tensor<Scalar,4> class_tf_ising::MPO_reduced_view(double single_site_energy) const {
    if (single_site_energy == 0){return MPO;}
    Eigen::Tensor<Scalar,4> temp  = MPO;
    temp.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-(g_field+random_field) * sz - single_site_energy * I);
    return temp;
}


std::unique_ptr<class_hamiltonian_base> class_tf_ising::clone() const {
    return std::make_unique<class_tf_ising>(*this);
}


void class_tf_ising::print_parameter_names() const {
    std::cout
            << std::setprecision(10)
            << std::setw(16) << std::left << "MPO"
            << std::setw(16) << std::left << "J"
            << std::setw(16) << std::left << "g"
            << std::setw(16) << std::left << "r"
            << std::setw(16) << std::left << "w"
            << std::setw(16) << std::left << "e"
            << std::setw(16) << std::left << "d"
            << std::endl;

}

void class_tf_ising::print_parameter_values() const {
    std::cout
            << std::setprecision(10)
            << std::setw(16) << std::left << get_position()
            << std::setw(16) << std::left << J_coupling
            << std::setw(16) << std::left << g_field
            << std::setw(16) << std::left << random_field
            << std::setw(16) << std::left << randomness_strength
            << std::setw(16) << std::left << energy_reduced
            << std::setw(16) << std::left << spin_dim
            << std::endl;

}


std::vector<double> class_tf_ising::get_all_parameters() const {
    return {(double)position,J_coupling,g_field,random_field,randomness_strength,energy_reduced,(double)spin_dim};
}

