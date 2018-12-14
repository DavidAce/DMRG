//
// Created by david on 2018-07-04.
//

#include "class_tf_ising.h"
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_random_numbers.h>
#include <iomanip>
#include <general/nmspc_math.h>

using namespace qm::spinOneHalf;
using Scalar = std::complex<double>;


class_tf_ising::class_tf_ising(): class_hamiltonian_base(){
    extent4     = {1, 1, spin_dim, spin_dim};
    extent2     = {spin_dim, spin_dim};
    r_rnd_field = rn::uniform_double(-w_rnd_strength,w_rnd_strength);
    build_mpo();
    qm::spinOneHalf::SX = qm::gen_manybody_spin(sx, 2);
    qm::spinOneHalf::SY = qm::gen_manybody_spin(sy, 2);
    qm::spinOneHalf::SZ = qm::gen_manybody_spin(sz, 2);
    qm::spinOneHalf::II = qm::gen_manybody_spin(I , 2);
}



void class_tf_ising::build_mpo()
/*! Builds the MPO hamiltonian as a rank 4 tensor. Notation following Schollwöck (2010)

 * H = - Σ J sz_{i} sz_{i+1} +  g_{i} sx_{i}
 *
 *  |      I        0   0   |
 *  |     sz        0   0   |
 *  | -(g+r)*sx     sz  I   |
 *
 *        2
 *        |
 *    0---H---1
 *        |
 *        3
 *
 */
{
    MPO.resize(3, 3, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(I);
    MPO.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(sz);
    MPO.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-(g_mag_field + r_rnd_field) * sx);
    MPO.slice(Eigen::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-J_coupling * sz);
    MPO.slice(Eigen::array<long, 4>{2, 2, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(I);
}

void class_tf_ising::randomize_hamiltonian(){
    r_rnd_field = rn::uniform_double(-w_rnd_strength,w_rnd_strength);
    MPO.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-(g_mag_field+r_rnd_field) * sx);
}

Eigen::Tensor<Scalar,4> class_tf_ising::MPO_reduced_view() const {
    if (e_reduced == 0){return MPO;}
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar,4> class_tf_ising::MPO_reduced_view(double site_energy) const {
    if (site_energy == 0){return MPO;}
    Eigen::Tensor<Scalar,4> temp  = MPO;
    temp.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-(g_mag_field+r_rnd_field) * sx - site_energy * I);
    return temp;
}

Eigen::MatrixXcd class_tf_ising::single_site_hamiltonian(
        int position,
        int sites,
        std::vector<Eigen::MatrixXcd> &SX,
        std::vector<Eigen::MatrixXcd> &SY[[maybe_unused]],
        std::vector<Eigen::MatrixXcd> &SZ)
        const
{
    int i = Math::mod(position,     sites);
    int j = Math::mod(position + 1, sites);
    return -(J_coupling * SZ[i] * SZ[j] + g_mag_field * 0.5*(SX[i]+SX[j])) ;
}


std::unique_ptr<class_hamiltonian_base> class_tf_ising::clone() const {return std::make_unique<class_tf_ising>(*this);}
void   class_tf_ising::set_reduced_energy(double site_energy)             {e_reduced = site_energy;}
size_t class_tf_ising::get_spin_dimension()                         const {return spin_dim;}
//double class_tf_ising::get_energy_reduced()                         const {return e_reduced;}
//double class_tf_ising::get_random_field()                           const {return r_rnd_field;}
//double class_tf_ising::get_randomness_strength()                    const {return w_rnd_strength;}
//

void class_tf_ising::print_parameter_names() const {
    std::cout
            << std::setprecision(10)
            << std::setw(16) << std::left << "MPO"
            << std::setw(16) << std::left << "J"
            << std::setw(16) << std::left << "g"
            << std::setw(16) << std::left << "r"
            << std::setw(16) << std::left << "lambda"
            << std::setw(16) << std::left << "e"
            << std::setw(16) << std::left << "d"
            << std::endl;
}

void class_tf_ising::print_parameter_values() const {
    std::cout
            << std::setprecision(10)
            << std::setw(16) << std::left << get_position()
            << std::setw(16) << std::left << J_coupling
            << std::setw(16) << std::left << g_mag_field
            << std::setw(16) << std::left << r_rnd_field
            << std::setw(16) << std::left << w_rnd_strength
            << std::setw(16) << std::left << e_reduced
            << std::setw(16) << std::left << spin_dim
            << std::endl;
}


std::vector<std::string> class_tf_ising::get_parameter_names() const {
    return {"position",
            "J_coupling",
            "g_mag_field",
            "r_rnd_field",
            "w_rnd_strength",
            "e_reduced",
            "spin_dim"
            };
}



std::vector<double> class_tf_ising::get_parameter_values() const {
    return {(double)get_position(),
            J_coupling,
            g_mag_field,
            r_rnd_field,
            w_rnd_strength,
            e_reduced,
            (double)spin_dim
            };
}


void class_tf_ising::set_non_local_parameters(std::vector<std::vector<double>> &chain_parameters){

}
