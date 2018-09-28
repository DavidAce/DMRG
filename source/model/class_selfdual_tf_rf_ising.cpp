//
// Created by david on 2018-07-06.
//

#include "class_selfdual_tf_rf_ising.h"
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_random_numbers.h>
#include <iomanip>
#include <general/nmspc_math.h>

using namespace qm::spinOneHalf;
using Scalar = std::complex<double>;


class_selfdual_tf_rf_ising::class_selfdual_tf_rf_ising(): class_hamiltonian_base(){
    extent4     = {1, 1, spin_dim, spin_dim};
    extent2     = {spin_dim, spin_dim};
    J_rnd       = rn::log_normal(J_log_mean,J_sigma);
    h_rnd       = rn::log_normal(h_log_mean,h_sigma);
    build_mpo();
}



void class_selfdual_tf_rf_ising::build_mpo()
/*! Builds the MPO hamiltonian as a rank 4 tensor. Notation following Schollwöck (2010)
 *
 * H = - Σ J_{i} sz_{i} sz_{i+1} +  h_{i} sx_{i} + l*(h sx_i sx_{i+1} + J sz_{i} sz_{i+2})
 *
 *  |     I                 0           0              0            0   |
 *  |     sz                0           0              0            0   |
 *  |     sx                0           0              0            0   |
 *  |     0                 I           0              0            0   |
 *  | -(h_rnd)*sx       -J_rnd*sz   -l*h_mean*sx     -l*J_mean*sz     I   |
 *
 *        2
 *        |
 *    0---H---1
 *        |
 *        3
 *
 */
{
    MPO.resize(5, 5, spin_dim, spin_dim);
    MPO.setZero();
    MPO.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(I);
    MPO.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(sz);
    MPO.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(sx);
    MPO.slice(Eigen::array<long, 4>{3, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(I);
    MPO.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-h_rnd * sx);
    MPO.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-J_rnd * sz);
    MPO.slice(Eigen::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-(lambda*h_mean) * sx);
    MPO.slice(Eigen::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-(lambda*J_mean) * sz);
    MPO.slice(Eigen::array<long, 4>{4, 4, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(I);
}

void class_selfdual_tf_rf_ising::randomize_hamiltonian(){
    J_rnd       = rn::log_normal(J_log_mean,J_sigma);
    h_rnd       = rn::log_normal(h_log_mean,h_sigma);
    MPO.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-h_rnd * sx);
    MPO.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-J_rnd * sz);
}

Eigen::Tensor<Scalar,4> class_selfdual_tf_rf_ising::MPO_reduced_view() const {
    if (e_reduced == 0){return MPO;}
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar,4> class_selfdual_tf_rf_ising::MPO_reduced_view(double site_energy) const {
    if (site_energy == 0){return MPO;}
    Eigen::Tensor<Scalar,4> temp  = MPO;
    temp.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-h_rnd * sx - site_energy * I);
    return temp;
}

Eigen::MatrixXcd class_selfdual_tf_rf_ising::single_site_hamiltonian(
        int position,
        int sites,
        std::vector<Eigen::MatrixXcd> &SX,
        std::vector<Eigen::MatrixXcd> &SY[[maybe_unused]],
        std::vector<Eigen::MatrixXcd> &SZ)
        const
{
    int i = Math::mod(position,     sites);
    int j = Math::mod(position + 1, sites);
    int k = Math::mod(position + 2, sites);
    return - (J_rnd * SZ[i]*SZ[j] + h_rnd*0.5*(SX[i]+SX[j]) + lambda*(h_mean * SX[i]*SX[j] + J_mean*SZ[i]*SZ[k]));
}


std::unique_ptr<class_hamiltonian_base> class_selfdual_tf_rf_ising::clone() const {return std::make_unique<class_selfdual_tf_rf_ising>(*this);}


void   class_selfdual_tf_rf_ising::set_reduced_energy(double site_energy)         {e_reduced = site_energy;}
int    class_selfdual_tf_rf_ising::get_spin_dimension()                        const {return spin_dim;}
//double class_selfdual_tf_rf_ising::get_energy_reduced()                        const {return e_reduced;}
//double class_selfdual_tf_rf_ising::get_random_field()                          const {return h_rnd;}
//double class_selfdual_tf_rf_ising::get_randomness_strength()                   const {return w_rnd_strength;}


void class_selfdual_tf_rf_ising::print_parameter_names() const {
    std::cout
            << std::setprecision(10)
            << std::setw(16) << std::left << "MPO #"
            << std::setw(16) << std::left << "J_rnd"
            << std::setw(16) << std::left << "h_rnd"
            << std::setw(16) << std::left << "J_log_mean"
            << std::setw(16) << std::left << "h_log_mean"
            << std::setw(16) << std::left << "J_sigma"
            << std::setw(16) << std::left << "h_sigma"
            << std::setw(16) << std::left << "lambda"
            << std::setw(16) << std::left << "e_reduced"
            << std::setw(16) << std::left << "spin_dim"
            << std::endl;
}

void class_selfdual_tf_rf_ising::print_parameter_values() const {
    std::cout
            << std::setprecision(10)
            << std::setw(16) << std::left << get_position()
            << std::setw(16) << std::left << J_rnd
            << std::setw(16) << std::left << h_rnd
            << std::setw(16) << std::left << J_log_mean
            << std::setw(16) << std::left << h_log_mean
            << std::setw(16) << std::left << J_sigma
            << std::setw(16) << std::left << h_sigma
            << std::setw(16) << std::left << lambda
            << std::setw(16) << std::left << e_reduced
            << std::setw(16) << std::left << spin_dim
            << std::endl;
}

std::vector<std::string> class_selfdual_tf_rf_ising::get_parameter_names() const {
    return {"position",
            "J_rnd",
            "h_rnd",
            "J_log_mean",
            "h_log_mean",
            "J_sigma",
            "h_sigma",
            "lambda",
            "e_reduced",
            "spin_dim"
    };
}



std::vector<double> class_selfdual_tf_rf_ising::get_parameter_values() const {
    return {(double)get_position(),
            J_rnd,
            h_rnd,
            J_log_mean,
            h_log_mean,
            J_sigma,
            h_sigma,
            lambda,
            e_reduced,
            (double)spin_dim
            };
}

