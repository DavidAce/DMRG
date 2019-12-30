//
// Created by david on 2018-07-04.
//

#include "class_tf_ising.h"
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <math/nmspc_random.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <simulation/nmspc_settings.h>

using namespace qm::spinOneHalf;
using Scalar = std::complex<double>;


class_tf_ising::class_tf_ising(size_t position_, std::string logName): class_model_base(position_,logName){
    spin_dim            = settings::model::tf_ising::d;
    J_coupling          = settings::model::tf_ising::J;
    g_mag_field         = settings::model::tf_ising::g;
    w_rnd_strength      = settings::model::tf_ising::w;

    extent4     = {1, 1, spin_dim, spin_dim};
    extent2     = {spin_dim, spin_dim};
    r_rnd_field = rn::uniform_double_box(-w_rnd_strength, w_rnd_strength);

    qm::spinOneHalf::SX = qm::gen_manybody_spin(sx, 2);
    qm::spinOneHalf::SY = qm::gen_manybody_spin(sy, 2);
    qm::spinOneHalf::SZ = qm::gen_manybody_spin(sz, 2);
    qm::spinOneHalf::II = qm::gen_manybody_spin(Id , 2);

    all_mpo_parameters_have_been_set = true; //There are no full latice parameters on this model!
    build_mpo();
}


void   class_tf_ising::set_hamiltonian(const Eigen::Tensor<Scalar,4> & MPO_, std::vector<double> & parameters) {
    mpo_internal = MPO_;
    set_hamiltonian(parameters);
    auto mpo1 = Eigen::Map<const Eigen::VectorXcd>(MPO_ .data(), MPO_ .size());
    auto mpo2 = Eigen::Map<const Eigen::VectorXcd>(MPO().data(), MPO().size());
    assert(mpo1 == mpo2 and "MPO mismatch!");
    if (mpo1 != mpo2) throw std::runtime_error("MPO mismatch");
}

void   class_tf_ising::set_hamiltonian(const std::vector<double> & parameters) {
    auto temp = Eigen::Map<const Eigen::VectorXd>(parameters.data(),parameters.size());
    set_hamiltonian(temp);
}


void   class_tf_ising::set_hamiltonian(const Eigen::MatrixXd & all_parameters, int position) {
    set_hamiltonian (all_parameters.row(position));
}


void   class_tf_ising::set_hamiltonian(const Eigen::VectorXd & parameters) {
    if((int)parameters.size() != num_params ) throw std::runtime_error("Wrong number of parameters given to initialize this model");
    position        = parameters(0);
    J_coupling      = parameters(1);
    g_mag_field     = parameters(2);
    r_rnd_field     = parameters(3);
    w_rnd_strength  = parameters(4);
    e_reduced       = parameters(5);
    spin_dim        = (int)parameters(6);
    all_mpo_parameters_have_been_set = true;
    build_mpo();
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
    if (not all_mpo_parameters_have_been_set) throw std::runtime_error("Improperly built MPO: Full lattice parameters haven't been set yet.");
    mpo_internal.resize(3, 3, spin_dim, spin_dim);
    mpo_internal.setZero();
    mpo_internal.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(g_mag_field + r_ptb_field + r_rnd_field) * sx - e_reduced * Id);
    mpo_internal.slice(Eigen::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-J_coupling * sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 2, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
}

void class_tf_ising::randomize_hamiltonian(){
    r_rnd_field = rn::uniform_double_box(-w_rnd_strength, w_rnd_strength);
    if(all_mpo_parameters_have_been_set or mpo_internal.size()>3){
        mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(g_mag_field+r_ptb_field + r_rnd_field) * sx - e_reduced * Id);
    }
}

void class_tf_ising::perturb_hamiltonian(double amplitude){
    r_ptb_field  = amplitude * r_rnd_field * rn::uniform_double_box(-1, 1);
    if(all_mpo_parameters_have_been_set or mpo_internal.size()>3){
        mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(g_mag_field+r_ptb_field + r_rnd_field) * sx - e_reduced * Id);
    }
    if (amplitude == 0.0 and is_perturbed())
        throw std::runtime_error("MPO(" + std::to_string(get_position()) + ": Should have unperturbed!");
}

bool class_tf_ising::is_perturbed() const {
    return r_ptb_field != 0.0;
}

Eigen::Tensor<Scalar,4> class_tf_ising::MPO_reduced_view() const {
    if (e_reduced == 0){return MPO();}
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar,4> class_tf_ising::MPO_reduced_view(double site_energy) const {
    if (site_energy == 0){return MPO();}
    Eigen::Tensor<Scalar,4> temp  = MPO();
    temp.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(g_mag_field + (r_rnd_field + r_ptb_field)) * sx - site_energy * Id);
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
    int i = math::mod(position,     sites);
    int j = math::mod(position + 1, sites);
    return -(J_coupling * SZ[i] * SZ[j] + g_mag_field * 0.5*(SX[i]+SX[j])) ;
}


std::unique_ptr<class_model_base> class_tf_ising::clone() const {return std::make_unique<class_tf_ising>(*this);}

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

//
//std::vector<double> class_tf_ising::get_random_parameter_values() const {
//    return {(double)get_position(),
//            J_rnd,
//            h_rnd,
//            J_log_mean,
//            h_log_mean,
//            J_avg,
//            h_avg,
//            J_sigma,
//            h_sigma,
//            lambda,
//            delta,
//            e_reduced,
//            (double)spin_dim
//    };
//}

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


void class_tf_ising::set_full_lattice_parameters([[maybe_unused]] const std::vector<std::vector<double>> chain_parameters, bool reverse){
    all_mpo_parameters_have_been_set = true;
    // Calculate average J_rnd on the whole state
    all_mpo_parameters_have_been_set = true;
    std::list<double> r_rnd_list;
    if(reverse){
        for (auto &params : chain_parameters){
            r_rnd_list.emplace_front(params[3]);
        }
    }else{
        for (auto &params : chain_parameters){
            r_rnd_list.push_back(params[3]);
        }
    }

    r_rnd_field  = *std::next(r_rnd_list.begin(), get_position());

}
