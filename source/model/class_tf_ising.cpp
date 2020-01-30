//
// Created by david on 2018-07-04.
//

#include "class_tf_ising.h"
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <simulation/nmspc_settings.h>

using namespace qm::spinOneHalf;
using Scalar = std::complex<double>;

class_tf_ising::class_tf_ising(size_t position_) : class_model_base(position_) {
    log      = Logger::setLogger("tf-ising");
    spin_dim = settings::model::tf_ising::d;
    J_nn     = settings::model::tf_ising::J;
    h_field  = settings::model::tf_ising::g;
    h_sigma  = settings::model::tf_ising::w;

    extent4 = {1, 1, spin_dim, spin_dim};
    extent2 = {spin_dim, spin_dim};
    h_rnd   = rn::uniform_double_box(-h_sigma, h_sigma);

    qm::spinOneHalf::SX = qm::gen_manybody_spin(sx, 2);
    qm::spinOneHalf::SY = qm::gen_manybody_spin(sy, 2);
    qm::spinOneHalf::SZ = qm::gen_manybody_spin(sz, 2);
    qm::spinOneHalf::II = qm::gen_manybody_spin(Id, 2);

    all_mpo_parameters_have_been_set = true; // There are no full lattice parameters on this model so we set it true immediately!
    class_tf_ising::build_mpo();
}

double class_tf_ising::get_field() const { return h_field + std::pow(h_ptb + h_rnd, 1 - beta); }
double class_tf_ising::get_coupling() const { return std::pow(J_nn, 1 - alpha); }

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
    if(not all_mpo_parameters_have_been_set) throw std::runtime_error("Improperly built MPO: Full lattice parameters haven't been set yet.");
    mpo_internal.resize(3, 3, spin_dim, spin_dim);
    mpo_internal.setZero();
    mpo_internal.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    mpo_internal.slice(Eigen::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 2, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
}

void class_tf_ising::randomize_hamiltonian() {
    h_rnd = rn::uniform_double_box(-h_sigma, h_sigma);
    if(all_mpo_parameters_have_been_set or mpo_internal.size() > 3) {
        mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    }
}

void class_tf_ising::set_coupling_damping(double alpha_) {
    alpha                            = alpha_;
}
void class_tf_ising::set_field_damping(double beta_) {
    beta                             = beta_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    }
}

void class_tf_ising::set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) {
    switch(ptbMode) {
        case PerturbMode::ABSOLUTE: {
            h_ptb = field_ptb;
            break;
        }
        case PerturbMode::PERCENTAGE: {
            h_ptb = h_rnd * field_ptb;
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_ABSOLUTE: {
            h_ptb = rn::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_PERCENTAGE: {
            h_ptb = h_rnd * rn::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
    }
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    }
    if(coupling_ptb == 0.0 and field_ptb == 0 and is_perturbed())
        throw std::runtime_error("MPO(" + std::to_string(get_position()) + ": Should have become unperturbed!");
}

bool class_tf_ising::is_perturbed() const { return h_ptb != 0.0; }
bool class_tf_ising::is_damped() const { return alpha != 0.0 or beta != 0.0; }

Eigen::Tensor<Scalar, 4> class_tf_ising::MPO_reduced_view() const {
    if(e_reduced == 0) {
        return MPO();
    }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> class_tf_ising::MPO_reduced_view(double site_energy) const {
    if(site_energy == 0) {
        return MPO();
    }
    Eigen::Tensor<Scalar, 4> temp                                           = MPO();
    temp.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - site_energy * Id);
    return temp;
}

Eigen::MatrixXcd class_tf_ising::single_site_hamiltonian(int position, int sites, std::vector<Eigen::MatrixXcd> &SX,
                                                         std::vector<Eigen::MatrixXcd> &SY [[maybe_unused]], std::vector<Eigen::MatrixXcd> &SZ) const {
    int i = math::mod(position, sites);
    int j = math::mod(position + 1, sites);
    return -(J_nn * SZ[i] * SZ[j] + h_field * 0.5 * (SX[i] + SX[j]));
}

std::unique_ptr<class_model_base> class_tf_ising::clone() const { return std::make_unique<class_tf_ising>(*this); }

size_t class_tf_ising::get_spin_dimension() const { return spin_dim; }

Eigen::Tensor<Scalar, 1> class_tf_ising::get_MPO_edge_left() const {
    Eigen::Tensor<Scalar, 1> ledge(3);
    ledge.setZero();
    ledge(2) = 1;
    return ledge;
}
Eigen::Tensor<Scalar, 1> class_tf_ising::get_MPO_edge_right() const {
    Eigen::Tensor<Scalar, 1> redge(3);
    redge.setZero();
    redge(0) = 1;
    return redge;
}

void class_tf_ising::set_parameters(const Parameters &parameters) {
    spin_dim     = get_val<size_t>(parameters, "spin_dim");
    position     = get_val<size_t>(parameters, "position");
    distribution = get_val<std::string>(parameters, "distribution");
    parity_sep   = get_val<bool>(parameters, "parity_sep");
    J_nn         = get_val<double>(parameters, "J_nn");
    J_nnn        = get_val<double>(parameters, "J_nnn");
    h_field      = get_val<double>(parameters, "h_field");
    h_rnd        = get_val<double>(parameters, "h_rnd");
    h_mean       = get_val<double>(parameters, "h_mean");
    h_sigma      = get_val<double>(parameters, "h_sigma");
    alpha        = get_val<double>(parameters, "alpha");
    beta         = get_val<double>(parameters, "beta");
    e_reduced    = get_val<double>(parameters, "e_reduced");
    psfactor     = get_val<double>(parameters, "psfactor");
    if(J_nnn != 0.0) throw std::runtime_error("Use of [J_nnn] - Next-nearest neighbor coupling - is not implemented yet");
    all_mpo_parameters_have_been_set = true;
}

class_tf_ising::Parameters class_tf_ising::get_parameters() const {
    /* clang-format off */
    Parameters parameters;
    parameters.push_back({"spin_dim", spin_dim});
    parameters.push_back({"position", get_position()});
    parameters.push_back({"distribution", distribution});
    parameters.push_back({"parity_sep", parity_sep});
    parameters.push_back({"J_nn", J_nn});
    parameters.push_back({"J_nnn", J_nnn});
    parameters.push_back({"h_field", h_field});
    parameters.push_back({"h_rnd", h_rnd});
    parameters.push_back({"h_mean", h_mean});
    parameters.push_back({"h_sigma", h_sigma});
    parameters.push_back({"alpha", alpha});
    parameters.push_back({"beta", beta});
    parameters.push_back({"e_reduced", e_reduced});
    parameters.push_back({"psfactor", psfactor});
    return parameters;
    /* clang-format on */
}

void class_tf_ising::register_h5_parameters() {}


void class_tf_ising::set_full_lattice_parameters([[maybe_unused]] std::vector<Parameters> lattice_parameters, bool reverse) {
    if(reverse) {
        std::reverse(lattice_parameters.begin(), lattice_parameters.end());
        for(size_t pos = 0; pos < lattice_parameters.size(); pos++)
        find_val(lattice_parameters[pos], "position") = pos;
    }
    find_val(lattice_parameters.back(), "J_nn")  = 0.0;
    find_val(lattice_parameters.back(), "J_nnn") = 0.0;
    // Recompute J_avg and h_avg from given J_rnd and h_rnd on all sites
    double J_sum = 0;
    double h_sum = 0;
    for(auto &site_params : lattice_parameters) {
        double J_nn_    = get_val<double>(site_params,"J_nn");
        double J_nnn_   = get_val<double>(site_params,"J_nnn");
        double h_field_ = get_val<double>(site_params,"h_field");
        double h_rnd_   = get_val<double>(site_params,"h_rnd");
        double h_ptb_   = get_val<double>(site_params,"h_ptb");
        double alpha_   = get_val<double>(site_params,"alpha");
        double beta_    = get_val<double>(site_params,"beta");
        J_sum += std::pow(J_nn_, 1.0 - alpha_) + J_nnn_;
        h_sum += h_field_ + std::pow(h_rnd_ + h_ptb_, 1 - beta_);
    }

    if(parity_sep or get_val<bool>(lattice_parameters[get_position()],"parity_sep")) {
        find_val(lattice_parameters[get_position()],"psfactor") = J_sum + h_sum;
    }
    set_parameters(lattice_parameters[get_position()]);
}
