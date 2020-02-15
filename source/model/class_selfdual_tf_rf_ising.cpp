//
// Created by david on 2018-07-06.
//

#include "class_selfdual_tf_rf_ising.h"
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <simulation/nmspc_settings.h>

using namespace qm::spinOneHalf;
using Scalar = std::complex<double>;

class_selfdual_tf_rf_ising::class_selfdual_tf_rf_ising(size_t position_) : class_model_base(position_) {
    log          = Logger::setLogger("selfdual-tf-rf-ising");
    spin_dim     = settings::model::selfdual_tf_rf_ising::d; /*!< Spin dimension */
    J_mean       = settings::model::selfdual_tf_rf_ising::J_mean;
    h_mean       = settings::model::selfdual_tf_rf_ising::h_mean;
    J_sigma      = settings::model::selfdual_tf_rf_ising::J_sigma;
    h_sigma      = settings::model::selfdual_tf_rf_ising::h_sigma;
    lambda       = settings::model::selfdual_tf_rf_ising::lambda;
    delta        = J_mean - h_mean;
    parity_sep   = settings::model::selfdual_tf_rf_ising::parity_sep;
    distribution = settings::model::selfdual_tf_rf_ising::distribution;
    e_reduced    = 0.0;
    extent4      = {1, 1, spin_dim, spin_dim};
    extent2      = {spin_dim, spin_dim};
    class_selfdual_tf_rf_ising::randomize_hamiltonian();
}

double class_selfdual_tf_rf_ising::get_coupling() const { return std::pow(J_rnd + J_ptb, 1 - alpha); }
double class_selfdual_tf_rf_ising::get_field() const { return std::pow(h_rnd + h_ptb, 1 - beta); }
double class_selfdual_tf_rf_ising::get_coupling(double J_rnd_, double J_ptb_, double alpha_) const { return std::pow(J_rnd_ + J_ptb_, 1 - alpha_); }
double class_selfdual_tf_rf_ising::get_field(double h_rnd_, double h_ptb_, double beta_) const { return std::pow(h_rnd_ + h_ptb_, 1 - beta_); }

void class_selfdual_tf_rf_ising::set_parameters(const Parameters &parameters) {
    spin_dim                         = get_val<size_t>(parameters, "spin_dim");
    position                         = get_val<size_t>(parameters, "position");
    distribution                     = get_val<std::string>(parameters, "distribution");
    parity_sep                       = get_val<bool>(parameters, "parity_sep");
    J_rnd                            = get_val<double>(parameters, "J_rnd");
    J_ptb                            = get_val<double>(parameters, "J_ptb");
    h_rnd                            = get_val<double>(parameters, "h_rnd");
    h_ptb                            = get_val<double>(parameters, "h_ptb");
    J_avg                            = get_val<double>(parameters, "J_avg");
    h_avg                            = get_val<double>(parameters, "h_avg");
    J_mean                           = get_val<double>(parameters, "J_mean");
    h_mean                           = get_val<double>(parameters, "h_mean");
    J_sigma                          = get_val<double>(parameters, "J_sigma");
    h_sigma                          = get_val<double>(parameters, "h_sigma");
    lambda                           = get_val<double>(parameters, "lambda");
    delta                            = get_val<double>(parameters, "delta");
    alpha                            = get_val<double>(parameters, "alpha");
    beta                             = get_val<double>(parameters, "beta");
    e_reduced                        = get_val<double>(parameters, "e_reduced");
    psfactor                         = get_val<double>(parameters, "psfactor");
    all_mpo_parameters_have_been_set = true;
    build_mpo();
}

class_selfdual_tf_rf_ising::Parameters class_selfdual_tf_rf_ising::get_parameters() const {
    /* clang-format off */
    Parameters parameters;
    parameters.push_back({"position", get_position()});
    parameters.push_back({"J_rnd", J_rnd});
    parameters.push_back({"h_rnd", h_rnd});
    parameters.push_back({"J_ptb", J_ptb});
    parameters.push_back({"h_ptb", h_ptb});
    parameters.push_back({"J_avg", J_avg});
    parameters.push_back({"h_avg", h_avg});
    parameters.push_back({"J_mean", J_mean});
    parameters.push_back({"h_mean", h_mean});
    parameters.push_back({"J_sigma", J_sigma});
    parameters.push_back({"h_sigma", h_sigma});
    parameters.push_back({"lambda", lambda});
    parameters.push_back({"delta", delta});
    parameters.push_back({"alpha", alpha});
    parameters.push_back({"beta", beta});
    parameters.push_back({"e_reduced", e_reduced});
    parameters.push_back({"spin_dim", spin_dim});
    parameters.push_back({"distribution", distribution});
    parameters.push_back({"parity_sep", parity_sep});
    parameters.push_back({"psfactor", psfactor});
    return parameters;
    /* clang-format on */
}

void class_selfdual_tf_rf_ising::register_h5_parameters() {}
//
//    H5ParameterType = H5Tcreate(H5T_COMPOUND, sizeof(Particle));
//    H5Tinsert(MY_HDF5_PARTICLE_TYPE, "x", HOFFSET(Particle, x), H5T_NATIVE_DOUBLE);
//    H5Tinsert(MY_HDF5_PARTICLE_TYPE, "y", HOFFSET(Particle, y), H5T_NATIVE_DOUBLE);
//    H5Tinsert(MY_HDF5_PARTICLE_TYPE, "z", HOFFSET(Particle, z), H5T_NATIVE_DOUBLE);
//    H5Tinsert(MY_HDF5_PARTICLE_TYPE, "t", HOFFSET(Particle, t), H5T_NATIVE_DOUBLE);
//    H5Tinsert(MY_HDF5_PARTICLE_TYPE, "name", HOFFSET(Particle, name), MY_HDF5_NAME_TYPE);
//
//
//}

void class_selfdual_tf_rf_ising::build_mpo()
/*! Builds the MPO hamiltonian as a rank 4 tensor. Notation following Schollwöck (2010)
 *
 * H = - Σ J_{i} sz_{i} sz_{i+1} +  h_{i} sx_{i} + l*(h_avg sx_i sx_{i+1} + J_avg sz_{i} sz_{i+2})
 *
 *  |     I                 0           0              0            0   |
 *  |     sz                0           0              0            0   |
 *  |     sx                0           0              0            0   |
 *  |     0                 I           0              0            0   |
 *  | -(h_rnd)*sx       -J_rnd*sz   -l*h_avg*sx     -l*J_avg*sz     I   |
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
    if(parity_sep) {
        mpo_internal.resize(6, 6, spin_dim, spin_dim);
        mpo_internal.setZero();
        mpo_internal.slice(Eigen::array<long, 4>{5, 5, 0, 0}, extent4).reshape(extent2) =
            Textra::MatrixTensorMap(sx); // Multiply the psfactor on the edge! Not on each MPO!
    } else {
        mpo_internal.resize(5, 5, spin_dim, spin_dim);
        mpo_internal.setZero();
    }

    mpo_internal.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sx);
    mpo_internal.slice(Eigen::array<long, 4>{3, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(lambda * h_avg) * sx);
    mpo_internal.slice(Eigen::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(lambda * J_avg) * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 4, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);

    if(Textra::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw std::runtime_error("MPO at position " + std::to_string(get_position()) + " has been constructed with NAN's");
    }
}

void class_selfdual_tf_rf_ising::randomize_hamiltonian() {
    if(distribution == "normal") {
        J_rnd = rn::normal(J_mean, J_sigma);
        h_rnd = rn::normal(h_mean, h_sigma);
    } else if(distribution == "lognormal") {
        J_rnd = rn::log_normal(J_mean, J_sigma);
        h_rnd = rn::log_normal(h_mean, h_sigma);
    } else if(distribution == "uniform") {
        J_rnd = rn::uniform_double_box(J_mean - J_sigma / 2.0, J_mean + J_sigma / 2.0);
        h_rnd = rn::uniform_double_box(h_mean - h_sigma / 2.0, h_mean + J_sigma / 2.0);
    } else {
        throw std::runtime_error("Wrong distribution given. Expected one of <normal>, <lognormal>, <uniform>");
    }

    all_mpo_parameters_have_been_set = false;
}

void class_selfdual_tf_rf_ising::set_coupling_damping(double alpha_) {
    alpha = alpha_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    }
}
void class_selfdual_tf_rf_ising::set_field_damping(double beta_) {
    beta = beta_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    }
}

void class_selfdual_tf_rf_ising::set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) {
    switch(ptbMode) {
        case PerturbMode::ABSOLUTE: {
            J_ptb = coupling_ptb;
            h_ptb = field_ptb;
            break;
        }
        case PerturbMode::PERCENTAGE: {
            J_ptb = J_rnd * coupling_ptb;
            h_ptb = h_rnd * field_ptb;
            break;
        }

        case PerturbMode::UNIFORM_RANDOM_ABSOLUTE: {
            J_ptb = rn::uniform_double_box(-coupling_ptb, coupling_ptb);
            h_ptb = rn::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_PERCENTAGE: {
            J_ptb = J_rnd * rn::uniform_double_box(-coupling_ptb, coupling_ptb);
            h_ptb = h_rnd * rn::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
    }
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
        mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    }
    if(coupling_ptb == 0.0 and field_ptb == 0 and is_perturbed())
        throw std::runtime_error("MPO(" + std::to_string(get_position()) + ": Should have become unperturbed!");
}

bool class_selfdual_tf_rf_ising::is_perturbed() const { return J_ptb != 0.0 or h_ptb != 0.0; }
bool class_selfdual_tf_rf_ising::is_damped() const { return alpha != 0.0 or beta != 0.0; }

Eigen::Tensor<Scalar, 4> class_selfdual_tf_rf_ising::MPO_reduced_view() const {
    if(e_reduced == 0) {
        return MPO();
    }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> class_selfdual_tf_rf_ising::MPO_reduced_view(double site_energy) const {
    if(site_energy == 0) {
        return MPO();
    }
    Eigen::Tensor<Scalar, 4> temp                                           = MPO();
    temp.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - site_energy * Id);
    temp.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    return temp;
}

Eigen::Tensor<Scalar, 1> class_selfdual_tf_rf_ising::get_MPO_edge_left() const {
    if(parity_sep) {
        Eigen::Tensor<Scalar, 1> ledge(6);
        ledge.setZero();
        ledge(4) = 1;
        ledge(5) = psfactor;
        return ledge;
    } else {
        Eigen::Tensor<Scalar, 1> ledge(5);
        ledge.setZero();
        ledge(4) = 1;
        return ledge;
    }
}

Eigen::Tensor<Scalar, 1> class_selfdual_tf_rf_ising::get_MPO_edge_right() const {
    if(parity_sep) {
        Eigen::Tensor<Scalar, 1> redge(6);
        redge.setZero();
        redge(0) = 1;
        redge(5) = 1;
        return redge;
    } else {
        Eigen::Tensor<Scalar, 1> redge(5);
        redge.setZero();
        redge(0) = 1;
        return redge;
    }
}

Eigen::MatrixXcd class_selfdual_tf_rf_ising::single_site_hamiltonian(int position, int sites, std::vector<Eigen::MatrixXcd> &SX,
                                                                     std::vector<Eigen::MatrixXcd> &SY [[maybe_unused]],
                                                                     std::vector<Eigen::MatrixXcd> &SZ) const {
    int i = math::mod(position, sites);
    int j = math::mod(position + 1, sites);
    int k = math::mod(position + 2, sites);
    return -(J_rnd * SZ[i] * SZ[j] + h_rnd * 0.5 * (SX[i] + SX[j]) + lambda * (h_avg * SX[i] * SX[j] + J_avg * SZ[i] * SZ[k]));
}

std::unique_ptr<class_model_base> class_selfdual_tf_rf_ising::clone() const { return std::make_unique<class_selfdual_tf_rf_ising>(*this); }

size_t class_selfdual_tf_rf_ising::get_spin_dimension() const { return spin_dim; }

void class_selfdual_tf_rf_ising::set_full_lattice_parameters(std::vector<Parameters> all_parameters, bool reverse) {
    if(reverse) {
        std::reverse(all_parameters.begin(), all_parameters.end());
        for(size_t pos = 0; pos < all_parameters.size(); pos++) {
            find_val(all_parameters[pos], "position") = pos;
            if(pos < all_parameters.size() - 1) {
                find_val(all_parameters[pos], "J_rnd") = find_val(all_parameters[pos + 1], "J_rnd");
                find_val(all_parameters[pos], "J_ptb") = find_val(all_parameters[pos + 1], "J_ptb");
            } else {
                find_val(all_parameters[pos], "J_rnd") = 0.0;
                find_val(all_parameters[pos], "J_ptb") = 0.0;
            }
        }
    } else {
        find_val(all_parameters.back(), "J_rnd") = 0.0;
        find_val(all_parameters.back(), "J_ptb") = 0.0;
    }

    // Recompute J_avg and h_avg from given J_rnd and h_rnd on all sites
    double J_avg_ = 0;
    double h_avg_ = 0;
    for(auto &site_params : all_parameters) {
        J_avg_ += get_val<double>(site_params, "J_rnd");
        h_avg_ += get_val<double>(site_params, "h_rnd");
    }
    J_avg_ /= (all_parameters.size() - 1);
    h_avg_ /= (all_parameters.size());
    for(auto &site_params : all_parameters) {
        find_val(site_params, "J_avg") = J_avg_;
        find_val(site_params, "h_avg") = h_avg_;
    }
    if(parity_sep or get_val<bool>(all_parameters[get_position()], "parity_sep")) {
        find_val(all_parameters[get_position()], "psfactor") = all_parameters.size() * (J_avg_ + h_avg_) * (1.0 + lambda);
    }
    set_parameters(all_parameters[get_position()]);
}
