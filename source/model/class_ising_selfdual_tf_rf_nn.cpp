//
// Created by david on 2018-07-06.
//

#include "class_ising_selfdual_tf_rf_nn.h"
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <h5pp/h5pp.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <simulation/nmspc_settings.h>

using namespace qm::spinOneHalf;
using Scalar = std::complex<double>;

class_ising_selfdual_tf_rf_nn::class_ising_selfdual_tf_rf_nn(std::string_view model_type_, size_t position_) : class_model_base(model_type_, position_) {
    pm.J_mean   = settings::model::ising_selfdual_tf_rf_nn::J_mean;
    pm.h_mean   = settings::model::ising_selfdual_tf_rf_nn::h_mean;
    pm.J_stdv   = settings::model::ising_selfdual_tf_rf_nn::J_stdv;
    pm.h_stdv   = settings::model::ising_selfdual_tf_rf_nn::h_stdv;
    pm.lambda   = settings::model::ising_selfdual_tf_rf_nn::lambda;
    pm.delta    = pm.J_mean - pm.h_mean;
    pm.spin_dim = settings::model::ising_selfdual_tf_rf_nn::spin_dim;
    std::strcpy(pm.distribution, settings::model::ising_selfdual_tf_rf_nn::distribution.c_str());
    parity_sep = settings::model::ising_selfdual_tf_rf_nn::parity_sep;

    extent4 = {1, 1, pm.spin_dim, pm.spin_dim};
    extent2 = {pm.spin_dim, pm.spin_dim};
    class_ising_selfdual_tf_rf_nn::randomize_hamiltonian();
    h5tb_ising_selfdual_tf_rf_nn::register_table_type();
}

double class_ising_selfdual_tf_rf_nn::get_coupling() const { return std::pow(pm.J_rand + pm.J_pert, 1 - alpha); }
double class_ising_selfdual_tf_rf_nn::get_field() const { return std::pow(pm.h_rand + pm.h_pert, 1 - beta); }
double class_ising_selfdual_tf_rf_nn::get_coupling(double J_rnd_, double J_ptb_, double alpha_) const { return std::pow(J_rnd_ + J_ptb_, 1 - alpha_); }
double class_ising_selfdual_tf_rf_nn::get_field(double h_rnd_, double h_ptb_, double beta_) const { return std::pow(h_rnd_ + h_ptb_, 1 - beta_); }

void class_ising_selfdual_tf_rf_nn::set_parameters(TableMap &parameters) {
    pm.J_rand   = std::any_cast<double>(parameters["J_rand"]);
    pm.J_pert   = std::any_cast<double>(parameters["J_pert"]);
    pm.h_rand   = std::any_cast<double>(parameters["h_rand"]);
    pm.h_pert   = std::any_cast<double>(parameters["h_pert"]);
    pm.J_avrg   = std::any_cast<double>(parameters["J_avrg"]);
    pm.h_avrg   = std::any_cast<double>(parameters["h_avrg"]);
    pm.J_mean   = std::any_cast<double>(parameters["J_mean"]);
    pm.h_mean   = std::any_cast<double>(parameters["h_mean"]);
    pm.J_stdv   = std::any_cast<double>(parameters["J_stdv"]);
    pm.h_stdv   = std::any_cast<double>(parameters["h_stdv"]);
    pm.lambda   = std::any_cast<double>(parameters["lambda"]);
    pm.delta    = std::any_cast<double>(parameters["delta"]);
    pm.spin_dim = std::any_cast<size_t>(parameters["spin_dim"]);
    std::strcpy(pm.distribution, std::any_cast<std::string>(parameters["distribution"]).c_str());
    all_mpo_parameters_have_been_set = true;
    build_mpo();
}

class_ising_selfdual_tf_rf_nn::TableMap class_ising_selfdual_tf_rf_nn::get_parameters() const {
    /* clang-format off */
    TableMap parameters;
    parameters["J_rand"] = pm.J_rand;
    parameters["h_rand"] = pm.h_rand;
    parameters["J_pert"] = pm.J_pert;
    parameters["h_pert"] = pm.h_pert;
    parameters["J_avrg"] = pm.J_avrg;
    parameters["h_avrg"] = pm.h_avrg;
    parameters["J_mean"] = pm.J_mean;
    parameters["h_mean"] = pm.h_mean;
    parameters["J_stdv"] = pm.J_stdv;
    parameters["h_stdv"] = pm.h_stdv;
    parameters["lambda"] = pm.lambda;
    parameters["delta"] =  pm.delta;
    parameters["spin_dim"] = pm.spin_dim;
    parameters["distribution"] = std::string(pm.distribution);
    return parameters;
    /* clang-format on */
}

void class_ising_selfdual_tf_rf_nn::build_mpo()
/*! Builds the MPO hamiltonian as a rank 4 tensor. Notation following Schollwöck (2010)
 *
 * H = - Σ J_{i} sz_{i} sz_{i+1} +  h_{i} sx_{i} + l*(h_avrg sx_i sx_{i+1} + J_avrg sz_{i} sz_{i+2})
 *
 *  |     I                 0           0              0            0   |
 *  |     sz                0           0              0            0   |
 *  |     sx                0           0              0            0   |
 *  |     0                 I           0              0            0   |
 *  | -(h_rnd)*sx       -J_rnd*sz   -l*h_av g*sx    -l*J_avrg*sz    I   |
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
        mpo_internal.resize(6, 6, pm.spin_dim, pm.spin_dim);
        mpo_internal.setZero();
        mpo_internal.slice(Eigen::array<long, 4>{5, 5, 0, 0}, extent4).reshape(extent2) =
            Textra::MatrixTensorMap(sx); // Multiply the psfactor on the edge! Not on each MPO!
    } else {
        mpo_internal.resize(5, 5, pm.spin_dim, pm.spin_dim);
        mpo_internal.setZero();
    }

    mpo_internal.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sx);
    mpo_internal.slice(Eigen::array<long, 4>{3, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(pm.lambda * pm.h_avrg) * sx);
    mpo_internal.slice(Eigen::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(pm.lambda * pm.J_avrg) * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 4, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);

    if(Textra::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw std::runtime_error("MPO at position " + std::to_string(get_position()) + " has been constructed with NAN's");
    }
}

void class_ising_selfdual_tf_rf_nn::randomize_hamiltonian() {
    if(std::string(pm.distribution) == "normal") {
        pm.J_rand = rn::normal(pm.J_mean, pm.J_stdv);
        pm.h_rand = rn::normal(pm.h_mean, pm.h_stdv);
    } else if(std::string(pm.distribution) == "lognormal") {
        pm.J_rand = rn::log_normal(pm.J_mean, pm.J_stdv);
        pm.h_rand = rn::log_normal(pm.h_mean, pm.h_stdv);
    } else if(std::string(pm.distribution) == "uniform") {
        pm.J_rand = rn::uniform_double_box(pm.J_mean - pm.J_stdv / 2.0, pm.J_mean + pm.J_stdv / 2.0);
        pm.h_rand = rn::uniform_double_box(pm.h_mean - pm.h_stdv / 2.0, pm.h_mean + pm.h_stdv / 2.0);
    } else {
        throw std::runtime_error("Wrong distribution given. Expected one of <normal>, <lognormal>, <uniform>");
    }

    all_mpo_parameters_have_been_set = false;
}

void class_ising_selfdual_tf_rf_nn::set_coupling_damping(double alpha_) {
    alpha = alpha_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    }
}
void class_ising_selfdual_tf_rf_nn::set_field_damping(double beta_) {
    beta = beta_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    }
}

void class_ising_selfdual_tf_rf_nn::set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) {
    switch(ptbMode) {
        case PerturbMode::ABSOLUTE: {
            pm.J_pert = coupling_ptb;
            pm.h_pert = field_ptb;
            break;
        }
        case PerturbMode::PERCENTAGE: {
            pm.J_pert = pm.J_rand * coupling_ptb;
            pm.h_pert = pm.h_rand * field_ptb;
            break;
        }

        case PerturbMode::UNIFORM_RANDOM_ABSOLUTE: {
            pm.J_pert = rn::uniform_double_box(-coupling_ptb, coupling_ptb);
            pm.h_pert = rn::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_PERCENTAGE: {
            pm.J_pert = pm.J_rand * rn::uniform_double_box(-coupling_ptb, coupling_ptb);
            pm.h_pert = pm.h_rand * rn::uniform_double_box(-field_ptb, field_ptb);
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

bool class_ising_selfdual_tf_rf_nn::is_perturbed() const { return pm.J_pert != 0.0 or pm.h_pert != 0.0; }

Eigen::Tensor<Scalar, 4> class_ising_selfdual_tf_rf_nn::MPO_reduced_view() const {
    if(e_reduced == 0) {
        return MPO();
    }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> class_ising_selfdual_tf_rf_nn::MPO_reduced_view(double site_energy) const {
    if(site_energy == 0) {
        return MPO();
    }
    Eigen::Tensor<Scalar, 4> temp                                           = MPO();
    temp.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - site_energy * Id);
    temp.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    return temp;
}

Eigen::Tensor<Scalar, 1> class_ising_selfdual_tf_rf_nn::get_MPO_edge_left() const {
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

Eigen::Tensor<Scalar, 1> class_ising_selfdual_tf_rf_nn::get_MPO_edge_right() const {
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

Eigen::MatrixXcd class_ising_selfdual_tf_rf_nn::single_site_hamiltonian(size_t position, size_t sites, std::vector<Eigen::MatrixXcd> &SX,
                                                                     std::vector<Eigen::MatrixXcd> &SY [[maybe_unused]],
                                                                     std::vector<Eigen::MatrixXcd> &SZ) const {
    auto i = math::mod(position, sites);
    auto j = math::mod(position + 1, sites);
    auto k = math::mod(position + 2, sites);
    return -(pm.J_rand * SZ[i] * SZ[j] + pm.h_rand * 0.5 * (SX[i] + SX[j]) + pm.lambda * (pm.h_avrg * SX[i] * SX[j] + pm.J_avrg * SZ[i] * SZ[k]));
}

std::unique_ptr<class_model_base> class_ising_selfdual_tf_rf_nn::clone() const { return std::make_unique<class_ising_selfdual_tf_rf_nn>(*this); }

size_t class_ising_selfdual_tf_rf_nn::get_spin_dimension() const { return pm.spin_dim; }

void class_ising_selfdual_tf_rf_nn::set_averages(std::vector<TableMap> all_parameters, bool reverse) {
    if(reverse) {
        // We need to reverse the parameters, and move them one step
        std::reverse(all_parameters.begin(), all_parameters.end());
        for(size_t pos = 0; pos < all_parameters.size(); pos++) {
            all_parameters[pos]["position"] = pos;
            if(pos < all_parameters.size() - 1) {
                all_parameters[pos]["J_rnd"] = pos < all_parameters.size() - 1 ? all_parameters[pos + 1]["J_rnd"] : 0.0;
                all_parameters[pos]["J_ptb"] = pos < all_parameters.size() - 1 ? all_parameters[pos + 1]["J_ptb"] : 0.0;
            }
        }
    } else {
        all_parameters.back()["J_rnd"] = 0.0;
        all_parameters.back()["J_ptb"] = 0.0;
    }

    // Recompute J_avrg and pm.h_avrg from given pm.J_rand and pm.h_rand on all sites
    double J_avrg_ = 0;
    double h_avrg_ = 0;
    for(auto &site_param : all_parameters) {
        J_avrg_ += std::any_cast<double>(site_param["J_rand"]);
        h_avrg_ += std::any_cast<double>(site_param["h_rand"]);
    }
    J_avrg_ /= static_cast<double>(all_parameters.size() - 1);
    h_avrg_ /= static_cast<double>(all_parameters.size());

    for(auto &site_params : all_parameters) {
        site_params["J_avrg"] = J_avrg_;
        site_params["h_avrg"] = h_avrg_;
    }
    if(parity_sep) psfactor = (J_avrg_ + h_avrg_) * (1.0 + pm.lambda) * static_cast<double>(all_parameters.size());
    set_parameters(all_parameters[get_position()]);
}

void class_ising_selfdual_tf_rf_nn::write_hamiltonian(h5pp::File &file, const std::string &model_prefix) const {
    std::string ham_path = model_prefix + "/Hamiltonian";
    if(not file.linkExists(ham_path)) file.createTable(h5tb_ising_selfdual_tf_rf_nn::h5_type, ham_path, "Selfdual Ising");
    file.appendTableEntries(pm, ham_path);
}

void class_ising_selfdual_tf_rf_nn::read_hamiltonian(const h5pp::File &file, const std::string &model_prefix) {
    std::string ham_path = model_prefix + "/Hamiltonian";
    if(file.linkExists(ham_path)) {
        pm                               = file.readTableEntries<h5tb_ising_selfdual_tf_rf_nn::table>(ham_path, position);
        all_mpo_parameters_have_been_set = true;
        build_mpo();
    }else{
        throw std::runtime_error(fmt::format("Could not load MPO. Table [{}] does not exist", ham_path));
    }
    // We can use the mpo's on file here to check everything is correct
    std::string mpo_dset = model_prefix + "/mpo/H_" + std::to_string(get_position());
    if(file.linkExists(mpo_dset)) {
        if(Textra::Tensor_to_Vector(MPO()) != Textra::Tensor_to_Vector(file.readDataset<Eigen::Tensor<Scalar,4>>(mpo_dset)))
            throw std::runtime_error("Built MPO does not match the MPO on file");
    }

    // Check that we are on the same point of the phase diagram
    if(std::abs(pm.J_mean - settings::model::ising_selfdual_tf_rf_nn::J_mean) > 1e-6)   throw std::runtime_error("J_mean != settings::model::ising_selfdual_tf_rf_nn::J_mean");
    if(std::abs(pm.h_mean - settings::model::ising_selfdual_tf_rf_nn::h_mean) > 1e-6)   throw std::runtime_error("h_mean != settings::model::ising_selfdual_tf_rf_nn::h_mean");
    if(std::abs(pm.J_stdv - settings::model::ising_selfdual_tf_rf_nn::J_stdv) > 1e-6)   throw std::runtime_error("J_stdv != settings::model::ising_selfdual_tf_rf_nn::J_stdv");
    if(std::abs(pm.h_stdv - settings::model::ising_selfdual_tf_rf_nn::h_stdv) > 1e-6)   throw std::runtime_error("h_stdv != settings::model::ising_selfdual_tf_rf_nn::h_stdv");
    if(std::abs(pm.lambda - settings::model::ising_selfdual_tf_rf_nn::lambda) > 1e-6)   throw std::runtime_error("lambda != settings::model::ising_selfdual_tf_rf_nn::lambda");
    if(pm.distribution != settings::model::ising_selfdual_tf_rf_nn::distribution) throw std::runtime_error("distribution != settings::model::ising_selfdual_tf_rf_nn::distribution");
}

