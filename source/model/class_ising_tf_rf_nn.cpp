//
// Created by david on 2018-07-04.
//

#include "class_ising_tf_rf_nn.h"
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <h5pp/h5pp.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <simulation/nmspc_settings.h>

using namespace qm::spinOneHalf;
using Scalar = std::complex<double>;

class_ising_tf_rf_nn::class_ising_tf_rf_nn(std::string_view model_type_, size_t position_) : class_model_base(model_type_,position_) {
    pm.J1       = settings::model::ising_tf_rf_nn::J1;
    pm.J2       = settings::model::ising_tf_rf_nn::J2;
    pm.h_tran   = settings::model::ising_tf_rf_nn::h_tran;
    pm.h_mean  = settings::model::ising_tf_rf_nn::h_mean;
    pm.h_stdv   = settings::model::ising_tf_rf_nn::h_stdv;
    pm.spin_dim = settings::model::ising_tf_rf_nn::spin_dim;
    std::strcpy(pm.distribution, settings::model::ising_tf_rf_nn::distribution.c_str());
    extent4 = {1, 1, pm.spin_dim, pm.spin_dim};
    extent2 = {pm.spin_dim, pm.spin_dim};

    qm::spinOneHalf::SX = qm::gen_manybody_spin(sx, 2);
    qm::spinOneHalf::SY = qm::gen_manybody_spin(sy, 2);
    qm::spinOneHalf::SZ = qm::gen_manybody_spin(sz, 2);
    qm::spinOneHalf::II = qm::gen_manybody_spin(Id, 2);

    class_ising_tf_rf_nn::randomize_hamiltonian();
    h5tb_ising_tf_rf_nn::register_table_type();
    all_mpo_parameters_have_been_set = true; // There are no full lattice parameters on this model so we set it true immediately!
    class_ising_tf_rf_nn::build_mpo();
}

double class_ising_tf_rf_nn::get_field() const { return pm.h_tran + std::pow(pm.h_pert + pm.h_rand, 1 - beta); }
double class_ising_tf_rf_nn::get_coupling() const { return std::pow(pm.J1, 1 - alpha); }

void class_ising_tf_rf_nn::build_mpo()
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
    mpo_internal.resize(3, 3, pm.spin_dim, pm.spin_dim);
    mpo_internal.setZero();
    mpo_internal.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    mpo_internal.slice(Eigen::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 2, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
}

void class_ising_tf_rf_nn::randomize_hamiltonian() {
    if(std::string(pm.distribution) == "normal") {
        pm.h_rand = rn::normal(pm.h_mean, pm.h_stdv);
    } else if(std::string(pm.distribution) == "lognormal") {
        pm.h_rand = rn::log_normal(pm.h_mean, pm.h_stdv);
    } else if(std::string(pm.distribution) == "uniform") {
        pm.h_rand = rn::uniform_double_box(pm.h_mean - pm.h_stdv / 2.0, pm.h_mean + pm.h_stdv / 2.0);
    } else {
        throw std::runtime_error("Wrong distribution given. Expected one of <normal>, <lognormal>, <uniform>");
    }
    all_mpo_parameters_have_been_set = false;
}

void class_ising_tf_rf_nn::set_coupling_damping(double alpha_) { alpha = alpha_; }
void class_ising_tf_rf_nn::set_field_damping(double beta_) {
    beta = beta_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    }
}

void class_ising_tf_rf_nn::set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) {
    switch(ptbMode) {
        case PerturbMode::ABSOLUTE: {
            pm.h_pert = field_ptb;
            break;
        }
        case PerturbMode::PERCENTAGE: {
            pm.h_pert = pm.h_rand * field_ptb;
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_ABSOLUTE: {
            pm.h_pert = rn::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_PERCENTAGE: {
            pm.h_pert = pm.h_rand * rn::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
    }
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    }
    if(coupling_ptb == 0.0 and field_ptb == 0 and is_perturbed())
        throw std::runtime_error("MPO(" + std::to_string(get_position()) + ": Should have become unperturbed!");
}

bool class_ising_tf_rf_nn::is_perturbed() const { return pm.h_pert != 0.0; }

Eigen::Tensor<Scalar, 4> class_ising_tf_rf_nn::MPO_reduced_view() const {
    if(e_reduced == 0) {
        return MPO();
    }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> class_ising_tf_rf_nn::MPO_reduced_view(double site_energy) const {
    if(site_energy == 0) {
        return MPO();
    }
    Eigen::Tensor<Scalar, 4> temp                                           = MPO();
    temp.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - site_energy * Id);
    return temp;
}

Eigen::MatrixXcd class_ising_tf_rf_nn::single_site_hamiltonian(size_t position, size_t sites, std::vector<Eigen::MatrixXcd> &SX,
                                                         std::vector<Eigen::MatrixXcd> &SY [[maybe_unused]], std::vector<Eigen::MatrixXcd> &SZ) const {
    auto i = math::mod(position, sites);
    auto j = math::mod(position + 1, sites);
    return -(pm.J1 * SZ[i] * SZ[j] + pm.h_tran * 0.5 * (SX[i] + SX[j]));
}

std::unique_ptr<class_model_base> class_ising_tf_rf_nn::clone() const { return std::make_unique<class_ising_tf_rf_nn>(*this); }

size_t class_ising_tf_rf_nn::get_spin_dimension() const { return pm.spin_dim; }

Eigen::Tensor<Scalar, 1> class_ising_tf_rf_nn::get_MPO_edge_left() const {
    Eigen::Tensor<Scalar, 1> ledge(3);
    ledge.setZero();
    ledge(2) = 1;
    return ledge;
}
Eigen::Tensor<Scalar, 1> class_ising_tf_rf_nn::get_MPO_edge_right() const {
    Eigen::Tensor<Scalar, 1> redge(3);
    redge.setZero();
    redge(0) = 1;
    return redge;
}

void class_ising_tf_rf_nn::set_parameters(TableMap &parameters) {
    pm.J1       = std::any_cast<double>(parameters["J1"]);
    pm.J2       = std::any_cast<double>(parameters["J2"]);
    pm.h_tran   = std::any_cast<double>(parameters["h_tran"]);
    pm.h_mean   = std::any_cast<double>(parameters["h_mean"]);
    pm.h_stdv   = std::any_cast<double>(parameters["h_stdv"]);
    pm.h_rand   = std::any_cast<double>(parameters["h_rand"]);
    pm.h_pert   = std::any_cast<double>(parameters["h_pert"]);
    pm.spin_dim = std::any_cast<size_t>(parameters["spin_dim"]);
    std::strcpy(pm.distribution, std::any_cast<std::string>(parameters["distribution"]).c_str());
    if(pm.J2 != 0.0) throw std::runtime_error("Use of [J2] - Next-nearest neighbor coupling - is not implemented yet");
    all_mpo_parameters_have_been_set = true;
}

class_ising_tf_rf_nn::TableMap class_ising_tf_rf_nn::get_parameters() const {
    /* clang-format off */
    TableMap parameters;
    parameters["J1"]            = pm.J1;
    parameters["J2"]            = pm.J2;
    parameters["h_tran"]        = pm.h_tran;
    parameters["h_mean"]        = pm.h_mean;
    parameters["h_stdv"]        = pm.h_stdv;
    parameters["h_rand"]        = pm.h_rand;
    parameters["h_pert"]        = pm.h_pert;
    parameters["spin_dim"]      = pm.spin_dim;
    parameters["distribution"]  = std::string(pm.distribution);
    return parameters;
    /* clang-format on */
}

void class_ising_tf_rf_nn::set_averages([[maybe_unused]] std::vector<TableMap> lattice_parameters, bool reverse) {
    if(reverse) {
        std::reverse(lattice_parameters.begin(), lattice_parameters.end());
        for(size_t pos = 0; pos < lattice_parameters.size(); pos++) lattice_parameters[pos]["position"] = pos;
    }
    lattice_parameters.back()["J1"]  = 0.0;
    lattice_parameters.back()["J2"] = 0.0;
    // Recompute J_avrg and h_avrg from given J_rnrd and h_rnd on all sites
    double J_sum = 0;
    double h_sum = 0;
    for(auto &site_params : lattice_parameters) {
        auto J1_    = std::any_cast<double>(site_params["J1"]);
        auto J2_    = std::any_cast<double>(site_params["J2"]);
        auto h_tran_ = std::any_cast<double>(site_params["h_tran"]);
        auto h_rand_   = std::any_cast<double>(site_params["h_rand"]);
        J_sum += J1_ + J2_;
        h_sum += h_tran_ + h_rand_;
    }

    if(parity_sep) psfactor = J_sum + h_sum;
    set_parameters(lattice_parameters[get_position()]);
}




void class_ising_tf_rf_nn::write_hamiltonian(h5pp::File &file, const std::string &model_prefix) const {
    std::string ham_path = model_prefix + "/Hamiltonian";
    if(not file.linkExists(ham_path)) file.createTable(h5tb_ising_tf_rf_nn::h5_type, ham_path, "Transverse-field Ising");
    file.appendTableEntries(pm, ham_path);
}

void class_ising_tf_rf_nn::read_hamiltonian(const h5pp::File &file, const std::string &model_prefix) {
    std::string ham_path = model_prefix + "/Hamiltonian";
    if(file.linkExists(ham_path)) {
        pm                               = file.readTableEntries<h5tb_ising_tf_rf_nn::table>(ham_path, position);
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
    if(std::abs(pm.J1 - settings::model::ising_tf_rf_nn::J1) > 1e-6)   throw std::runtime_error("J1 != settings::model::ising_tf_rf_nn::J1");
    if(std::abs(pm.J2 - settings::model::ising_tf_rf_nn::J2) > 1e-6)   throw std::runtime_error("J2 != settings::model::ising_tf_rf_nn::J2");
    if(std::abs(pm.h_tran - settings::model::ising_tf_rf_nn::h_tran) > 1e-6) throw std::runtime_error("h_tran != settings::model::ising_tf_rf_nn::h_tran");
    if(std::abs(pm.h_mean - settings::model::ising_tf_rf_nn::h_mean) > 1e-6) throw std::runtime_error("h_mean != settings::model::ising_tf_rf_nn::h_mean");
    if(std::abs(pm.h_stdv - settings::model::ising_tf_rf_nn::h_stdv) > 1e-6) throw std::runtime_error("h_stdv != settings::model::ising_tf_rf_nn::h_stdv");
}