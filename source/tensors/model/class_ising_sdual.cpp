//
// Created by david on 2018-07-06.
//

#include "class_ising_sdual.h"
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <h5pp/h5pp.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <config/nmspc_settings.h>

using namespace qm::spinOneHalf;
using Scalar = std::complex<double>;

class_ising_sdual::class_ising_sdual(ModelType model_type_, size_t position_) : class_mpo_site(model_type_, position_) {
    h5tb.param.J_mean   = settings::model::ising_sdual::J_mean;
    h5tb.param.h_mean   = settings::model::ising_sdual::h_mean;
    h5tb.param.J_stdv   = settings::model::ising_sdual::J_stdv;
    h5tb.param.h_stdv   = settings::model::ising_sdual::h_stdv;
    h5tb.param.lambda   = settings::model::ising_sdual::lambda;
    h5tb.param.delta    = h5tb.param.J_mean - h5tb.param.h_mean;
    h5tb.param.spin_dim = settings::model::ising_sdual::spin_dim;
    std::strcpy(h5tb.param.distribution, settings::model::ising_sdual::distribution.c_str());
    parity_sep = settings::model::ising_sdual::parity_sep;

    extent4 = {1, 1, h5tb.param.spin_dim, h5tb.param.spin_dim};
    extent2 = {h5tb.param.spin_dim, h5tb.param.spin_dim};
    class_ising_sdual::randomize_hamiltonian();
    h5tb_ising_sdual::register_table_type();
}

double class_ising_sdual::get_coupling() const { return std::pow(h5tb.param.J_rand + h5tb.param.J_pert, 1 - alpha); }
double class_ising_sdual::get_field() const { return std::pow(h5tb.param.h_rand + h5tb.param.h_pert, 1 - beta); }
double class_ising_sdual::get_coupling(double J_rnd_, double J_ptb_, double alpha_) const { return std::pow(J_rnd_ + J_ptb_, 1 - alpha_); }
double class_ising_sdual::get_field(double h_rnd_, double h_ptb_, double beta_) const { return std::pow(h_rnd_ + h_ptb_, 1 - beta_); }
void class_ising_sdual::print_parameter_names() const { h5tb.print_parameter_names(); }
void class_ising_sdual::print_parameter_values() const { h5tb.print_parameter_values(); }

void class_ising_sdual::set_parameters(TableMap &parameters) {
    h5tb.param.J_mean   = std::any_cast<double>(parameters["J_mean"]);
    h5tb.param.J_stdv   = std::any_cast<double>(parameters["J_stdv"]);
    h5tb.param.J_rand   = std::any_cast<double>(parameters["J_rand"]);
    h5tb.param.J_avrg   = std::any_cast<double>(parameters["J_avrg"]);
    h5tb.param.J_pert   = std::any_cast<double>(parameters["J_pert"]);
    h5tb.param.h_mean   = std::any_cast<double>(parameters["h_mean"]);
    h5tb.param.h_stdv   = std::any_cast<double>(parameters["h_stdv"]);
    h5tb.param.h_rand   = std::any_cast<double>(parameters["h_rand"]);
    h5tb.param.h_avrg   = std::any_cast<double>(parameters["h_avrg"]);
    h5tb.param.h_pert   = std::any_cast<double>(parameters["h_pert"]);
    h5tb.param.lambda   = std::any_cast<double>(parameters["lambda"]);
    h5tb.param.delta    = std::any_cast<double>(parameters["delta"]);
    h5tb.param.spin_dim = std::any_cast<long>(parameters["spin_dim"]);
    std::strcpy(h5tb.param.distribution, std::any_cast<std::string>(parameters["distribution"]).c_str());
    all_mpo_parameters_have_been_set = true;
    build_mpo();
}

class_ising_sdual::TableMap class_ising_sdual::get_parameters() const {
    /* clang-format off */
    TableMap parameters;
    parameters["J_mean"]   = h5tb.param.J_mean;
    parameters["J_stdv"]   = h5tb.param.J_stdv;
    parameters["J_rand"]   = h5tb.param.J_rand;
    parameters["J_avrg"]   = h5tb.param.J_avrg;
    parameters["J_pert"]   = h5tb.param.J_pert;
    parameters["h_mean"]   = h5tb.param.h_mean;
    parameters["h_stdv"]   = h5tb.param.h_stdv;
    parameters["h_rand"]   = h5tb.param.h_rand;
    parameters["h_avrg"]   = h5tb.param.h_avrg;
    parameters["h_pert"]   = h5tb.param.h_pert;
    parameters["lambda"]   = h5tb.param.lambda;
    parameters["delta"]    = h5tb.param.delta;
    parameters["spin_dim"] = h5tb.param.spin_dim;
    parameters["distribution"] = std::string(h5tb.param.distribution);
    return parameters;
    /* clang-format on */
}

void class_ising_sdual::build_mpo()
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
        mpo_internal.resize(6, 6, h5tb.param.spin_dim, h5tb.param.spin_dim);
        mpo_internal.setZero();
        mpo_internal.slice(Eigen::array<long, 4>{5, 5, 0, 0}, extent4).reshape(extent2) =
            Textra::MatrixTensorMap(sx); // Multiply the psfactor on the edge! Not on each MPO!
    } else {
        mpo_internal.resize(5, 5, h5tb.param.spin_dim, h5tb.param.spin_dim);
        mpo_internal.setZero();
    }

    mpo_internal.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sx);
    mpo_internal.slice(Eigen::array<long, 4>{3, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(h5tb.param.lambda * h5tb.param.h_avrg) * sx);
    mpo_internal.slice(Eigen::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(h5tb.param.lambda * h5tb.param.J_avrg) * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 4, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);

    if(Textra::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw std::runtime_error("MPO at position " + std::to_string(get_position()) + " has been constructed with NAN's");
    }
}

void class_ising_sdual::randomize_hamiltonian() {
    if(std::string(h5tb.param.distribution) == "normal") {
        h5tb.param.J_rand = rn::normal(h5tb.param.J_mean, h5tb.param.J_stdv);
        h5tb.param.h_rand = rn::normal(h5tb.param.h_mean, h5tb.param.h_stdv);
    } else if(std::string(h5tb.param.distribution) == "lognormal") {
        h5tb.param.J_rand = rn::log_normal(h5tb.param.J_mean, h5tb.param.J_stdv);
        h5tb.param.h_rand = rn::log_normal(h5tb.param.h_mean, h5tb.param.h_stdv);
    } else if(std::string(h5tb.param.distribution) == "uniform") {
        h5tb.param.J_rand = rn::uniform_double_box(h5tb.param.J_mean - h5tb.param.J_stdv / 2.0, h5tb.param.J_mean + h5tb.param.J_stdv / 2.0);
        h5tb.param.h_rand = rn::uniform_double_box(h5tb.param.h_mean - h5tb.param.h_stdv / 2.0, h5tb.param.h_mean + h5tb.param.h_stdv / 2.0);
    } else {
        throw std::runtime_error("Wrong distribution given. Expected one of <normal>, <lognormal>, <uniform>");
    }

    all_mpo_parameters_have_been_set = false;
}

void class_ising_sdual::set_coupling_damping(double alpha_) {
    alpha = alpha_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    }
}
void class_ising_sdual::set_field_damping(double beta_) {
    beta = beta_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    }
}

void class_ising_sdual::set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) {
    switch(ptbMode) {
        case PerturbMode::ABSOLUTE: {
            h5tb.param.J_pert = coupling_ptb;
            h5tb.param.h_pert = field_ptb;
            break;
        }
        case PerturbMode::PERCENTAGE: {
            h5tb.param.J_pert = h5tb.param.J_rand * coupling_ptb;
            h5tb.param.h_pert = h5tb.param.h_rand * field_ptb;
            break;
        }

        case PerturbMode::UNIFORM_RANDOM_ABSOLUTE: {
            h5tb.param.J_pert = rn::uniform_double_box(-coupling_ptb, coupling_ptb);
            h5tb.param.h_pert = rn::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_PERCENTAGE: {
            h5tb.param.J_pert = h5tb.param.J_rand * rn::uniform_double_box(-coupling_ptb, coupling_ptb);
            h5tb.param.h_pert = h5tb.param.h_rand * rn::uniform_double_box(-field_ptb, field_ptb);
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

bool class_ising_sdual::is_perturbed() const { return h5tb.param.J_pert != 0.0 or h5tb.param.h_pert != 0.0; }

Eigen::Tensor<Scalar, 4> class_ising_sdual::MPO_reduced_view() const {
    if(e_reduced == 0) {
        return MPO();
    }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> class_ising_sdual::MPO_reduced_view(double site_energy) const {
    if(site_energy == 0) {
        return MPO();
    }
    Eigen::Tensor<Scalar, 4> temp                                           = MPO();
    temp.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - site_energy * Id);
    temp.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    return temp;
}

Eigen::Tensor<Scalar, 1> class_ising_sdual::get_MPO_edge_left() const {
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

Eigen::Tensor<Scalar, 1> class_ising_sdual::get_MPO_edge_right() const {
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

Eigen::MatrixXcd class_ising_sdual::single_site_hamiltonian(size_t position, size_t sites, std::vector<Eigen::MatrixXcd> &SX,
                                                            std::vector<Eigen::MatrixXcd> &SY [[maybe_unused]], std::vector<Eigen::MatrixXcd> &SZ) const {
    auto i = math::mod(position, sites);
    auto j = math::mod(position + 1, sites);
    auto k = math::mod(position + 2, sites);
    return -(h5tb.param.J_rand * SZ[i] * SZ[j] + h5tb.param.h_rand * 0.5 * (SX[i] + SX[j]) +
             h5tb.param.lambda * (h5tb.param.h_avrg * SX[i] * SX[j] + h5tb.param.J_avrg * SZ[i] * SZ[k]));
}

std::unique_ptr<class_mpo_site> class_ising_sdual::clone() const { return std::make_unique<class_ising_sdual>(*this); }

long class_ising_sdual::get_spin_dimension() const { return h5tb.param.spin_dim; }

void class_ising_sdual::set_averages(std::vector<TableMap> all_parameters, bool reverse) {
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

    // Recompute J_avrg and pm.param.h_avrg from given pm.param.J_rand and pm.param.h_rand on all sites
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
    if(parity_sep) psfactor = (J_avrg_ + h_avrg_) * (1.0 + h5tb.param.lambda) * static_cast<double>(all_parameters.size());
    set_parameters(all_parameters[get_position()]);
}

void class_ising_sdual::write_hamiltonian(h5pp::File &file, const std::string &model_prefix) const {
    std::string ham_path = model_prefix + "/Hamiltonian";
    if(not file.linkExists(ham_path)) file.createTable(h5tb_ising_sdual::h5_type, ham_path, "Selfdual Ising");
    file.appendTableEntries(h5tb, ham_path);
}

void class_ising_sdual::read_hamiltonian(const h5pp::File &file, const std::string &model_prefix) {
    std::string ham_path = model_prefix + "/Hamiltonian";
    if(file.linkExists(ham_path)) {
        h5tb.param                       = file.readTableEntries<h5tb_ising_sdual::table>(ham_path, position);
        all_mpo_parameters_have_been_set = true;
        build_mpo();
    } else {
        throw std::runtime_error(fmt::format("Could not load MPO. Table [{}] does not exist", ham_path));
    }
    // We can use the mpo's on file here to check everything is correct
    std::string mpo_dset = model_prefix + "/mpo/H_" + std::to_string(get_position());
    if(file.linkExists(mpo_dset)) {
        if(Textra::Tensor_to_Vector(MPO()) != Textra::Tensor_to_Vector(file.readDataset<Eigen::Tensor<Scalar, 4>>(mpo_dset)))
            throw std::runtime_error("Built MPO does not match the MPO on file");
    }

    // Check that we are on the same point of the phase diagram
    if(std::abs(h5tb.param.J_mean - settings::model::ising_sdual::J_mean) > 1e-6) throw std::runtime_error("J_mean != settings::model::ising_sdual::J_mean");
    if(std::abs(h5tb.param.h_mean - settings::model::ising_sdual::h_mean) > 1e-6) throw std::runtime_error("h_mean != settings::model::ising_sdual::h_mean");
    if(std::abs(h5tb.param.J_stdv - settings::model::ising_sdual::J_stdv) > 1e-6) throw std::runtime_error("J_stdv != settings::model::ising_sdual::J_stdv");
    if(std::abs(h5tb.param.h_stdv - settings::model::ising_sdual::h_stdv) > 1e-6) throw std::runtime_error("h_stdv != settings::model::ising_sdual::h_stdv");
    if(std::abs(h5tb.param.lambda - settings::model::ising_sdual::lambda) > 1e-6) throw std::runtime_error("lambda != settings::model::ising_sdual::lambda");
    if(h5tb.param.distribution != settings::model::ising_sdual::distribution)
        throw std::runtime_error("distribution != settings::model::ising_sdual::distribution");
}
