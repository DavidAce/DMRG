//
// Created by david on 2018-07-04.
//

#include "class_ising_tf_rf.h"
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <h5pp/h5pp.h>
#include <iomanip>
#include <math/num.h>
#include <math/rnd.h>
#include <physics/nmspc_quantum_mechanics.h>

using namespace qm::spinHalf;
using Scalar = std::complex<double>;

class_ising_tf_rf::class_ising_tf_rf(ModelType model_type_, size_t position_) : class_mpo_site(model_type_, position_) {
    h5tb.param.J1       = settings::model::ising_tf_rf::J1;
    h5tb.param.J2       = settings::model::ising_tf_rf::J2;
    h5tb.param.h_tran   = settings::model::ising_tf_rf::h_tran;
    h5tb.param.h_mean   = settings::model::ising_tf_rf::h_mean;
    h5tb.param.h_stdv   = settings::model::ising_tf_rf::h_stdv;
    h5tb.param.spin_dim = settings::model::ising_tf_rf::spin_dim;
    copy_c_str(settings::model::ising_tf_rf::distribution, h5tb.param.distribution);
    extent4 = {1, 1, h5tb.param.spin_dim, h5tb.param.spin_dim};
    extent2 = {h5tb.param.spin_dim, h5tb.param.spin_dim};

    qm::spinHalf::SX = qm::gen_manybody_spins(sx, 2);
    qm::spinHalf::SY = qm::gen_manybody_spins(sy, 2);
    qm::spinHalf::SZ = qm::gen_manybody_spins(sz, 2);
    qm::spinHalf::II = qm::gen_manybody_spins(id, 2);

    class_ising_tf_rf::randomize_hamiltonian();
    h5tb_ising_tf_rf::register_table_type();
    all_mpo_parameters_have_been_set = true; // There are no full lattice parameters on this model so we set it true immediately!
    class_ising_tf_rf::build_mpo();
}

double class_ising_tf_rf::get_field() const { return h5tb.param.h_tran + std::pow(h5tb.param.h_pert + h5tb.param.h_rand, 1 - beta); }
double class_ising_tf_rf::get_coupling() const { return std::pow(h5tb.param.J1, 1 - alpha); }
void   class_ising_tf_rf::print_parameter_names() const { h5tb.print_parameter_names(); }
void   class_ising_tf_rf::print_parameter_values() const { h5tb.print_parameter_values(); }

void class_ising_tf_rf::set_parameters(TableMap &parameters) {
    h5tb.param.J1       = std::any_cast<double>(parameters["J1"]);
    h5tb.param.J2       = std::any_cast<double>(parameters["J2"]);
    h5tb.param.h_tran   = std::any_cast<double>(parameters["h_tran"]);
    h5tb.param.h_mean   = std::any_cast<double>(parameters["h_mean"]);
    h5tb.param.h_stdv   = std::any_cast<double>(parameters["h_stdv"]);
    h5tb.param.h_rand   = std::any_cast<double>(parameters["h_rand"]);
    h5tb.param.h_pert   = std::any_cast<double>(parameters["h_pert"]);
    h5tb.param.spin_dim = std::any_cast<long>(parameters["spin_dim"]);
    copy_c_str(std::any_cast<std::string>(parameters["distribution"]), h5tb.param.distribution);
    if(h5tb.param.J2 != 0.0) throw std::runtime_error("Use of [J2] - Next-nearest neighbor coupling - is not implemented yet");
    all_mpo_parameters_have_been_set = true;
    build_mpo();
}

class_ising_tf_rf::TableMap class_ising_tf_rf::get_parameters() const {
    /* clang-format off */
    TableMap parameters;
    parameters["J1"]            = h5tb.param.J1;
    parameters["J2"]            = h5tb.param.J2;
    parameters["h_tran"]        = h5tb.param.h_tran;
    parameters["h_mean"]        = h5tb.param.h_mean;
    parameters["h_stdv"]        = h5tb.param.h_stdv;
    parameters["h_rand"]        = h5tb.param.h_rand;
    parameters["h_pert"]        = h5tb.param.h_pert;
    parameters["spin_dim"]      = h5tb.param.spin_dim;
    parameters["distribution"]  = std::string(h5tb.param.distribution);
    return parameters;
    /* clang-format on */
}

void class_ising_tf_rf::build_mpo()
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
    mpo_internal.resize(3, 3, h5tb.param.spin_dim, h5tb.param.spin_dim);
    mpo_internal.setZero();
    mpo_internal.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(id);
    mpo_internal.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_field() * sx - e_reduced * id);
    mpo_internal.slice(Eigen::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_coupling() * sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 2, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(id);
    if(Textra::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw std::runtime_error(fmt::format("MPO at position {} has NAN's", get_position()));
    }
    build_mpo_squared();
}

Eigen::Tensor<Scalar, 1> class_ising_tf_rf::get_MPO_edge_left() const {
    Eigen::Tensor<Scalar, 1> ledge(3);
    ledge.setZero();
    ledge(2) = 1;
    return ledge;
}
Eigen::Tensor<Scalar, 1> class_ising_tf_rf::get_MPO_edge_right() const {
    Eigen::Tensor<Scalar, 1> redge(3);
    redge.setZero();
    redge(0) = 1;
    return redge;
}

Eigen::Tensor<Scalar, 1> class_ising_tf_rf::get_MPO2_edge_left() const {
    auto edge = get_MPO_edge_left();
    auto dim  = edge.dimension(0);
    return edge.contract(edge, Textra::idx()).reshape(Textra::array1{dim * dim});
}
Eigen::Tensor<Scalar, 1> class_ising_tf_rf::get_MPO2_edge_right() const {
    auto edge = get_MPO_edge_right();
    auto dim  = edge.dimension(0);
    return edge.contract(edge, Textra::idx()).reshape(Textra::array1{dim * dim});
}


void class_ising_tf_rf::randomize_hamiltonian() {
    if(std::string(h5tb.param.distribution) == "normal") {
        h5tb.param.h_rand = rnd::normal(h5tb.param.h_mean, h5tb.param.h_stdv);
    } else if(std::string(h5tb.param.distribution) == "lognormal") {
        h5tb.param.h_rand = rnd::log_normal(h5tb.param.h_mean, h5tb.param.h_stdv);
    } else if(std::string(h5tb.param.distribution) == "uniform") {
        h5tb.param.h_rand = rnd::uniform_double_box(h5tb.param.h_mean - h5tb.param.h_stdv / 2.0, h5tb.param.h_mean + h5tb.param.h_stdv / 2.0);
    } else {
        throw std::runtime_error("Wrong distribution given. Expected one of <normal>, <lognormal>, <uniform>");
    }
    all_mpo_parameters_have_been_set = false;
    mpo_squared                      = std::nullopt;
}

void class_ising_tf_rf::set_coupling_damping(double alpha_) { alpha = alpha_; }
void class_ising_tf_rf::set_field_damping(double beta_) {
    beta = beta_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_field() * sx - e_reduced * id);
        mpo_squared                                                                     = std::nullopt;
    }
}

void class_ising_tf_rf::set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) {
    switch(ptbMode) {
        case PerturbMode::ABSOLUTE: {
            h5tb.param.h_pert = field_ptb;
            break;
        }
        case PerturbMode::PERCENTAGE: {
            h5tb.param.h_pert = h5tb.param.h_rand * field_ptb;
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_ABSOLUTE: {
            h5tb.param.h_pert = rnd::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_PERCENTAGE: {
            h5tb.param.h_pert = h5tb.param.h_rand * rnd::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
    }
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_field() * sx - e_reduced * id);
        mpo_squared                                                                     = std::nullopt;
    }
    if(coupling_ptb == 0.0 and field_ptb == 0 and is_perturbed())
        throw std::runtime_error(fmt::format("MPO({}): Should have become unperturbed!", get_position()));
}

bool class_ising_tf_rf::is_perturbed() const { return h5tb.param.h_pert != 0.0; }


Eigen::Tensor<Scalar, 4> class_ising_tf_rf::MPO_nbody_view(const std::vector<size_t> &nbody_terms) const {
    // This function returns a view of the MPO including only n-body terms.
    // For instance, if nbody_terms == {2,3}, this would exclude on-site terms.
    if(nbody_terms.empty()) return MPO();
    double J1 = 0, J2 = 0;
    for(auto &&n : nbody_terms){
        if(n == 1) J1 = 1.0;
        if(n == 2) J2 = 1.0;
    }
    Eigen::Tensor<Scalar, 4> MPO_nbody = MPO();
    MPO_nbody.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-J1*get_field() * sx - e_reduced * id);
    MPO_nbody.slice(Eigen::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-J2*get_coupling() * sz);
    return MPO_nbody;
}


Eigen::Tensor<Scalar, 4> class_ising_tf_rf::MPO_reduced_view() const {
    if(e_reduced == 0) { return MPO(); }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> class_ising_tf_rf::MPO_reduced_view(double site_energy) const {
    if(site_energy == 0) { return MPO(); }
    Eigen::Tensor<Scalar, 4> temp                                           = MPO();
    temp.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_field() * sx - site_energy * id);
    return temp;
}

std::unique_ptr<class_mpo_site> class_ising_tf_rf::clone() const { return std::make_unique<class_ising_tf_rf>(*this); }

long class_ising_tf_rf::get_spin_dimension() const { return h5tb.param.spin_dim; }

void class_ising_tf_rf::set_averages([[maybe_unused]] std::vector<TableMap> lattice_parameters, bool infinite, bool reverse) {
    if(reverse) {
        std::reverse(lattice_parameters.begin(), lattice_parameters.end());
        for(size_t pos = 0; pos < lattice_parameters.size(); pos++) lattice_parameters[pos]["position"] = pos;
    }
    if(not infinite){
        lattice_parameters.back()["J1"] = 0.0;
        lattice_parameters.back()["J2"] = 0.0;
        lattice_parameters.end()[-2]["J2"] = 0.0;
    }
    // Recompute J_avrg and h_avrg from given J_rnrd and h_rnd on all sites
    double J_sum = 0;
    double h_sum = 0;
    for(auto &site_params : lattice_parameters) {
        auto J1_     = std::any_cast<double>(site_params["J1"]);
        auto J2_     = std::any_cast<double>(site_params["J2"]);
        auto h_tran_ = std::any_cast<double>(site_params["h_tran"]);
        auto h_rand_ = std::any_cast<double>(site_params["h_rand"]);
        J_sum += J1_ + J2_;
        h_sum += h_tran_ + h_rand_;
    }

    if(parity_sep) psfactor = J_sum + h_sum;
    set_parameters(lattice_parameters[get_position()]);
}

void class_ising_tf_rf::save_hamiltonian(h5pp::File &file, const std::string &table_path) const {
    if(not file.linkExists(table_path)) file.createTable(h5tb_ising_tf_rf::h5_type, table_path, "Transverse-field Ising");
    file.appendTableRecords(h5tb, table_path);
    // Position 0 is also responsible for writing attributes
    if(position.value() != 0) return;

    file.writeAttribute(h5tb.param.J1, "J1", table_path);
    file.writeAttribute(h5tb.param.J2, "J2", table_path);
    file.writeAttribute(h5tb.param.h_mean, "h_mean", table_path);
    file.writeAttribute(h5tb.param.h_stdv, "h_stdv", table_path);
    file.writeAttribute(h5tb.param.h_tran, "h_tran", table_path);
    file.writeAttribute(h5tb.param.distribution, "distribution", table_path);
    file.writeAttribute(h5tb.param.spin_dim, "spin_dim", table_path);
}

void class_ising_tf_rf::load_hamiltonian(const h5pp::File &file, const std::string &table_path) {
    std::string ham_prefix = table_path + "/hamiltonian";
    if(file.linkExists(ham_prefix)) {
        h5tb.param                       = file.readTableRecords<h5tb_ising_tf_rf::table>(ham_prefix, position);
        all_mpo_parameters_have_been_set = true;
        build_mpo();
    } else {
        throw std::runtime_error(fmt::format("Could not load MPO. Table [{}] does not exist", ham_prefix));
    }
    // We can use the mpo's on file here to check everything is correct
    std::string mpo_dset = fmt::format("{}/mpo/H_{}", table_path, get_position());
    if(file.linkExists(mpo_dset)) {
        if(Textra::Tensor_to_Vector(MPO()) != Textra::Tensor_to_Vector(file.readDataset<Eigen::Tensor<Scalar, 4>>(mpo_dset)))
            throw std::runtime_error("Built MPO does not match the MPO on file");
    }

    // Check that we are on the same point of the phase diagram
    if(std::abs(h5tb.param.J1 - settings::model::ising_tf_rf::J1) > 1e-6) throw std::runtime_error("J1 != settings::model::ising_tf_rf::J1");
    if(std::abs(h5tb.param.J2 - settings::model::ising_tf_rf::J2) > 1e-6) throw std::runtime_error("J2 != settings::model::ising_tf_rf::J2");
    if(std::abs(h5tb.param.h_tran - settings::model::ising_tf_rf::h_tran) > 1e-6) throw std::runtime_error("h_tran != settings::model::ising_tf_rf::h_tran");
    if(std::abs(h5tb.param.h_mean - settings::model::ising_tf_rf::h_mean) > 1e-6) throw std::runtime_error("h_mean != settings::model::ising_tf_rf::h_mean");
    if(std::abs(h5tb.param.h_stdv - settings::model::ising_tf_rf::h_stdv) > 1e-6) throw std::runtime_error("h_stdv != settings::model::ising_tf_rf::h_stdv");
}