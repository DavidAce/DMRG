//
// Created by david on 2018-07-04.
//

#include "class_lbit.h"
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <h5pp/h5pp.h>
#include <iomanip>
#include <math/num.h>
#include <math/rnd.h>
#include <physics/nmspc_quantum_mechanics.h>

using namespace qm::spinHalf;
using Scalar = std::complex<double>;

class_lbit::class_lbit(ModelType model_type_, size_t position_) : class_mpo_site(model_type_, position_) {
    h5tb.param.J1_mean  = settings::model::lbit::J1_mean;
    h5tb.param.J2_mean  = settings::model::lbit::J2_mean;
    h5tb.param.J3_mean  = settings::model::lbit::J3_mean;
    h5tb.param.J1_wdth  = settings::model::lbit::J1_wdth;
    h5tb.param.J2_wdth  = settings::model::lbit::J2_wdth;
    h5tb.param.J3_wdth  = settings::model::lbit::J3_wdth;
    h5tb.param.spin_dim = settings::model::lbit::spin_dim;
    copy_c_str(settings::model::lbit::distribution,h5tb.param.distribution);
    extent4 = {1, 1, h5tb.param.spin_dim, h5tb.param.spin_dim};
    extent2 = {h5tb.param.spin_dim, h5tb.param.spin_dim};
    h5tb_lbit::register_table_type();
    all_mpo_parameters_have_been_set = false; // There are no full lattice parameters but we set it to false so we remember to call randomize on all sites in every model type
}

double class_lbit::get_field() const { return  + std::pow(h5tb.param.J1_rand, 1 - beta); }
double class_lbit::get_coupling() const { return std::pow(h5tb.param.J2_rand + h5tb.param.J3_rand, 1 - alpha); }
void   class_lbit::print_parameter_names() const { h5tb_lbit::print_parameter_names(); }
void   class_lbit::print_parameter_values() const { h5tb.print_parameter_values(); }

void class_lbit::set_parameters(TableMap &parameters) {
    h5tb.param.J1_rand  = std::any_cast<double>(parameters["J1_rand"]);
    h5tb.param.J2_rand  = std::any_cast<double>(parameters["J2_rand"]);
    h5tb.param.J3_rand  = std::any_cast<double>(parameters["J3_rand"]);
    h5tb.param.J1_wdth  = std::any_cast<double>(parameters["J1_wdth"]);
    h5tb.param.J2_wdth  = std::any_cast<double>(parameters["J2_wdth"]);
    h5tb.param.J3_wdth  = std::any_cast<double>(parameters["J3_wdth"]);
    h5tb.param.spin_dim = std::any_cast<long>(parameters["spin_dim"]);
    std::strcpy(h5tb.param.distribution, std::any_cast<std::string>(parameters["distribution"]).c_str());
    all_mpo_parameters_have_been_set = true;
    build_mpo();
}

class_lbit::TableMap class_lbit::get_parameters() const {
    /* clang-format off */
    TableMap parameters;
    parameters["J1_rand"]            = h5tb.param.J1_rand;
    parameters["J2_rand"]            = h5tb.param.J2_rand;
    parameters["J3_rand"]            = h5tb.param.J3_rand;
    parameters["J1_wdth"]            = h5tb.param.J1_wdth;
    parameters["J2_wdth"]            = h5tb.param.J2_wdth;
    parameters["J3_wdth"]            = h5tb.param.J3_wdth;
    parameters["spin_dim"]      = h5tb.param.spin_dim;
    parameters["distribution"]  = std::string(h5tb.param.distribution);
    return parameters;
    /* clang-format on */
}

void class_lbit::build_mpo()
/*! Builds the MPO hamiltonian as a rank 4 tensor. Notation following Schollwöck (2010)

 * H = Σ J1_rand n_i  + Σ J2_rand n_{i} n_{i+1} + Σ J3_rand n_{i} n_{i+1} n_{i+2}
 *
 * where n_i = 0.5 * (1 + σ_i^z),  and σ_i^z is the diagonal 2x2 pauli matrix
 *
 *       2            |     I      0     0     0   |
 *       |            |     n      0     0     0   |
 *   0---H---1    =   |     0      n     0     0   |
 *       |            |   J1_rand*n   J2_rand*n  J3_rand*n    I   |
 *       3
 *
 *        2
 *        |
 *    0---H---1
 *        |
 *        3
 */
{
    if(not all_mpo_parameters_have_been_set) throw std::runtime_error("Improperly built MPO: Full lattice parameters haven't been set yet.");
    Eigen::Tensor<Scalar,2> n = Textra::MatrixToTensor(0.5*(id+sz));
    Eigen::Tensor<Scalar,2> i = Textra::MatrixTensorMap(id);
    mpo_internal.resize(4, 4, h5tb.param.spin_dim, h5tb.param.spin_dim);
    mpo_internal.setZero();
    mpo_internal.slice(Textra::array4{0, 0, 0, 0}, extent4).reshape(extent2) = i;
    mpo_internal.slice(Textra::array4{1, 0, 0, 0}, extent4).reshape(extent2) = n;
    mpo_internal.slice(Textra::array4{2, 1, 0, 0}, extent4).reshape(extent2) = n;
    mpo_internal.slice(Textra::array4{3, 0, 0, 0}, extent4).reshape(extent2) = h5tb.param.J1_rand * n - e_reduced * i;
    mpo_internal.slice(Textra::array4{3, 1, 0, 0}, extent4).reshape(extent2) = h5tb.param.J2_rand * n;
    mpo_internal.slice(Textra::array4{3, 2, 0, 0}, extent4).reshape(extent2) = h5tb.param.J3_rand * n;
    mpo_internal.slice(Textra::array4{3, 3, 0, 0}, extent4).reshape(extent2) = i;
    if(Textra::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw std::runtime_error(fmt::format("MPO at position {} has NAN's", get_position()));
    }
    unique_id = std::nullopt;
    build_mpo_squared();
}


Eigen::Tensor<Scalar, 1> class_lbit::get_MPO_edge_left() const {
    Eigen::Tensor<Scalar, 1> ledge(4);
    ledge.setZero();
    ledge(3) = 1;
    return ledge;
}

Eigen::Tensor<Scalar, 1> class_lbit::get_MPO_edge_right() const {
    Eigen::Tensor<Scalar, 1> redge(4);
    redge.setZero();
    redge(0) = 1;
    return redge;
}

Eigen::Tensor<Scalar, 1> class_lbit::get_MPO2_edge_left() const {
    auto edge = get_MPO_edge_left();
    auto dim  = edge.dimension(0);
    return edge.contract(edge, Textra::idx()).reshape(Textra::array1{dim * dim});
}
Eigen::Tensor<Scalar, 1> class_lbit::get_MPO2_edge_right() const {
    auto edge = get_MPO_edge_right();
    auto dim  = edge.dimension(0);
    return edge.contract(edge, Textra::idx()).reshape(Textra::array1{dim * dim});
}



void class_lbit::randomize_hamiltonian() {
    if(std::string(h5tb.param.distribution) == "normal") {
        h5tb.param.J1_rand = rnd::normal(settings::model::lbit::J1_mean, settings::model::lbit::J1_wdth);
        h5tb.param.J2_rand = rnd::normal(settings::model::lbit::J2_mean, settings::model::lbit::J2_wdth);
        h5tb.param.J3_rand = rnd::normal(settings::model::lbit::J3_mean, settings::model::lbit::J3_wdth);
    } else if(std::string(h5tb.param.distribution) == "lognormal") {
        h5tb.param.J1_rand =  rnd::log_normal(settings::model::lbit::J1_mean, settings::model::lbit::J1_wdth);
        h5tb.param.J2_rand =  rnd::log_normal(settings::model::lbit::J2_mean, settings::model::lbit::J2_wdth);
        h5tb.param.J3_rand =  rnd::log_normal(settings::model::lbit::J3_mean, settings::model::lbit::J3_wdth);
    } else if(std::string(h5tb.param.distribution) == "uniform") {
        h5tb.param.J1_rand = settings::model::lbit::J1_mean + rnd::uniform_double_box(settings::model::lbit::J1_wdth);
        h5tb.param.J2_rand = settings::model::lbit::J2_mean + rnd::uniform_double_box(settings::model::lbit::J2_wdth);
        h5tb.param.J3_rand = settings::model::lbit::J3_mean + rnd::uniform_double_box(settings::model::lbit::J3_wdth);
    }else if(std::string(h5tb.param.distribution) == "constant") {
            h5tb.param.J1_rand = settings::model::lbit::J1_mean;
            h5tb.param.J2_rand = settings::model::lbit::J2_mean;
            h5tb.param.J3_rand = settings::model::lbit::J3_mean;
    } else {
        throw std::runtime_error("Wrong distribution given. Expected one of <normal>, <lognormal>, <uniform> or <constant>");
    }
    all_mpo_parameters_have_been_set = false;
    mpo_squared                      = std::nullopt;
}

void class_lbit::set_coupling_damping(double alpha_) { alpha = alpha_; }
void class_lbit::set_field_damping(double beta_) {
    beta = beta_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_field() * sx - e_reduced * id);
        mpo_squared                                                                     = std::nullopt;
        unique_id                                                                       = std::nullopt;
        unique_id_sq                                                                    = std::nullopt;
    }
}

void class_lbit::set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) {
    switch(ptbMode) {
        case PerturbMode::ABSOLUTE: {
            h5tb.param.J1_pert = field_ptb;
            h5tb.param.J2_pert = coupling_ptb;
            h5tb.param.J3_pert = coupling_ptb;
            break;
        }
        case PerturbMode::PERCENTAGE: {
            h5tb.param.J1_pert = h5tb.param.J1_pert * field_ptb;
            h5tb.param.J2_pert = h5tb.param.J2_pert * coupling_ptb;
            h5tb.param.J3_pert = h5tb.param.J3_pert * coupling_ptb;
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_ABSOLUTE: {
            h5tb.param.J1_pert = rnd::uniform_double_box(-field_ptb, field_ptb);
            h5tb.param.J2_pert = rnd::uniform_double_box(-coupling_ptb, coupling_ptb);
            h5tb.param.J3_pert = rnd::uniform_double_box(-coupling_ptb, coupling_ptb);
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_PERCENTAGE: {
            h5tb.param.J1_pert = h5tb.param.J1_pert * rnd::uniform_double_box(-field_ptb, field_ptb);
            h5tb.param.J2_pert = h5tb.param.J2_pert * rnd::uniform_double_box(-coupling_ptb, coupling_ptb);
            h5tb.param.J3_pert = h5tb.param.J3_pert * rnd::uniform_double_box(-coupling_ptb, coupling_ptb);
            break;
        }
    }
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_field() * sx - e_reduced * id);
        mpo_squared                                                                     = std::nullopt;
        unique_id                                                                       = std::nullopt;
        unique_id_sq                                                                    = std::nullopt;
    }
    if(coupling_ptb == 0.0 and field_ptb == 0 and is_perturbed())
        throw std::runtime_error(fmt::format("MPO({}): Should have become unperturbed!", get_position()));
}

bool class_lbit::is_perturbed() const { return h5tb.param.J1_pert != 0.0 or h5tb.param.J2_pert != 0.0 or h5tb.param.J3_pert != 0.0; }


Eigen::Tensor<Scalar, 4> class_lbit::MPO_nbody_view(const std::vector<size_t> &nbody_terms) const {
    // This function returns a view of the MPO including only n-body terms.
    // For instance, if nbody_terms == {2,3}, this would exclude on-site terms.
    if(nbody_terms.empty()) return MPO();
    double J1 = 0;
    double J2 = 0;
    double J3 = 0;
    for(const auto & n : nbody_terms){
        if(n == 1) J1 = 1;
        if(n == 2) J2 = 1;
        if(n == 3) J3 = 1;
    }
    Eigen::Tensor<Scalar, 4> MPO_nbody = MPO();
    Eigen::Tensor<Scalar,2> n = Textra::MatrixToTensor(0.5*(id+sz));
    Eigen::Tensor<Scalar,2> i = Textra::MatrixTensorMap(id);
    MPO_nbody.slice(Textra::array4{3, 0, 0, 0}, extent4).reshape(extent2) = J1 * h5tb.param.J1_rand * n - e_reduced * i;
    MPO_nbody.slice(Textra::array4{3, 1, 0, 0}, extent4).reshape(extent2) = J2 * h5tb.param.J2_rand * n;
    MPO_nbody.slice(Textra::array4{3, 2, 0, 0}, extent4).reshape(extent2) = J3 * h5tb.param.J3_rand * n;
    return MPO_nbody;
}

Eigen::Tensor<Scalar, 4> class_lbit::MPO_reduced_view() const {
    if(e_reduced == 0) { return MPO(); }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> class_lbit::MPO_reduced_view(double site_energy) const {
    if(site_energy == 0) { return MPO(); }
    Eigen::Tensor<Scalar, 4> temp                                           = MPO();
    temp.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_field() * sx - site_energy * id);
    return temp;
}

std::unique_ptr<class_mpo_site> class_lbit::clone() const { return std::make_unique<class_lbit>(*this); }

long class_lbit::get_spin_dimension() const { return h5tb.param.spin_dim; }

void class_lbit::set_averages([[maybe_unused]] std::vector<TableMap> lattice_parameters, bool infinite, bool reverse) {
    if(reverse) {
        std::reverse(lattice_parameters.begin(), lattice_parameters.end());
        for(size_t pos = 0; pos < lattice_parameters.size(); pos++) lattice_parameters[pos]["position"] = pos;
    }
    if(not infinite){
        lattice_parameters.back()["J2_rand"] = 0.0;
        lattice_parameters.back()["J3_rand"] = 0.0;
        lattice_parameters.end()[-2]["J3_rand"] = 0.0;
    }

    double J_sum = 0;
    for(auto &site_params : lattice_parameters) {
        auto J1_     = std::any_cast<double>(site_params["J1_rand"]);
        auto J2_     = std::any_cast<double>(site_params["J2_rand"]);
        auto J3_     = std::any_cast<double>(site_params["J3_rand"]);
        J_sum += J1_ + J2_ + J3_;
    }
    if(parity_sep) psfactor = J_sum;
    set_parameters(lattice_parameters[get_position()]);
}

void class_lbit::save_hamiltonian(h5pp::File &file, const std::string &table_path) const {
    if(not file.linkExists(table_path)) file.createTable(h5tb_lbit::h5_type, table_path, "Transverse-field Ising");
    file.appendTableRecords(h5tb, table_path);
    // Position 0 is also responsible for writing attributes
    if(position.value() != 0) return;

    file.writeAttribute(h5tb.param.J1_mean, "J1_mean", table_path);
    file.writeAttribute(h5tb.param.J2_mean, "J2_mean", table_path);
    file.writeAttribute(h5tb.param.J3_mean, "J3_mean", table_path);
    file.writeAttribute(h5tb.param.J1_wdth, "J1_wdth", table_path);
    file.writeAttribute(h5tb.param.J2_wdth, "J2_wdth", table_path);
    file.writeAttribute(h5tb.param.J3_wdth, "J3_wdth", table_path);
    file.writeAttribute(h5tb.param.distribution, "distribution", table_path);
    file.writeAttribute(h5tb.param.spin_dim, "spin_dim", table_path);
}

void class_lbit::load_hamiltonian(const h5pp::File &file, const std::string &model_prefix) {
    auto ham_table = fmt::format("{}/hamiltonian", model_prefix);
    if(file.linkExists(ham_table)) {
        h5tb.param                       = file.readTableRecords<h5tb_lbit::table>(ham_table, position);
        all_mpo_parameters_have_been_set = true;
        build_mpo();
    } else {
        throw std::runtime_error(fmt::format("Could not load MPO. Table [{}] does not exist", ham_table));
    }

    // Check that we are on the same point of the phase diagram
    using namespace settings::model::lbit;
    if(std::abs(h5tb.param.J1_mean - J1_mean) > 1e-6) throw std::runtime_error(fmt::format("J1_mean {:.16f} != {:.16f} lbit::J1_mean", h5tb.param.J1_rand, J1_mean));
    if(std::abs(h5tb.param.J2_mean - J2_mean) > 1e-6) throw std::runtime_error(fmt::format("J2_mean {:.16f} != {:.16f} lbit::J2_mean", h5tb.param.J2_rand, J2_mean));
    if(std::abs(h5tb.param.J3_mean - J3_mean) > 1e-6) throw std::runtime_error(fmt::format("J3_mean {:.16f} != {:.16f} lbit::J3_mean", h5tb.param.J3_rand, J3_mean));
    if(std::abs(h5tb.param.J1_wdth - J1_wdth) > 1e-6) throw std::runtime_error(fmt::format("J1_wdth {:.16f} != {:.16f} lbit::J1_wdth", h5tb.param.J1_wdth, J1_wdth));
    if(std::abs(h5tb.param.J2_wdth - J2_wdth) > 1e-6) throw std::runtime_error(fmt::format("J2_wdth {:.16f} != {:.16f} lbit::J2_wdth", h5tb.param.J2_wdth, J2_wdth));
    if(std::abs(h5tb.param.J3_wdth - J3_wdth) > 1e-6) throw std::runtime_error(fmt::format("J3_wdth {:.16f} != {:.16f} lbit::J3_wdth", h5tb.param.J3_wdth, J3_wdth));

    // We can use the mpo's on file here to check everything is correct
    std::string mpo_dset = fmt::format("{}/mpo/H_{}", model_prefix, get_position());
    if(file.linkExists(mpo_dset)) {
        if(Textra::Tensor_to_Vector(MPO()) != Textra::Tensor_to_Vector(file.readDataset<Eigen::Tensor<Scalar, 4>>(mpo_dset)))
            throw std::runtime_error("Built MPO does not match the MPO on file");
    }
}