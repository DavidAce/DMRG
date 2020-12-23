//
// Created by david on 2018-07-06.
//

#include "class_ising_sdual.h"
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <h5pp/h5pp.h>
#include <iomanip>
#include <math/num.h>
#include <math/rnd.h>
#include <physics/nmspc_quantum_mechanics.h>

using namespace qm::spinHalf;
using Scalar = std::complex<double>;

double delta_to_J_mean(double delta) { return delta > 0 ? 1.0 : std::exp(delta); }

double delta_to_h_mean(double delta) { return delta > 0 ? std::exp(-delta) : 1.0; }

class_ising_sdual::class_ising_sdual(ModelType model_type_, size_t position_) : class_mpo_site(model_type_, position_) {
    h5tb.param.lambda = settings::model::ising_sdual::lambda;
    h5tb.param.delta  = settings::model::ising_sdual::delta;
    h5tb.param.J_mean = delta_to_J_mean(h5tb.param.delta);
    h5tb.param.h_mean = delta_to_h_mean(h5tb.param.delta);
    h5tb.param.J_stdv = settings::model::ising_sdual::J_stdv;
    h5tb.param.h_stdv = settings::model::ising_sdual::h_stdv;

    // Sanity check on delta, J_mean, h_mean
    double delta_check = std::log(h5tb.param.J_mean) - std::log(h5tb.param.h_mean);
    if(std::abs(h5tb.param.delta - delta_check) > 1e-10)
        throw std::logic_error(
            fmt::format("Error when transforming delta to (J_mean, h_mean): delta {:.12f} != {:.16f} delta_check", h5tb.param.delta, delta_check));

    h5tb.param.spin_dim = settings::model::ising_sdual::spin_dim;
    copy_c_str(settings::model::ising_sdual::distribution,h5tb.param.distribution);
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
void   class_ising_sdual::print_parameter_names() const { h5tb.print_parameter_names(); }
void   class_ising_sdual::print_parameter_values() const { h5tb.print_parameter_values(); }

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
    copy_c_str(std::any_cast<std::string>(parameters["distribution"]),h5tb.param.distribution);
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
 *  | -(h_rnd)*sx       -J_rnd*sz   -l*h_avrg*sx    -l*J_avrg*sz    I   |
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

    mpo_internal.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(id);
    mpo_internal.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sx);
    mpo_internal.slice(Eigen::array<long, 4>{3, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_field() * sx - e_reduced * id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_coupling() * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-(h5tb.param.lambda * h5tb.param.h_avrg) * sx);
    mpo_internal.slice(Eigen::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-(h5tb.param.lambda * h5tb.param.J_avrg) * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 4, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(id);
    if(Textra::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw std::runtime_error(fmt::format("MPO at position {} has NAN's", get_position()));
    }
    unique_id    = std::nullopt;
    build_mpo_squared();
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
Eigen::Tensor<Scalar, 1> class_ising_sdual::get_MPO2_edge_left() const {
    auto edge = get_MPO_edge_left();
    auto dim  = edge.dimension(0);
    return edge.contract(edge, Textra::idx()).reshape(Textra::array1{dim * dim});
}
Eigen::Tensor<Scalar, 1> class_ising_sdual::get_MPO2_edge_right() const {
    auto edge = get_MPO_edge_right();
    auto dim  = edge.dimension(0);
    return edge.contract(edge, Textra::idx()).reshape(Textra::array1{dim * dim});
}


void class_ising_sdual::randomize_hamiltonian() {
    if(std::string_view(h5tb.param.distribution) == "normal") {
        h5tb.param.J_rand = rnd::normal(h5tb.param.J_mean, h5tb.param.J_stdv);
        h5tb.param.h_rand = rnd::normal(h5tb.param.h_mean, h5tb.param.h_stdv);
    } else if(std::string_view(h5tb.param.distribution) == "lognormal") {
        h5tb.param.J_rand = rnd::log_normal(h5tb.param.J_mean, h5tb.param.J_stdv);
        h5tb.param.h_rand = rnd::log_normal(h5tb.param.h_mean, h5tb.param.h_stdv);
    } else if(std::string_view(h5tb.param.distribution) == "uniform") {
        h5tb.param.J_rand = rnd::uniform_double_box(h5tb.param.J_mean - h5tb.param.J_stdv / 2.0, h5tb.param.J_mean + h5tb.param.J_stdv / 2.0);
        h5tb.param.h_rand = rnd::uniform_double_box(h5tb.param.h_mean - h5tb.param.h_stdv / 2.0, h5tb.param.h_mean + h5tb.param.h_stdv / 2.0);
    } else if(std::string_view(h5tb.param.distribution) == "constant"){
        h5tb.param.J_rand = h5tb.param.J_mean;
        h5tb.param.h_rand = h5tb.param.h_mean;
    } else {
        throw std::runtime_error("Wrong distribution given. Expected one of <normal>, <lognormal>, <uniform>");
    }

    all_mpo_parameters_have_been_set = false;
    mpo_squared                      = std::nullopt;
    unique_id                        = std::nullopt;
    unique_id_sq                     = std::nullopt;

}

void class_ising_sdual::set_coupling_damping(double alpha_) {
    alpha = alpha_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_coupling() * sz);
        mpo_squared                                                                     = std::nullopt;
        unique_id                                                                       = std::nullopt;
        unique_id_sq                                                                    = std::nullopt;
    }
}
void class_ising_sdual::set_field_damping(double beta_) {
    beta = beta_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_field() * sx - e_reduced * id);
        mpo_squared                                                                     = std::nullopt;
        unique_id                                                                       = std::nullopt;
        unique_id_sq                                                                    = std::nullopt;
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
            h5tb.param.J_pert = rnd::uniform_double_box(-coupling_ptb, coupling_ptb);
            h5tb.param.h_pert = rnd::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_PERCENTAGE: {
            h5tb.param.J_pert = h5tb.param.J_rand * rnd::uniform_double_box(-coupling_ptb, coupling_ptb);
            h5tb.param.h_pert = h5tb.param.h_rand * rnd::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
    }
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_field() * sx - e_reduced * id);
        mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_coupling() * sz);
        mpo_squared                                                                     = std::nullopt;
        unique_id                                                                       = std::nullopt;
        unique_id_sq                                                                    = std::nullopt;
    }
    if(coupling_ptb == 0.0 and field_ptb == 0 and is_perturbed())
        throw std::runtime_error(fmt::format("MPO({}): Should have become unperturbed!", get_position()));
}

bool class_ising_sdual::is_perturbed() const { return h5tb.param.J_pert != 0.0 or h5tb.param.h_pert != 0.0; }


Eigen::Tensor<Scalar, 4> class_ising_sdual::MPO_nbody_view(const std::vector<size_t> &nbody_terms) const {
    // This function returns a view of the MPO including only n-body terms.
    // For instance, if nbody_terms == {2,3}, this would exclude on-site terms.
    if(nbody_terms.empty()) return MPO();
    double J1 = 0, J2 = 0;
    for(auto &&n : nbody_terms){
        if(n == 1) J1 = 1.0;
        if(n == 2) J2 = 1.0;
    }
    Eigen::Tensor<Scalar, 4> MPO_nbody = MPO();
    MPO_nbody.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-J1*get_field() * sx - e_reduced * id);
    MPO_nbody.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-J2*get_coupling() * sz);
    MPO_nbody.slice(Eigen::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-J2*(h5tb.param.lambda * h5tb.param.h_avrg) * sx);
    MPO_nbody.slice(Eigen::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-J2*(h5tb.param.lambda * h5tb.param.J_avrg) * sz);
    return MPO_nbody;
}


Eigen::Tensor<Scalar, 4> class_ising_sdual::MPO_reduced_view() const {
    if(e_reduced == 0) { return MPO(); }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> class_ising_sdual::MPO_reduced_view(double site_energy) const {
    if(site_energy == 0) { return MPO(); }
    Eigen::Tensor<Scalar, 4> temp                                           = MPO();
    temp.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixToTensor(-get_field() * sx - site_energy * id);
    return temp;
}


std::unique_ptr<class_mpo_site> class_ising_sdual::clone() const { return std::make_unique<class_ising_sdual>(*this); }

long class_ising_sdual::get_spin_dimension() const { return h5tb.param.spin_dim; }

void class_ising_sdual::set_averages(std::vector<TableMap> all_parameters, bool infinite, bool reverse) {
    if(reverse) {
        // We need to reverse the parameters, and move them one step
        std::reverse(all_parameters.begin(), all_parameters.end());
        for(size_t pos = 0; pos < all_parameters.size(); pos++) {
            all_parameters[pos]["position"] = pos;
            if(pos < all_parameters.size() - 1) {
                if(infinite){
                    all_parameters[pos]["J_rand"] = pos < all_parameters.size() - 1 ? all_parameters[pos + 1]["J_rand"] : all_parameters[0]["J_rand"];
                    all_parameters[pos]["J_pert"] = pos < all_parameters.size() - 1 ? all_parameters[pos + 1]["J_pert"] : all_parameters[0]["J_pert"];
                }else {
                    all_parameters[pos]["J_rand"] = pos < all_parameters.size() - 1 ? all_parameters[pos + 1]["J_rand"] : 0.0;
                    all_parameters[pos]["J_pert"] = pos < all_parameters.size() - 1 ? all_parameters[pos + 1]["J_pert"] : 0.0;
                }

            }
        }
    } else {
        if(not infinite){
            all_parameters.back()["J_rand"] = 0.0;
            all_parameters.back()["J_pert"] = 0.0;
        }
    }

    // Recompute J_avrg and pm.param.h_avrg from given pm.param.J_rand and pm.param.h_rand on all sites
    double J_avrg_ = 0;
    double h_avrg_ = 0;
    if(infinite){
        J_avrg_ = h5tb.param.J_mean;
        h_avrg_ = h5tb.param.h_mean;
    }else{
        for(auto &site_param : all_parameters) {
            J_avrg_ += std::any_cast<double>(site_param["J_rand"]);
            h_avrg_ += std::any_cast<double>(site_param["h_rand"]);
        }
        J_avrg_ /= static_cast<double>(all_parameters.size() - 1);
        h_avrg_ /= static_cast<double>(all_parameters.size());
    }

    for(auto &site_params : all_parameters) {
        site_params["J_avrg"] = J_avrg_;
        site_params["h_avrg"] = h_avrg_;
    }
    if(parity_sep) psfactor = (J_avrg_ + h_avrg_) * (1.0 + h5tb.param.lambda) * static_cast<double>(all_parameters.size());
    set_parameters(all_parameters[get_position()]);
}

void class_ising_sdual::save_hamiltonian(h5pp::File &file, const std::string &hamiltonian_table_path) const {
    if(not file.linkExists(hamiltonian_table_path)) file.createTable(h5tb_ising_sdual::h5_type, hamiltonian_table_path, "Selfdual Ising");
    file.appendTableRecords(h5tb, hamiltonian_table_path);
    // Position 0 is also responsible for writing attributes
    if(position.value() != 0) return;
    file.writeAttribute(h5tb.param.J_mean, "J_mean", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.J_stdv, "J_stdv", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.J_avrg, "J_avrg", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.h_mean, "h_mean", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.h_stdv, "h_stdv", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.h_avrg, "h_avrg", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.lambda, "lambda", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.delta, "delta", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.distribution, "distribution", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.spin_dim, "spin_dim", hamiltonian_table_path);
}

void class_ising_sdual::load_hamiltonian(const h5pp::File &file, const std::string &hamiltonian_table_path) {
    if(file.linkExists(hamiltonian_table_path)) {
        h5tb.param                       = file.readTableRecords<h5tb_ising_sdual::table>(hamiltonian_table_path, position);
        all_mpo_parameters_have_been_set = true;
        build_mpo();
    } else
        throw std::runtime_error(fmt::format("Could not load MPO. Table [{}] does not exist", hamiltonian_table_path));

    // Check that we are on the same point of the phase diagram
    if(std::abs(h5tb.param.delta - settings::model::ising_sdual::delta) > 1e-6)
        throw std::runtime_error(
            fmt::format("delta {:.16f} != {:.16f} settings::model::ising_sdual::delta", h5tb.param.delta, settings::model::ising_sdual::delta));
    if(std::abs(h5tb.param.lambda - settings::model::ising_sdual::lambda) > 1e-6)
        throw std::runtime_error(
            fmt::format("lambda {:.16f} != {:.16f} settings::model::ising_sdual::lambda", h5tb.param.lambda, settings::model::ising_sdual::lambda));
    if(std::abs(h5tb.param.J_stdv - settings::model::ising_sdual::J_stdv) > 1e-6)
        throw std::runtime_error(
            fmt::format("J_stdv {:.16f} != {:.16f} settings::model::ising_sdual::J_stdv", h5tb.param.J_stdv, settings::model::ising_sdual::J_stdv));
    if(std::abs(h5tb.param.h_stdv - settings::model::ising_sdual::h_stdv) > 1e-6)
        throw std::runtime_error(
            fmt::format("h_stdv {:.16f} != {:.16f} settings::model::ising_sdual::h_stdv", h5tb.param.h_stdv, settings::model::ising_sdual::h_stdv));
    if(h5tb.param.distribution != settings::model::ising_sdual::distribution)
        throw std::runtime_error(fmt::format("distribution {} != {} settings::model::ising_sdual::distribution", h5tb.param.distribution,
                                             settings::model::ising_sdual::distribution));

    // Sanity check on delta, J_mean, h_mean
    double delta_check = std::log(h5tb.param.J_mean) - std::log(h5tb.param.h_mean);
    if(std::abs(h5tb.param.delta - delta_check) > 1e-10)
        throw std::logic_error(
                fmt::format("Error when transforming delta to (J_mean, h_mean): delta {:.12f} != {:.16f} delta_check", h5tb.param.delta, delta_check));


}
