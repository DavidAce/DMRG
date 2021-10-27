#include "IsingSdual.h"
#include <config/settings.h>
#include <debug/exceptions.h>
#include <h5pp/h5pp.h>
#include <iomanip>
#include <math/num.h>
#include <math/rnd.h>
#include <math/tenx.h>
#include <qm/spin.h>

double delta_to_J_mean(double delta) { return delta > 0 ? 1.0 : std::exp(delta); }

double delta_to_h_mean(double delta) { return delta > 0 ? std::exp(-delta) : 1.0; }

IsingSdual::IsingSdual(ModelType model_type_, size_t position_) : MpoSite(model_type_, position_) {
    h5tb.param.lambda = settings::model::ising_sdual::lambda;
    h5tb.param.delta  = settings::model::ising_sdual::delta;
    h5tb.param.J_mean = delta_to_J_mean(h5tb.param.delta);
    h5tb.param.h_mean = delta_to_h_mean(h5tb.param.delta);
    h5tb.param.J_stdv = settings::model::ising_sdual::J_stdv;
    h5tb.param.h_stdv = settings::model::ising_sdual::h_stdv;

    // Sanity check on delta, J_mean, h_mean
    double delta_check = std::log(h5tb.param.J_mean) - std::log(h5tb.param.h_mean);
    if(std::abs(h5tb.param.delta - delta_check) > 1e-10)
        throw except::logic_error("error when transforming delta to (J_mean, h_mean): delta {:.12f} != {:.16f} delta_check", h5tb.param.delta, delta_check);

    h5tb.param.spin_dim = settings::model::ising_sdual::spin_dim;
    copy_c_str(settings::model::ising_sdual::distribution, h5tb.param.distribution);
    parity_sep = settings::model::ising_sdual::parity_sep;

    extent4 = {1, 1, h5tb.param.spin_dim, h5tb.param.spin_dim};
    extent2 = {h5tb.param.spin_dim, h5tb.param.spin_dim};
    h5tb_ising_sdual::register_table_type();
}

double IsingSdual::get_coupling() const { return h5tb.param.J_rand + h5tb.param.J_pert; }
double IsingSdual::get_field() const { return h5tb.param.h_rand + h5tb.param.h_pert; }
void   IsingSdual::print_parameter_names() const { h5tb.print_parameter_names(); }
void   IsingSdual::print_parameter_values() const { h5tb.print_parameter_values(); }

void IsingSdual::set_parameters(TableMap &parameters) {
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
    copy_c_str(std::any_cast<std::string>(parameters["distribution"]), h5tb.param.distribution);
    all_mpo_parameters_have_been_set = true;
}

IsingSdual::TableMap IsingSdual::get_parameters() const {
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

void IsingSdual::build_mpo()
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
    using namespace qm::spin::half;
    tools::log->debug("mpo({}): building sdual mpo", get_position());
    if(not all_mpo_parameters_have_been_set)
        throw except::runtime_error("mpo({}): can't build mpo: full lattice parameters haven't been set yet.", get_position());
    if(parity_sep) {
        mpo_internal.resize(6, 6, h5tb.param.spin_dim, h5tb.param.spin_dim);
        mpo_internal.setZero();
        mpo_internal.slice(std::array<long, 4>{5, 5, 0, 0}, extent4).reshape(extent2) =
            tenx::TensorMap(sx); // Multiply the psfactor on the edge! Not on each MPO!
    } else {
        mpo_internal.resize(5, 5, h5tb.param.spin_dim, h5tb.param.spin_dim);
        mpo_internal.setZero();
    }

    mpo_internal.slice(std::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(id);
    mpo_internal.slice(std::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(sz);
    mpo_internal.slice(std::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(sx);
    mpo_internal.slice(std::array<long, 4>{3, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(id);
    mpo_internal.slice(std::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-get_field() * sx - e_reduced * id);
    mpo_internal.slice(std::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-get_coupling() * sz);
    mpo_internal.slice(std::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-(h5tb.param.lambda * h5tb.param.h_avrg) * sx);
    mpo_internal.slice(std::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-(h5tb.param.lambda * h5tb.param.J_avrg) * sz);
    mpo_internal.slice(std::array<long, 4>{4, 4, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(id);
    if(tenx::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw except::runtime_error("mpo({}): found nan", get_position());
    }
    unique_id = std::nullopt;
}

void IsingSdual::randomize_hamiltonian() {
    if(std::string_view(h5tb.param.distribution) == "normal") {
        h5tb.param.J_rand = rnd::normal(h5tb.param.J_mean, h5tb.param.J_stdv);
        h5tb.param.h_rand = rnd::normal(h5tb.param.h_mean, h5tb.param.h_stdv);
    } else if(std::string_view(h5tb.param.distribution) == "lognormal") {
        h5tb.param.J_rand = rnd::log_normal(h5tb.param.J_mean, h5tb.param.J_stdv);
        h5tb.param.h_rand = rnd::log_normal(h5tb.param.h_mean, h5tb.param.h_stdv);
    } else if(std::string_view(h5tb.param.distribution) == "uniform") {
        h5tb.param.J_rand = rnd::uniform_double_box(h5tb.param.J_mean - h5tb.param.J_stdv / 2.0, h5tb.param.J_mean + h5tb.param.J_stdv / 2.0);
        h5tb.param.h_rand = rnd::uniform_double_box(h5tb.param.h_mean - h5tb.param.h_stdv / 2.0, h5tb.param.h_mean + h5tb.param.h_stdv / 2.0);
    } else if(std::string_view(h5tb.param.distribution) == "constant") {
        h5tb.param.J_rand = h5tb.param.J_mean;
        h5tb.param.h_rand = h5tb.param.h_mean;
    } else {
        throw except::runtime_error("wrong distribution [{}]: expected one of normal | lognormal | uniform | constant", h5tb.param.distribution);
    }

    all_mpo_parameters_have_been_set = false;
    mpo_squared                      = std::nullopt;
    unique_id                        = std::nullopt;
    unique_id_sq                     = std::nullopt;
}

void IsingSdual::set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) {
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
        using namespace qm::spin::half;
        mpo_internal.slice(std::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-get_field() * sx - e_reduced * id);
        mpo_internal.slice(std::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-get_coupling() * sz);
        mpo_squared                                                                   = std::nullopt;
        unique_id                                                                     = std::nullopt;
        unique_id_sq                                                                  = std::nullopt;
    }
    if(coupling_ptb == 0.0 and field_ptb == 0 and is_perturbed()) throw except::runtime_error("mpo({}): should have become unperturbed!", get_position());
}

bool IsingSdual::is_perturbed() const { return h5tb.param.J_pert != 0.0 or h5tb.param.h_pert != 0.0; }

Eigen::Tensor<MpoSite::cplx, 4> IsingSdual::MPO_nbody_view(std::optional<std::vector<size_t>>                  nbody,
                                                           [[maybe_unused]] std::optional<std::vector<size_t>> skip) const {
    // This function returns a view of the MPO including only n-body terms.
    // For instance, if nbody_terms == {2,3}, this would exclude on-site terms.
    // Next-nearest neighbor terms are counted as 3-body terms because 3 sites are involved: the skipped site counts

    if(not nbody) return MPO();
    double J1 = 0, J2 = 0, J3 = 0;
    for(const auto &n : nbody.value()) {
        if(n == 1) J1 = 1.0;
        if(n == 2) J2 = 1.0;
        if(n == 3) J3 = 1.0;
    }
    using namespace qm::spin::half;
    Eigen::Tensor<cplx, 4> MPO_nbody                                           = MPO();
    MPO_nbody.slice(std::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-J1 * get_field() * sx - e_reduced * id);
    MPO_nbody.slice(std::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-J2 * get_coupling() * sz);
    MPO_nbody.slice(std::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-J2 * (h5tb.param.lambda * h5tb.param.h_avrg) * sx);
    MPO_nbody.slice(std::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-J3 * (h5tb.param.lambda * h5tb.param.J_avrg) * sz);
    return MPO_nbody;
}

Eigen::Tensor<MpoSite::cplx, 4> IsingSdual::MPO_reduced_view() const {
    if(e_reduced == 0) { return MPO(); }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<MpoSite::cplx, 4> IsingSdual::MPO_reduced_view(double site_energy) const {
    using namespace qm::spin::half;
    Eigen::Tensor<cplx, 4> temp                                               = MPO();
    long                   row                                                = temp.dimension(0) - 1;
    long                   col                                                = 0;
    temp.slice(std::array<long, 4>{row, col, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-get_field() * sx - site_energy * id);
    return temp;
}

std::unique_ptr<MpoSite> IsingSdual::clone() const { return std::make_unique<IsingSdual>(*this); }

long IsingSdual::get_spin_dimension() const { return h5tb.param.spin_dim; }

void IsingSdual::set_averages(std::vector<TableMap> all_parameters, bool infinite, bool reverse) {
    if(reverse) {
        // We need to reverse the parameters, and move them one step
        std::reverse(all_parameters.begin(), all_parameters.end());
        for(size_t pos = 0; pos < all_parameters.size(); pos++) {
            all_parameters[pos]["position"] = pos;
            if(pos < all_parameters.size() - 1) {
                if(infinite) {
                    all_parameters[pos]["J_rand"] = pos < all_parameters.size() - 1 ? all_parameters[pos + 1]["J_rand"] : all_parameters[0]["J_rand"];
                    all_parameters[pos]["J_pert"] = pos < all_parameters.size() - 1 ? all_parameters[pos + 1]["J_pert"] : all_parameters[0]["J_pert"];
                } else {
                    all_parameters[pos]["J_rand"] = pos < all_parameters.size() - 1 ? all_parameters[pos + 1]["J_rand"] : 0.0;
                    all_parameters[pos]["J_pert"] = pos < all_parameters.size() - 1 ? all_parameters[pos + 1]["J_pert"] : 0.0;
                }
            }
        }
    } else {
        if(not infinite) {
            all_parameters.back()["J_rand"] = 0.0;
            all_parameters.back()["J_pert"] = 0.0;
        }
    }

    // Recompute J_avrg and pm.param.h_avrg from given pm.param.J_rand and pm.param.h_rand on all sites
    double J_avrg_ = 0;
    double h_avrg_ = 0;
    if(infinite) {
        J_avrg_ = h5tb.param.J_mean;
        h_avrg_ = h5tb.param.h_mean;
    } else {
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

void IsingSdual::save_hamiltonian(h5pp::File &file, std::string_view hamiltonian_table_path) const {
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

void IsingSdual::load_hamiltonian(const h5pp::File &file, std::string_view model_path) {
    auto ham_table = fmt::format("{}/hamiltonian", model_path);
    if(file.linkExists(ham_table)) {
        h5tb.param                       = file.readTableRecords<h5tb_ising_sdual::table>(ham_table, position);
        all_mpo_parameters_have_been_set = true;
    } else throw except::runtime_error("could not load mpo: table [{}] does not exist", ham_table);

    // Check that we are on the same point of the phase diagram
    using namespace settings::model::ising_sdual;
    if(std::abs(h5tb.param.delta - delta) > 1e-6) throw except::runtime_error("delta  {:.16f} != {:.16f} ising_sdual::delta", h5tb.param.delta, delta);
    if(std::abs(h5tb.param.lambda - lambda) > 1e-6) throw except::runtime_error("lambda {:.16f} != {:.16f} ising_sdual::lambda", h5tb.param.lambda, lambda);
    if(std::abs(h5tb.param.J_stdv - J_stdv) > 1e-6) throw except::runtime_error("J_stdv {:.16f} != {:.16f} ising_sdual::J_stdv", h5tb.param.J_stdv, J_stdv);
    if(std::abs(h5tb.param.h_stdv - h_stdv) > 1e-6) throw except::runtime_error("h_stdv {:.16f} != {:.16f} ising_sdual::h_stdv", h5tb.param.h_stdv, h_stdv);
    if(h5tb.param.distribution != distribution)
        throw except::runtime_error("distribution {} != {} ising_sdual::distribution", h5tb.param.distribution, distribution);

    // Sanity check on delta, J_mean, h_mean
    double delta_check = std::log(h5tb.param.J_mean) - std::log(h5tb.param.h_mean);
    if(std::abs(h5tb.param.delta - delta_check) > 1e-10)
        throw except::logic_error("Error when transforming delta to (J_mean, h_mean): delta {:.12f} != {:.16f} delta_check", h5tb.param.delta, delta_check);
}
