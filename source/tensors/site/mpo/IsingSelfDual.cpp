#include "IsingSelfDual.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "math/num.h"
#include "math/rnd.h"
#include "math/tenx.h"
#include "qm/spin.h"
#include <h5pp/h5pp.h>
#include <iomanip>

double delta_to_J_mean(double delta) { return std::min(1.0, std::exp(delta)); }

double delta_to_h_mean(double delta) { return std::min(1.0, std::exp(-delta)); }

IsingSelfDual::IsingSelfDual(ModelType model_type_, size_t position_) : MpoSite(model_type_, position_) {
    h5tb.param.lambda = settings::model::ising_sdual::lambda;
    h5tb.param.delta  = settings::model::ising_sdual::delta;
    h5tb.param.J_mean = delta_to_J_mean(h5tb.param.delta);
    h5tb.param.h_mean = delta_to_h_mean(h5tb.param.delta);
    // Sanity check on delta, J_mean, h_mean
    double delta_check = std::log(h5tb.param.J_mean) - std::log(h5tb.param.h_mean);
    if(std::abs(h5tb.param.delta - delta_check) > 1e-10)
        throw except::logic_error("error when transforming delta to (J_mean, h_mean): delta {:.12f} != {:.16f} delta_check", h5tb.param.delta, delta_check);

    h5tb.param.spin_dim = settings::model::ising_sdual::spin_dim;
    copy_c_str(settings::model::ising_sdual::distribution, h5tb.param.distribution);
    parity_sep = settings::model::ising_sdual::parity_sep;

    extent4 = {1, 1, h5tb.param.spin_dim, h5tb.param.spin_dim};
    extent2 = {h5tb.param.spin_dim, h5tb.param.spin_dim};
    h5tb_ising_selfdual::register_table_type();
}

double IsingSelfDual::get_coupling() const { return h5tb.param.J_rand; }
double IsingSelfDual::get_field() const { return h5tb.param.h_rand; }
void   IsingSelfDual::print_parameter_names() const { h5tb.print_parameter_names(); }
void   IsingSelfDual::print_parameter_values() const { h5tb.print_parameter_values(); }

void IsingSelfDual::set_parameters(TableMap &parameters) {
    h5tb.param.J_mean   = std::any_cast<double>(parameters["J_mean"]);
    h5tb.param.J_wdth   = std::any_cast<double>(parameters["J_wdth"]);
    h5tb.param.J_rand   = std::any_cast<double>(parameters["J_rand"]);
    h5tb.param.h_mean   = std::any_cast<double>(parameters["h_mean"]);
    h5tb.param.h_wdth   = std::any_cast<double>(parameters["h_wdth"]);
    h5tb.param.h_rand   = std::any_cast<double>(parameters["h_rand"]);
    h5tb.param.lambda   = std::any_cast<double>(parameters["lambda"]);
    h5tb.param.delta    = std::any_cast<double>(parameters["delta"]);
    h5tb.param.spin_dim = std::any_cast<long>(parameters["spin_dim"]);
    copy_c_str(std::any_cast<std::string>(parameters["distribution"]), h5tb.param.distribution);
    all_mpo_parameters_have_been_set = true;
}

IsingSelfDual::TableMap IsingSelfDual::get_parameters() const {
    /* clang-format off */
    TableMap parameters;
    parameters["J_mean"]   = h5tb.param.J_mean;
    parameters["J_wdth"]   = h5tb.param.J_wdth;
    parameters["J_rand"]   = h5tb.param.J_rand;
    parameters["h_mean"]   = h5tb.param.h_mean;
    parameters["h_wdth"]   = h5tb.param.h_wdth;
    parameters["h_rand"]   = h5tb.param.h_rand;
    parameters["lambda"]   = h5tb.param.lambda;
    parameters["delta"]    = h5tb.param.delta;
    parameters["spin_dim"] = h5tb.param.spin_dim;
    parameters["distribution"] = std::string(h5tb.param.distribution);
    return parameters;
    /* clang-format on */
}

void IsingSelfDual::build_mpo()
/*! Builds the MPO hamiltonian as a rank 4 tensor. Notation following Schollwöck (2010)
 *
 * H = Σ J_{i} σx_{i} σx_{i+1} +  h_{i} σz_{i} + λ (h_mean σz_i σz_{i+1} + J_mean σx_{i} σx_{i+2})
 *
 *  |     I            0          0           0           0   |
 *  |     σx           0          0           0           0   |
 *  |     σz           0          0           0           0   |
 *  |     0            I          0           0           0   |
 *  |   h_rand*σz  J_rand*σx  l*h_mean*σz  l*J_mean*σx    I   |
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
    tools::log->debug("mpo({}): building ising-selfdual mpo", get_position());
    if(not all_mpo_parameters_have_been_set)
        throw except::runtime_error("mpo({}): can't build mpo: full lattice parameters haven't been set yet.", get_position());

    mpo_internal.resize(5, 5, h5tb.param.spin_dim, h5tb.param.spin_dim);
    mpo_internal.setZero();
    if(parity_sep) {
        mpo_internal.resize(6, 6, h5tb.param.spin_dim, h5tb.param.spin_dim);
        mpo_internal.setZero();
        // Multiply the psfactor on the edge! Not on each MPO!
        mpo_internal.slice(std::array<long, 4>{5, 5, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(sz);
    }

    mpo_internal.slice(std::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(id);
    mpo_internal.slice(std::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(sx);
    mpo_internal.slice(std::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(sz);
    mpo_internal.slice(std::array<long, 4>{3, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(id);
    mpo_internal.slice(std::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(h5tb.param.h_rand * sz - e_shift * id);
    mpo_internal.slice(std::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(h5tb.param.J_rand * sx);
    mpo_internal.slice(std::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(h5tb.param.h_mean * h5tb.param.lambda * sz);
    mpo_internal.slice(std::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(h5tb.param.J_mean * h5tb.param.lambda * sx);
    mpo_internal.slice(std::array<long, 4>{4, 4, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(id);
    if(tenx::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw except::runtime_error("mpo({}): found nan", get_position());
    }
    unique_id = std::nullopt;
}

void IsingSelfDual::randomize_hamiltonian() {
    if(std::string_view(h5tb.param.distribution) == "normal") {
        h5tb.param.J_wdth = 1.0;
        h5tb.param.J_wdth = 1.0;
        h5tb.param.J_rand = rnd::normal(h5tb.param.J_mean, h5tb.param.J_wdth);
        h5tb.param.h_rand = rnd::normal(h5tb.param.h_mean, h5tb.param.J_wdth);
    } else if(std::string_view(h5tb.param.distribution) == "lognormal") {
        // See http://arxiv.org/abs/1711.00020
        h5tb.param.J_wdth = 1.0;
        h5tb.param.J_wdth = 1.0;
        h5tb.param.J_rand = rnd::log_normal(h5tb.param.J_mean, h5tb.param.J_wdth);
        h5tb.param.h_rand = rnd::log_normal(h5tb.param.h_mean, h5tb.param.J_wdth);
    } else if(std::string_view(h5tb.param.distribution) == "uniform") {
        h5tb.param.J_wdth = 2.0 * h5tb.param.J_mean;
        h5tb.param.J_wdth = 2.0 * h5tb.param.h_mean;
        h5tb.param.J_rand = rnd::uniform_double_box(0, h5tb.param.J_wdth);
        h5tb.param.h_rand = rnd::uniform_double_box(0, h5tb.param.h_wdth);
    } else if(std::string_view(h5tb.param.distribution) == "constant") {
        h5tb.param.J_wdth = 0.0;
        h5tb.param.J_wdth = 0.0;
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

Eigen::Tensor<MpoSite::cplx, 4> IsingSelfDual::MPO_nbody_view(std::optional<std::vector<size_t>>                  nbody,
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
    MPO_nbody.slice(std::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(J1 * h5tb.param.h_rand * sz - e_shift * id);
    MPO_nbody.slice(std::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(J2 * h5tb.param.J_rand * sx);
    MPO_nbody.slice(std::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(J2 * h5tb.param.h_mean * h5tb.param.lambda * sz);
    MPO_nbody.slice(std::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(J3 * h5tb.param.J_mean * h5tb.param.lambda * sx);
    return MPO_nbody;
}

Eigen::Tensor<MpoSite::cplx, 4> IsingSelfDual::MPO_shifted_view() const { return MPO_shifted_view(e_shift); }

Eigen::Tensor<MpoSite::cplx, 4> IsingSelfDual::MPO_shifted_view(double energy_shift_per_site) const {
    using namespace qm::spin::half;
    Eigen::Tensor<cplx, 4> temp                                           = MPO();
    temp.slice(std::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(get_field() * sz - energy_shift_per_site * id);
    return temp;
}

std::unique_ptr<MpoSite> IsingSelfDual::clone() const { return std::make_unique<IsingSelfDual>(*this); }

long IsingSelfDual::get_spin_dimension() const { return h5tb.param.spin_dim; }

void IsingSelfDual::set_averages(std::vector<TableMap> all_parameters, bool infinite) {
    if(not infinite) { all_parameters.back()["J_rand"] = 0.0; }
    if(parity_sep) psfactor = 2.0 * (h5tb.param.J_mean + h5tb.param.h_mean) * (1.0 + h5tb.param.lambda) * static_cast<double>(all_parameters.size());
    set_parameters(all_parameters[get_position()]);
}

void IsingSelfDual::save_hamiltonian(h5pp::File &file, std::string_view hamiltonian_table_path) const {
    if(not file.linkExists(hamiltonian_table_path)) file.createTable(h5tb_ising_selfdual::h5_type, hamiltonian_table_path, "Selfdual Ising");
    file.appendTableRecords(h5tb, hamiltonian_table_path);
    // Position 0 is also responsible for writing attributes
    if(position.value() != 0) return;
    file.writeAttribute(h5tb.param.J_mean, "J_mean", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.h_mean, "h_mean", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.lambda, "lambda", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.delta, "delta", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.distribution, "distribution", hamiltonian_table_path);
    file.writeAttribute(h5tb.param.spin_dim, "spin_dim", hamiltonian_table_path);
}

void IsingSelfDual::load_hamiltonian(const h5pp::File &file, std::string_view model_path) {
    auto ham_table = fmt::format("{}/hamiltonian", model_path);
    if(file.linkExists(ham_table)) {
        h5tb.param                       = file.readTableRecords<h5tb_ising_selfdual::table>(ham_table, position);
        all_mpo_parameters_have_been_set = true;
    } else
        throw except::runtime_error("could not load mpo: table [{}] does not exist", ham_table);

    // Check that we are on the same point of the phase diagram
    using namespace settings::model::ising_sdual;
    if(std::abs(h5tb.param.delta - delta) > 1e-6) throw except::runtime_error("delta  {:.16f} != {:.16f} ising_sdual::delta", h5tb.param.delta, delta);
    if(std::abs(h5tb.param.lambda - lambda) > 1e-6) throw except::runtime_error("lambda {:.16f} != {:.16f} ising_sdual::lambda", h5tb.param.lambda, lambda);
    if(h5tb.param.distribution != distribution)
        throw except::runtime_error("distribution {} != {} ising_sdual::distribution", h5tb.param.distribution, distribution);

    // Sanity check on delta, J_mean, h_mean
    double delta_check = std::log(h5tb.param.J_mean) - std::log(h5tb.param.h_mean);
    if(std::abs(h5tb.param.delta - delta_check) > 1e-10)
        throw except::logic_error("Error when transforming delta to (J_mean, h_mean): delta {:.12f} != {:.16f} delta_check", h5tb.param.delta, delta_check);
}
