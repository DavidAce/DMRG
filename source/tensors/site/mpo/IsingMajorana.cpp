#include "IsingMajorana.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "math/num.h"
#include "math/rnd.h"
#include "math/tenx.h"
#include "qm/spin.h"
#include <h5pp/h5pp.h>
#include <iomanip>

double delta_to_J_wdth(double delta) { return std::exp(delta / 2.0); }

double delta_to_h_wdth(double delta) { return std::exp(-delta / 2.0); }

IsingMajorana::IsingMajorana(ModelType model_type_, size_t position_) : MpoSite(model_type_, position_) {
    h5tb.param.g      = settings::model::ising_majorana::g;
    h5tb.param.delta  = settings::model::ising_majorana::delta;
    h5tb.param.J_wdth = delta_to_J_wdth(h5tb.param.delta);
    h5tb.param.h_wdth = delta_to_h_wdth(h5tb.param.delta);
    h5tb.param.J_mean = 0.5 * h5tb.param.J_wdth;
    h5tb.param.h_mean = 0.5 * h5tb.param.h_wdth;

    // Sanity check on delta, J_mean, h_mean
    double delta_check = std::log(h5tb.param.J_mean) - std::log(h5tb.param.h_mean);
    if(std::abs(h5tb.param.delta - delta_check) > 1e-10)
        throw except::logic_error("error when transforming delta to (J_mean, h_mean): delta {:.12f} != {:.16f} delta_check", h5tb.param.delta, delta_check);

    h5tb.param.spin_dim = settings::model::ising_majorana::spin_dim;
    copy_c_str(settings::model::ising_majorana::distribution, h5tb.param.distribution);
    parity_sep = settings::model::ising_majorana::parity_sep;

    extent4 = {1, 1, h5tb.param.spin_dim, h5tb.param.spin_dim};
    extent2 = {h5tb.param.spin_dim, h5tb.param.spin_dim};
    h5tb_ising_majorana::register_table_type();
}

double IsingMajorana::get_coupling() const { return h5tb.param.J_rand; }
double IsingMajorana::get_field() const { return h5tb.param.h_rand; }
void   IsingMajorana::print_parameter_names() const { h5tb.print_parameter_names(); }
void   IsingMajorana::print_parameter_values() const { h5tb.print_parameter_values(); }

void IsingMajorana::set_parameters(TableMap &parameters) {
    h5tb.param.J_mean   = std::any_cast<double>(parameters["J_mean"]);
    h5tb.param.J_wdth   = std::any_cast<double>(parameters["J_wdth"]);
    h5tb.param.J_rand   = std::any_cast<double>(parameters["J_rand"]);
    h5tb.param.h_mean   = std::any_cast<double>(parameters["h_mean"]);
    h5tb.param.h_wdth   = std::any_cast<double>(parameters["h_wdth"]);
    h5tb.param.h_rand   = std::any_cast<double>(parameters["h_rand"]);
    h5tb.param.g        = std::any_cast<double>(parameters["g"]);
    h5tb.param.delta    = std::any_cast<double>(parameters["delta"]);
    h5tb.param.spin_dim = std::any_cast<long>(parameters["spin_dim"]);
    copy_c_str(std::any_cast<std::string>(parameters["distribution"]), h5tb.param.distribution);
    all_mpo_parameters_have_been_set = true;
}

IsingMajorana::TableMap IsingMajorana::get_parameters() const {
    /* clang-format off */
    TableMap parameters;
    parameters["J_mean"]   = h5tb.param.J_mean;
    parameters["J_wdth"]   = h5tb.param.J_wdth;
    parameters["J_rand"]   = h5tb.param.J_rand;
    parameters["h_mean"]   = h5tb.param.h_mean;
    parameters["h_wdth"]   = h5tb.param.h_wdth;
    parameters["h_rand"]   = h5tb.param.h_rand;
    parameters["g"]        = h5tb.param.g;
    parameters["delta"]    = h5tb.param.delta;
    parameters["spin_dim"] = h5tb.param.spin_dim;
    parameters["distribution"] = std::string(h5tb.param.distribution);
    return parameters;
    /* clang-format on */
}

void IsingMajorana::build_mpo()
/*! Builds the MPO hamiltonian as a rank 4 tensor. Notation following Schollwöck (2010)
 *
 * H = Σ J_{i} σx_{i} σx_{i+1} + h_{i} σz_{i} + g*(σz_i σz_{i+1} + σx_{i} σx_{i+2})
 *
 *  |    I           0          0          0         0   |
 *  |    σx          0          0          0         0   |
 *  |    σz          0          0          0         0   |
 *  |    0           I          0          0         0   |
 *  | h_rand*σz  J_rand*σx     g*σz      g*σx        I   |
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
    tools::log->debug("mpo({}): building ising-majorana mpo", get_position());
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
    mpo_internal.slice(std::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(get_field() * sz - e_shift * id);
    mpo_internal.slice(std::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(get_coupling() * sx);
    mpo_internal.slice(std::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(h5tb.param.g * sz);
    mpo_internal.slice(std::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(h5tb.param.g * sx);
    mpo_internal.slice(std::array<long, 4>{4, 4, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(id);
    if(tenx::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw except::runtime_error("mpo({}): found nan", get_position());
    }
    unique_id = std::nullopt;
}

void IsingMajorana::randomize_hamiltonian() {
    if(std::string_view(h5tb.param.distribution) == "normal") {
        h5tb.param.J_rand = rnd::normal(h5tb.param.J_mean, h5tb.param.J_wdth);
        h5tb.param.h_rand = rnd::normal(h5tb.param.h_mean, h5tb.param.h_wdth);
    } else if(std::string_view(h5tb.param.distribution) == "lognormal") {
        h5tb.param.J_rand = rnd::log_normal(h5tb.param.J_mean, h5tb.param.J_wdth);
        h5tb.param.h_rand = rnd::log_normal(h5tb.param.h_mean, h5tb.param.h_wdth);
    } else if(std::string_view(h5tb.param.distribution) == "uniform") {
        h5tb.param.J_rand = rnd::uniform_double_box(0, h5tb.param.J_wdth);
        h5tb.param.h_rand = rnd::uniform_double_box(0, h5tb.param.h_wdth);
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

Eigen::Tensor<MpoSite::cplx, 4> IsingMajorana::MPO_nbody_view(std::optional<std::vector<size_t>>                  nbody,
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
    MPO_nbody.slice(std::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(J1 * get_field() * sz - e_shift * id);
    MPO_nbody.slice(std::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(J2 * get_coupling() * sx);
    MPO_nbody.slice(std::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(J2 * h5tb.param.g * sz);
    MPO_nbody.slice(std::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(J3 * h5tb.param.g * sx);
    return MPO_nbody;
}

Eigen::Tensor<MpoSite::cplx, 4> IsingMajorana::MPO_shifted_view() const { return MPO_shifted_view(e_shift); }

Eigen::Tensor<MpoSite::cplx, 4> IsingMajorana::MPO_shifted_view(double energy_shift_per_site) const {
    using namespace qm::spin::half;
    Eigen::Tensor<cplx, 4> temp                                           = MPO();
    temp.slice(std::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(get_field() * sz - energy_shift_per_site * id);
    return temp;
}

std::unique_ptr<MpoSite> IsingMajorana::clone() const { return std::make_unique<IsingMajorana>(*this); }

long IsingMajorana::get_spin_dimension() const { return h5tb.param.spin_dim; }

void IsingMajorana::set_averages(std::vector<TableMap> all_parameters, bool infinite) {
    if(not infinite) { all_parameters.back()["J_rand"] = 0.0; }
    //    if(parity_sep) psfactor = 0.5;
    if(parity_sep) psfactor = 2.0 * (h5tb.param.J_mean + h5tb.param.h_mean) * (1.0 + h5tb.param.g) * static_cast<double>(all_parameters.size());
    set_parameters(all_parameters[get_position()]);
}

void IsingMajorana::save_hamiltonian(h5pp::File &file, std::string_view hamiltonian_table_path) const {
    if(not file.linkExists(hamiltonian_table_path)) file.createTable(h5tb_ising_majorana::h5_type, hamiltonian_table_path, "Selfdual Ising");
    file.appendTableRecords(h5tb, hamiltonian_table_path);
    // Position 0 is also responsible for writing attributes
    if(position.value() != 0) return;
    file.writeAttribute(h5tb.param.J_mean, hamiltonian_table_path, "J_mean");
    file.writeAttribute(h5tb.param.J_wdth, hamiltonian_table_path, "J_wdth");
    file.writeAttribute(h5tb.param.h_mean, hamiltonian_table_path, "h_mean");
    file.writeAttribute(h5tb.param.h_wdth, hamiltonian_table_path, "h_wdth");
    file.writeAttribute(h5tb.param.g, hamiltonian_table_path, "g");
    file.writeAttribute(h5tb.param.delta, hamiltonian_table_path, "delta");
    file.writeAttribute(h5tb.param.distribution, hamiltonian_table_path, "distribution");
    file.writeAttribute(h5tb.param.spin_dim, hamiltonian_table_path, "spin_dim");
}

void IsingMajorana::load_hamiltonian(const h5pp::File &file, std::string_view model_path) {
    auto ham_table = fmt::format("{}/hamiltonian", model_path);
    if(file.linkExists(ham_table)) {
        h5tb.param                       = file.readTableRecords<h5tb_ising_majorana::table>(ham_table, position);
        all_mpo_parameters_have_been_set = true;
    } else
        throw except::runtime_error("could not load mpo: table [{}] does not exist", ham_table);

    // Check that we are on the same point of the phase diagram
    auto J_wdth = delta_to_J_wdth(h5tb.param.delta);
    auto h_wdth = delta_to_h_wdth(h5tb.param.delta);

    using namespace settings::model::ising_majorana;
    if(std::abs(h5tb.param.delta - delta) > 1e-6) throw except::runtime_error("delta  {:.16f} != {:.16f} ising_majorana::delta", h5tb.param.delta, delta);
    if(std::abs(h5tb.param.g - g) > 1e-6) throw except::runtime_error("g {:.16f} != {:.16f} ising_majorana::g", h5tb.param.g, g);
    if(std::abs(h5tb.param.J_wdth - J_wdth) > 1e-6) throw except::runtime_error("J_wdth {:.16f} != {:.16f} ising_majorana::J_wdth", h5tb.param.J_wdth, J_wdth);
    if(std::abs(h5tb.param.h_wdth - h_wdth) > 1e-6) throw except::runtime_error("h_wdth {:.16f} != {:.16f} ising_majorana::h_wdth", h5tb.param.h_wdth, h_wdth);
    if(h5tb.param.distribution != distribution)
        throw except::runtime_error("distribution {} != {} ising_majorana::distribution", h5tb.param.distribution, distribution);

    // Sanity check on delta, J_mean, h_mean
    double delta_check = std::log(h5tb.param.J_mean) - std::log(h5tb.param.h_mean);
    if(std::abs(h5tb.param.delta - delta_check) > 1e-10)
        throw except::logic_error("Error when transforming delta to (J_mean, h_mean): delta {:.12f} != {:.16f} delta_check", h5tb.param.delta, delta_check);
}
