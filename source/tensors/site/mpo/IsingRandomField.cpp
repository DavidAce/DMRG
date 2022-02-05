#include "IsingRandomField.h"
#include <config/settings.h>
#include <debug/exceptions.h>
#include <h5pp/h5pp.h>
#include <iomanip>
#include <math/num.h>
#include <math/rnd.h>
#include <math/tenx.h>
#include <qm/spin.h>

IsingRandomField::IsingRandomField(ModelType model_type_, size_t position_) : MpoSite(model_type_, position_) {
    h5tb.param.J1       = settings::model::ising_tf_rf::J1;
    h5tb.param.J2       = settings::model::ising_tf_rf::J2;
    h5tb.param.h_tran   = settings::model::ising_tf_rf::h_tran;
    h5tb.param.h_mean   = settings::model::ising_tf_rf::h_mean;
    h5tb.param.h_wdth   = settings::model::ising_tf_rf::h_wdth;
    h5tb.param.spin_dim = settings::model::ising_tf_rf::spin_dim;
    copy_c_str(settings::model::ising_tf_rf::distribution, h5tb.param.distribution);
    extent4 = {1, 1, h5tb.param.spin_dim, h5tb.param.spin_dim};
    extent2 = {h5tb.param.spin_dim, h5tb.param.spin_dim};
    using namespace qm::spin::half;
    qm::spin::half::SX = qm::spin::gen_manybody_spins(sx, 2);
    qm::spin::half::SY = qm::spin::gen_manybody_spins(sy, 2);
    qm::spin::half::SZ = qm::spin::gen_manybody_spins(sz, 2);
    qm::spin::half::II = qm::spin::gen_manybody_spins(id, 2);

    h5tb_ising_tf_rf::register_table_type();
    all_mpo_parameters_have_been_set =
        false; // There are no full lattice parameters but we set it to true here since the model is not supposed to be randomized per site
}

double IsingRandomField::get_field() const { return h5tb.param.h_tran + h5tb.param.h_rand; }
double IsingRandomField::get_coupling() const { return h5tb.param.J1; }
void   IsingRandomField::print_parameter_names() const { h5tb.print_parameter_names(); }
void   IsingRandomField::print_parameter_values() const { h5tb.print_parameter_values(); }

void IsingRandomField::set_parameters(TableMap &parameters) {
    h5tb.param.J1       = std::any_cast<double>(parameters["J1"]);
    h5tb.param.J2       = std::any_cast<double>(parameters["J2"]);
    h5tb.param.h_tran   = std::any_cast<double>(parameters["h_tran"]);
    h5tb.param.h_mean   = std::any_cast<double>(parameters["h_mean"]);
    h5tb.param.h_wdth   = std::any_cast<double>(parameters["h_wdth"]);
    h5tb.param.h_rand   = std::any_cast<double>(parameters["h_rand"]);
    h5tb.param.spin_dim = std::any_cast<long>(parameters["spin_dim"]);
    copy_c_str(std::any_cast<std::string>(parameters["distribution"]), h5tb.param.distribution);
    if(h5tb.param.J2 != 0.0) throw except::runtime_error("mpo({}): use of [J2] - Next-nearest neighbor coupling - is not implemented yet", get_position());
    all_mpo_parameters_have_been_set = true;
}

IsingRandomField::TableMap IsingRandomField::get_parameters() const {
    /* clang-format off */
    TableMap parameters;
    parameters["J1"]            = h5tb.param.J1;
    parameters["J2"]            = h5tb.param.J2;
    parameters["h_tran"]        = h5tb.param.h_tran;
    parameters["h_mean"]        = h5tb.param.h_mean;
    parameters["h_wdth"]        = h5tb.param.h_wdth;
    parameters["h_rand"]        = h5tb.param.h_rand;
    parameters["spin_dim"]      = h5tb.param.spin_dim;
    parameters["distribution"]  = std::string(h5tb.param.distribution);
    return parameters;
    /* clang-format on */
}

void IsingRandomField::build_mpo()
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
    using namespace qm::spin::half;
    tools::log->debug("mpo({}): building tf-rf ising mpo", get_position());
    if(not all_mpo_parameters_have_been_set)
        throw except::runtime_error("mpo({}): can't build mpo: full lattice parameters haven't been set yet.", get_position());
    mpo_internal.resize(3, 3, h5tb.param.spin_dim, h5tb.param.spin_dim);
    mpo_internal.setZero();
    mpo_internal.slice(std::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(id);
    mpo_internal.slice(std::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(sz);
    mpo_internal.slice(std::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-get_field() * sx - e_shift * id);
    mpo_internal.slice(std::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-get_coupling() * sz);
    mpo_internal.slice(std::array<long, 4>{2, 2, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(id);
    if(tenx::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw except::runtime_error("mpo({}): found nan", get_position());
    }
    unique_id = std::nullopt;
}

void IsingRandomField::randomize_hamiltonian() {
    if(std::string(h5tb.param.distribution) == "normal") {
        h5tb.param.h_rand = rnd::normal(h5tb.param.h_mean, h5tb.param.h_wdth);
    } else if(std::string(h5tb.param.distribution) == "lognormal") {
        h5tb.param.h_rand = rnd::log_normal(h5tb.param.h_mean, h5tb.param.h_wdth);
    } else if(std::string(h5tb.param.distribution) == "uniform") {
        h5tb.param.h_rand = rnd::uniform_double_box(h5tb.param.h_mean - h5tb.param.h_wdth / 2.0, h5tb.param.h_mean + h5tb.param.h_wdth / 2.0);
    } else {
        throw except::runtime_error("wrong distribution [{}]: expected one of normal | lognormal | uniform", h5tb.param.distribution);
    }
    all_mpo_parameters_have_been_set = false;
    mpo_squared                      = std::nullopt;
    unique_id                        = std::nullopt;
    unique_id_sq                     = std::nullopt;
}

Eigen::Tensor<MpoSite::cplx, 4> IsingRandomField::MPO_nbody_view(std::optional<std::vector<size_t>>                  nbody,
                                                                 [[maybe_unused]] std::optional<std::vector<size_t>> skip) const {
    // This function returns a view of the MPO including only n-body terms.
    // For instance, if nbody_terms == {2,3}, this would exclude on-site terms.
    if(not nbody) return MPO();
    double J1 = 0, J2 = 0.0;
    for(const auto &n : nbody.value()) {
        if(n == 1) J1 = 1.0;
        if(n == 2) J2 = 1.0;
    }
    using namespace qm::spin::half;
    Eigen::Tensor<cplx, 4> MPO_nbody                                           = MPO();
    MPO_nbody.slice(std::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-J1 * get_field() * sx - e_shift * id);
    MPO_nbody.slice(std::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-J2 * get_coupling() * sz);
    return MPO_nbody;
}

Eigen::Tensor<MpoSite::cplx, 4> IsingRandomField::MPO_shifted_view() const { return MPO_shifted_view(e_shift); }

Eigen::Tensor<MpoSite::cplx, 4> IsingRandomField::MPO_shifted_view(double site_energy) const {
    using namespace qm::spin::half;
    if(site_energy == 0) { return MPO(); }
    Eigen::Tensor<cplx, 4> temp                                           = MPO();
    temp.slice(std::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(-get_field() * sx - site_energy * id);
    return temp;
}

std::unique_ptr<MpoSite> IsingRandomField::clone() const { return std::make_unique<IsingRandomField>(*this); }

long IsingRandomField::get_spin_dimension() const { return h5tb.param.spin_dim; }

void IsingRandomField::set_averages([[maybe_unused]] std::vector<TableMap> lattice_parameters, bool infinite) {
    if(not infinite) {
        lattice_parameters.back()["J1"]    = 0.0;
        lattice_parameters.back()["J2"]    = 0.0;
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

void IsingRandomField::save_hamiltonian(h5pp::File &file, std::string_view table_path) const {
    if(not file.linkExists(table_path)) file.createTable(h5tb_ising_tf_rf::h5_type, table_path, "Transverse-field Ising");
    file.appendTableRecords(h5tb, table_path);
    // Position 0 is also responsible for writing attributes
    if(position.value() != 0) return;

    file.writeAttribute(h5tb.param.J1, "J1", table_path);
    file.writeAttribute(h5tb.param.J2, "J2", table_path);
    file.writeAttribute(h5tb.param.h_mean, "h_mean", table_path);
    file.writeAttribute(h5tb.param.h_wdth, "h_wdth", table_path);
    file.writeAttribute(h5tb.param.h_tran, "h_tran", table_path);
    file.writeAttribute(h5tb.param.distribution, "distribution", table_path);
    file.writeAttribute(h5tb.param.spin_dim, "spin_dim", table_path);
}

void IsingRandomField::load_hamiltonian(const h5pp::File &file, std::string_view model_prefix) {
    auto ham_table = fmt::format("{}/hamiltonian", model_prefix);
    if(file.linkExists(ham_table)) {
        h5tb.param                       = file.readTableRecords<h5tb_ising_tf_rf::table>(ham_table, position);
        all_mpo_parameters_have_been_set = true;
    } else {
        throw except::runtime_error("Could not load MPO. Table [{}] does not exist", ham_table);
    }

    // Check that we are on the same point of the phase diagram
    using namespace settings::model::ising_tf_rf;
    if(std::abs(h5tb.param.J1 - J1) > 1e-6) throw except::runtime_error("J1_rand {:.16f} != {:.16f} ising_tf_rf::J1_rand", h5tb.param.J1, J1);
    if(std::abs(h5tb.param.J2 - J2) > 1e-6) throw except::runtime_error("J2_rand {:.16f} != {:.16f} ising_tf_rf::J2_rand", h5tb.param.J2, J2);
    if(std::abs(h5tb.param.h_tran - h_tran) > 1e-6) throw except::runtime_error("h_tran {:.16f} != {:.16f} ising_tf_rf::h_tran", h5tb.param.h_tran, h_tran);
    if(std::abs(h5tb.param.h_mean - h_mean) > 1e-6) throw except::runtime_error("h_mean {:.16f} != {:.16f} ising_tf_rf::h_mean", h5tb.param.h_mean, h_mean);
    if(std::abs(h5tb.param.h_wdth - h_wdth) > 1e-6) throw except::runtime_error("h_wdth {:.16f} != {:.16f} ising_tf_rf::h_wdth", h5tb.param.h_wdth, h_wdth);
}