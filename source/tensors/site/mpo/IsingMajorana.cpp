#include "IsingMajorana.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "math/num.h"
#include "math/rnd.h"
#include "math/tenx.h"
#include "qm/spin.h"
#include "tools/common/log.h"
#include <h5pp/h5pp.h>

// double delta_to_J_mean(double delta) { return ;  }
// double delta_to_h_mean(double delta) { return ;  }
double delta_to_W_J(double delta) { return std::exp(delta / 2.0); }
double delta_to_W_h(double delta) { return std::exp(-delta / 2.0); }

IsingMajorana::IsingMajorana(ModelType model_type_, size_t position_) : MpoSite(model_type_, position_) {
    h5tb.param.g            = settings::model::ising_majorana::g;
    h5tb.param.delta        = settings::model::ising_majorana::delta;
    h5tb.param.spin_dim     = settings::model::ising_majorana::spin_dim;
    h5tb.param.distribution = settings::model::ising_majorana::distribution;
    extent4                 = {1, 1, h5tb.param.spin_dim, h5tb.param.spin_dim};
    extent2                 = {h5tb.param.spin_dim, h5tb.param.spin_dim};
}

double IsingMajorana::get_coupling() const { return h5tb.param.J_rand; }
double IsingMajorana::get_field() const { return h5tb.param.h_rand; }
void   IsingMajorana::print_parameter_names() const { h5tb.print_parameter_names(); }
void   IsingMajorana::print_parameter_values() const { h5tb.print_parameter_values(); }

void IsingMajorana::set_parameters(TableMap &parameters) {
    h5tb.param.g                     = std::any_cast<decltype(h5tb.param.g)>(parameters["g"]);
    h5tb.param.delta                 = std::any_cast<decltype(h5tb.param.delta)>(parameters["delta"]);
    h5tb.param.J_rand                = std::any_cast<decltype(h5tb.param.J_rand)>(parameters["J_rand"]);
    h5tb.param.h_rand                = std::any_cast<decltype(h5tb.param.h_rand)>(parameters["h_rand"]);
    h5tb.param.spin_dim              = std::any_cast<decltype(h5tb.param.spin_dim)>(parameters["spin_dim"]);
    h5tb.param.distribution          = std::any_cast<decltype(h5tb.param.distribution)>(parameters["distribution"]);
    all_mpo_parameters_have_been_set = true;
}

IsingMajorana::TableMap IsingMajorana::get_parameters() const {
    /* clang-format off */
    TableMap parameters;
    parameters["g"]             = h5tb.param.g;
    parameters["delta"]         = h5tb.param.delta;
    parameters["J_rand"]        = h5tb.param.J_rand;
    parameters["h_rand"]        = h5tb.param.h_rand;
    parameters["spin_dim"]      = h5tb.param.spin_dim;
    parameters["distribution"]  = h5tb.param.distribution;
    return parameters;
    /* clang-format on */
}

std::any IsingMajorana::get_parameter(const std::string_view name) const {
    /* clang-format off */
    if(name == "g")                 return  h5tb.param.g;
    else if(name == "delta")        return  h5tb.param.delta;
    else if(name == "J_rand")       return  h5tb.param.J_rand;
    else if(name == "h_rand")       return  h5tb.param.h_rand;
    else if(name == "spin_dim")     return  h5tb.param.spin_dim;
    else if(name == "distribution") return  h5tb.param.distribution;
    /* clang-format on */
    throw except::logic_error("Invalid parameter name for IsingMajorana model: {}", name);
}

void IsingMajorana::build_mpo()
/*! Builds the MPO hamiltonian as a rank 4 tensor.
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
 *  We use the definitions from https://link.aps.org/doi/10.1103/PhysRevResearch.4.L032016
 *  We take g > 0 and use box distributions for h_i and J_i in Box(0,W_h|J), where W_h|J are the
 *  disorder strengths of h_i and J_i.
 *  The disorder strengths are set to W_J = 1/W_h = W, giving us a single control parameter
 *  delta = <ln J_i> - <ln h_i> =  2ln(W).
 *
 */
{
    using namespace qm::spin::half;
    tools::log->debug("mpo({}): building ising-majorana mpo", get_position());
    if(not all_mpo_parameters_have_been_set)
        throw except::runtime_error("mpo({}): can't build mpo: full lattice parameters haven't been set yet.", get_position());
    mpo_internal.resize(5, 5, h5tb.param.spin_dim, h5tb.param.spin_dim);
    mpo_internal.setZero();
    mpo_internal.slice(std::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(id);
    mpo_internal.slice(std::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(sx);
    mpo_internal.slice(std::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(sz);
    mpo_internal.slice(std::array<long, 4>{3, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(id);
    mpo_internal.slice(std::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(get_field() * sz - e_shift * id);
    mpo_internal.slice(std::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(get_coupling() * sx);
    mpo_internal.slice(std::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(h5tb.param.g * sz);
    mpo_internal.slice(std::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(h5tb.param.g * sx);
    mpo_internal.slice(std::array<long, 4>{4, 4, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(id);

    if(parity_shift_sign_mpo != 0 and not parity_shift_axus_mpo.empty()) {
        // This redefines H --> H² + Q(σ), where
        //      * Q(σ) = 0.5 * ( I - prod(σ) ) = Proj(-σ), i.e. the "conjugate" projection operator (sign flipped).
        //      * σ is a pauli matrix (usually σ^z)
        // Observe that Q(σ)|ψ+-⟩ = (1 -+ 1) |ψ+-⟩
        // For ground state DMRG (fDMRG) we can add the projection on H directly, as
        //              H --> (H + Q(σ))
        //     such that
        //              (H + Q(σ)) |ψ+⟩ = (σ² + 0.5(1-1)) |ψ+⟩ = (E + 0) |ψ+⟩
        //              (H + Q(σ)) |ψ-⟩ = (σ² + 0.5(1+1)) |ψ-⟩ = (E + 1) |ψ-⟩
        auto d0 = mpo_internal.dimension(0);
        auto d1 = mpo_internal.dimension(1);
        auto d2 = mpo_internal.dimension(2);
        auto d3 = mpo_internal.dimension(3);
        auto pl = qm::spin::half::get_pauli(parity_shift_axus_mpo);

        Eigen::Tensor<cplx, 4> mpo_with_parity_shift_op(d0 + 2, d1 + 2, d2, d3);
        mpo_with_parity_shift_op.setZero();
        mpo_with_parity_shift_op.slice(tenx::array4{0, 0, 0, 0}, mpo_internal.dimensions())          = mpo_internal;
        mpo_with_parity_shift_op.slice(tenx::array4{d0, d1, 0, 0}, extent4).reshape(extent2)         = tenx::TensorMap(id);
        mpo_with_parity_shift_op.slice(tenx::array4{d0 + 1, d1 + 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(pl);
        mpo_internal                                                                                 = mpo_with_parity_shift_op;
    }

    if(tenx::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw except::runtime_error("mpo({}): found nan", get_position());
    }
    unique_id = std::nullopt;
}

void IsingMajorana::randomize_hamiltonian() {
    if(h5tb.param.distribution != "uniform") throw except::runtime_error("IsingMajorana expects a uniform distribution. Got: {}", h5tb.param.distribution);
    h5tb.param.J_rand = rnd::uniform_double_box(0, delta_to_W_J(h5tb.param.delta));
    h5tb.param.h_rand = rnd::uniform_double_box(0, delta_to_W_h(h5tb.param.delta));

    all_mpo_parameters_have_been_set = false;
    mpo_squared                      = std::nullopt;
    unique_id                        = std::nullopt;
    unique_id_sq                     = std::nullopt;
}

Eigen::Tensor<cplx, 4> IsingMajorana::MPO_nbody_view(std::optional<std::vector<size_t>> nbody, [[maybe_unused]] std::optional<std::vector<size_t>> skip) const {
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

Eigen::Tensor<cplx_t, 4> IsingMajorana::MPO_nbody_view_t([[maybe_unused]] std::optional<std::vector<size_t>> nbody,
                                                         [[maybe_unused]] std::optional<std::vector<size_t>> skip) const {
    throw except::runtime_error("IsingMajorana::MPO_nbody_view_t is not implemented");
}

Eigen::Tensor<cplx, 4> IsingMajorana::MPO_energy_shifted_view() const { return MPO_energy_shifted_view(e_shift); }

Eigen::Tensor<cplx, 4> IsingMajorana::MPO_energy_shifted_view(double energy_shift_per_site) const {
    using namespace qm::spin::half;
    Eigen::Tensor<cplx, 4> temp                                           = MPO();
    temp.slice(std::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = tenx::TensorCast(get_field() * sz - energy_shift_per_site * id);
    return temp;
}

std::unique_ptr<MpoSite> IsingMajorana::clone() const { return std::make_unique<IsingMajorana>(*this); }

long IsingMajorana::get_spin_dimension() const { return h5tb.param.spin_dim; }

void IsingMajorana::set_averages(std::vector<TableMap> all_parameters, bool infinite) {
    if(not infinite) { all_parameters.back()["J_rand"] = 0.0; }
    set_parameters(all_parameters[get_position()]);
}

void IsingMajorana::save_hamiltonian(h5pp::File &file, std::string_view hamiltonian_table_path) const {
    if(not file.linkExists(hamiltonian_table_path)) file.createTable(h5tb.get_h5_type(), hamiltonian_table_path, "Ising-Majorana");
    file.appendTableRecords(h5tb.param, hamiltonian_table_path);
}

void IsingMajorana::load_hamiltonian(const h5pp::File &file, std::string_view model_path) {
    auto ham_table = fmt::format("{}/hamiltonian", model_path);
    if(file.linkExists(ham_table)) {
        h5tb.param                       = file.readTableRecords<h5tb_ising_majorana::table>(ham_table, position);
        all_mpo_parameters_have_been_set = true;
    } else
        throw except::runtime_error("could not load mpo: table [{}] does not exist", ham_table);

    using namespace settings::model::ising_majorana;
    if(std::abs(h5tb.param.g - g) > 1e-6) throw except::runtime_error("g {:.16f} != {:.16f} ising_majorana::g", h5tb.param.g, g);
    if(std::abs(h5tb.param.delta - delta) > 1e-6) throw except::runtime_error("delta  {:.16f} != {:.16f} ising_majorana::delta", h5tb.param.delta, delta);

    build_mpo();
}
