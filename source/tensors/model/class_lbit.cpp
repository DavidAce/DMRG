//
// Created by david on 2018-07-04.
//

#include "class_lbit.h"
#include <config/nmspc_settings.h>
#include <general/nmspc_iter.h>
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
    h5tb.param.J2_base  = settings::model::lbit::J2_base;
    h5tb.param.J2_span  = settings::model::lbit::J2_span;
    h5tb.param.f_mixer  = settings::model::lbit::f_mixer;
    h5tb.param.u_layer  = settings::model::lbit::u_layer;
    h5tb.param.spin_dim = settings::model::lbit::spin_dim;
    if(h5tb.param.J2_span > 8) throw std::runtime_error(fmt::format("Maximum J2_span supported is 8 | Got: {}",h5tb.param.J2_span));
    copy_c_str(settings::model::lbit::distribution, h5tb.param.distribution);
    extent4 = {1, 1, h5tb.param.spin_dim, h5tb.param.spin_dim};
    extent2 = {h5tb.param.spin_dim, h5tb.param.spin_dim};
    h5tb_lbit::register_table_type();
    all_mpo_parameters_have_been_set =
        false; // There are no full lattice parameters but we set it to false so we remember to call randomize on all sites in every model type
}

// double class_lbit::get_field() const { return +std::pow(h5tb.param.J1_rand, 1 - beta); }
// double class_lbit::get_coupling() const { return std::pow(h5tb.param.J2_rand + h5tb.param.J3_rand, 1 - alpha); }
void class_lbit::print_parameter_names() const { h5tb_lbit::print_parameter_names(); }
void class_lbit::print_parameter_values() const { h5tb.print_parameter_values(); }

void class_lbit::set_parameters(TableMap &parameters) {
    h5tb.param.J1_rand  = std::any_cast<double>(parameters["J1_rand"]);
    h5tb.param.J2_rand  = std::any_cast<h5tb_lbit::J2Type>(parameters["J2_rand"]);
    h5tb.param.J3_rand  = std::any_cast<double>(parameters["J3_rand"]);
    h5tb.param.J1_wdth  = std::any_cast<double>(parameters["J1_wdth"]);
    h5tb.param.J2_wdth  = std::any_cast<double>(parameters["J2_wdth"]);
    h5tb.param.J3_wdth  = std::any_cast<double>(parameters["J3_wdth"]);
    h5tb.param.J2_base  = std::any_cast<double>(parameters["J2_base"]);
    h5tb.param.f_mixer  = std::any_cast<double>(parameters["f_mixer"]);
    h5tb.param.u_layer  = std::any_cast<size_t>(parameters["u_layer"]);
    h5tb.param.spin_dim = std::any_cast<long>(parameters["spin_dim"]);
    std::strcpy(h5tb.param.distribution, std::any_cast<std::string>(parameters["distribution"]).c_str());
    all_mpo_parameters_have_been_set = true;
    build_mpo();
}

class_lbit::TableMap class_lbit::get_parameters() const {
    /* clang-format off */
    TableMap parameters;
    parameters["J1_rand"]       = h5tb.param.J1_rand;
    parameters["J2_rand"]       = h5tb.param.J2_rand;
    parameters["J3_rand"]       = h5tb.param.J3_rand;
    parameters["J1_wdth"]       = h5tb.param.J1_wdth;
    parameters["J2_wdth"]       = h5tb.param.J2_wdth;
    parameters["J3_wdth"]       = h5tb.param.J3_wdth;
    parameters["J2_base"]       = h5tb.param.J2_base;
    parameters["f_mixer"]       = h5tb.param.f_mixer;
    parameters["u_layer"]       = h5tb.param.u_layer;
    parameters["spin_dim"]      = h5tb.param.spin_dim;
    parameters["distribution"]  = std::string(h5tb.param.distribution);
    return parameters;
    /* clang-format on */
}

void class_lbit::build_mpo()

/*
      Builds the MPO for the l-bit Hamiltonian as a rank-4 tensor.

      H =   Σ_i  J1(i)   * n_{i}
          + Σ_ij J2(i,j) * n_{i} * n_{i+j}
          + Σ_i  J3(i)   * n_{i} * n_{i+1} * n_{i+2}

      MPO:

            2            |    I     .     .     .     .     .     .     .     .    .    . |
            |            |    n     .     .     .     .     .     .     .     .    .    . |
        0---M---1    =   |    .     I     .     .     .     .     .     .     .    .    . |
            |            |    .     .     I     .     .     .     .     .     .    .    . |
            3            |    .     .     .     I     .     .     .     .     .    .    . |
                         |    .     .     .     .     I     .     .     .     .    .    . |
                         |    .     .     .     .     .     I     .     .     .    .    . |
                         |    .     .     .     .     .     .     I     .     .    .    . |
                         |    .     .     .     .     .     .     .     I     .    .    . |
                         |    .     n     .     .     .     .     .     .     .    .    . |
                         | J1*n J21*n J22*n J23*n J24*n J25*n J26*n J27*n J28*n J3*n    I |

       where
        *	i,j are site indices
        *   n = 0.5 * (1 + σ^z),
        *   σ^z is the diagonal 2x2 pauli matrix
        *   I is the 2x2 identity matrix
        *   J1,J2? and J3 are random 1,2 and 3-body couplings
        *   Jij = J21, J22... couples sites i,j at |i-j| <= 5
        *   Jij = J2_rand(i,j) = Uniform(-w+m,w+m) * b^-|i-j|, where
                * b is a suitable base, such as 5
                * w is the box width of the uniform distribution
                * m=mean is a constant offset of the distribution

       Built from the following finite-state-machine ([k] are matrix indices)

       I==[0]-------------J1*n(i)--------------[10]==I
           |                                    |
           |---n(i)---[1]----J21*n(i+1)---------|
                       | \                      |
                       |  I--[2]---J22*n(i+2)---|
                       |      I                 |
                       |     [3]---J23*n(i+3)---|
                       |      I                 |
                       |     [4]---J24*n(i+4)---|
                       |      I                 |
                       |     [5]---J25*n(i+5)---|
                       |      I                 |
                       |     [6]---J26*n(i+6)---|
                       |      I                 |
                       |     [7]---J27*n(i+7)---|
                       |      I                 |
                       |     [8]---J28*n(i+8)---|
                       |                        |
                       |-n(i+1)--[9]--J3*n(i+2)-|
  */

{
    if(not all_mpo_parameters_have_been_set) throw std::runtime_error("Improperly built MPO: Full lattice parameters haven't been set yet.");
    Eigen::Tensor<Scalar, 2> n = Textra::TensorCast(0.5 * (id + sz));
    Eigen::Tensor<Scalar, 2> i = Textra::TensorMap(id);
    mpo_internal.resize(11, 11, h5tb.param.spin_dim, h5tb.param.spin_dim);
    mpo_internal.setZero();
    mpo_internal.slice(Textra::array4{0, 0, 0, 0}, extent4).reshape(extent2)   = i;
    mpo_internal.slice(Textra::array4{1, 0, 0, 0}, extent4).reshape(extent2)   = n;
    mpo_internal.slice(Textra::array4{2, 1, 0, 0}, extent4).reshape(extent2)   = i;
    mpo_internal.slice(Textra::array4{3, 2, 0, 0}, extent4).reshape(extent2)   = i;
    mpo_internal.slice(Textra::array4{4, 3, 0, 0}, extent4).reshape(extent2)   = i;
    mpo_internal.slice(Textra::array4{5, 4, 0, 0}, extent4).reshape(extent2)   = i;
    mpo_internal.slice(Textra::array4{6, 5, 0, 0}, extent4).reshape(extent2)   = i;
    mpo_internal.slice(Textra::array4{7, 6, 0, 0}, extent4).reshape(extent2)   = i;
    mpo_internal.slice(Textra::array4{8, 7, 0, 0}, extent4).reshape(extent2)   = i;
    mpo_internal.slice(Textra::array4{9, 1, 0, 0}, extent4).reshape(extent2)   = n;
    mpo_internal.slice(Textra::array4{10, 0, 0, 0}, extent4).reshape(extent2)  = h5tb.param.J1_rand * n - e_reduced * i;
    for (auto && [r, J2r] : iter::enumerate(h5tb.param.J2_rand)){
        if(r == 0) continue;
        if(r > h5tb.param.J2_span) break;
        long rl = static_cast<long>(r);
        mpo_internal.slice(Textra::array4{10, rl, 0, 0}, extent4).reshape(extent2) = J2r * n;
    }



    mpo_internal.slice(Textra::array4{10, 1, 0, 0}, extent4).reshape(extent2)  = h5tb.param.J2_rand[1] * n;
    mpo_internal.slice(Textra::array4{10, 2, 0, 0}, extent4).reshape(extent2)  = h5tb.param.J2_rand[2] * n;
    mpo_internal.slice(Textra::array4{10, 3, 0, 0}, extent4).reshape(extent2)  = h5tb.param.J2_rand[3] * n;
    mpo_internal.slice(Textra::array4{10, 4, 0, 0}, extent4).reshape(extent2)  = h5tb.param.J2_rand[4] * n;
    mpo_internal.slice(Textra::array4{10, 5, 0, 0}, extent4).reshape(extent2)  = h5tb.param.J2_rand[5] * n;
    mpo_internal.slice(Textra::array4{10, 6, 0, 0}, extent4).reshape(extent2)  = h5tb.param.J2_rand[6] * n;
    mpo_internal.slice(Textra::array4{10, 7, 0, 0}, extent4).reshape(extent2)  = h5tb.param.J2_rand[7] * n;
    mpo_internal.slice(Textra::array4{10, 8, 0, 0}, extent4).reshape(extent2)  = h5tb.param.J2_rand[8] * n;
    mpo_internal.slice(Textra::array4{10, 9, 0, 0}, extent4).reshape(extent2)  = h5tb.param.J3_rand * n;
    mpo_internal.slice(Textra::array4{10, 10, 0, 0}, extent4).reshape(extent2) = i;
    if(Textra::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw std::runtime_error(fmt::format("MPO at position {} has NAN's", get_position()));
    }
    unique_id = std::nullopt;
    build_mpo_squared();
}

Eigen::Tensor<Scalar, 1> class_lbit::get_MPO_edge_left() const {
    auto                     ldim = mpo_internal.dimension(0);
    Eigen::Tensor<Scalar, 1> ledge(ldim);
    ledge.setZero();
    ledge(ldim - 1) = 1;
    return ledge;
}

Eigen::Tensor<Scalar, 1> class_lbit::get_MPO_edge_right() const {
    auto                     rdim = mpo_internal.dimension(1);
    Eigen::Tensor<Scalar, 1> redge(rdim);
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
    // J2(i,j) = J2_base^-|i-j| * Random_ij(J2_mean-J2_wdth_i/2, J2_mean+J2_wdth_i/2)
    using namespace settings::model::lbit;
    for(auto &&[j, J2] : iter::enumerate(h5tb.param.J2_rand)) J2 = std::pow(J2_base, -static_cast<double>(j)); // J2_base^-|i-j|

    if(std::string(h5tb.param.distribution) == "normal") {
        h5tb.param.J1_rand = rnd::normal(J1_mean, J1_wdth);
        h5tb.param.J3_rand = rnd::normal(J3_mean, J3_wdth);
        for(auto &J2 : h5tb.param.J2_rand) J2 *= rnd::normal(J2_mean, J2_wdth);
    } else if(std::string(h5tb.param.distribution) == "lognormal") {
        h5tb.param.J1_rand = rnd::log_normal(J1_mean, J1_wdth);
        h5tb.param.J3_rand = rnd::log_normal(J3_mean, J3_wdth);
        for(auto &J2 : h5tb.param.J2_rand) J2 *= rnd::log_normal(J2_mean, J2_wdth);
    } else if(std::string(h5tb.param.distribution) == "uniform") {
        h5tb.param.J1_rand = rnd::uniform_double_box(J1_mean - J1_wdth / 2.0, J1_mean + J1_wdth / 2.0);
        h5tb.param.J3_rand = rnd::uniform_double_box(J3_mean - J3_wdth / 2.0, J3_mean + J3_wdth / 2.0);
        for(auto &J2 : h5tb.param.J2_rand) J2 *= rnd::uniform_double_box(J2_mean - J2_wdth / 2.0, J2_mean + J2_wdth / 2.0);
    } else if(std::string(h5tb.param.distribution) == "constant") {
        h5tb.param.J1_rand = settings::model::lbit::J1_mean;
        h5tb.param.J3_rand = settings::model::lbit::J3_mean;
        for(auto &J2 : h5tb.param.J2_rand) J2 *= settings::model::lbit::J2_mean;
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
        Eigen::Tensor<Scalar, 2> n                                               = Textra::TensorCast(0.5 * (id + sz));
        Eigen::Tensor<Scalar, 2> i                                               = Textra::TensorMap(id);
        mpo_internal.slice(Textra::array4{3, 0, 0, 0}, extent4).reshape(extent2) = h5tb.param.J1_rand * n - e_reduced * i;
        mpo_squared                                                              = std::nullopt;
        unique_id                                                                = std::nullopt;
        unique_id_sq                                                             = std::nullopt;
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
        // double class_lbit::get_field() const { return +std::pow(h5tb.param.J1_rand, 1 - beta); }
        // double class_lbit::get_coupling() const { return std::pow(h5tb.param.J2_rand + h5tb.param.J3_rand, 1 - alpha); }

        Eigen::Tensor<Scalar, 2> n                                               = Textra::TensorCast(0.5 * (id + sz));
        Eigen::Tensor<Scalar, 2> i                                               = Textra::TensorMap(id);
        mpo_internal.slice(Textra::array4{3, 0, 0, 0}, extent4).reshape(extent2) = h5tb.param.J1_rand * n - e_reduced * i;
        mpo_squared                                                              = std::nullopt;
        unique_id                                                                = std::nullopt;
        unique_id_sq                                                             = std::nullopt;
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
    for(const auto &n : nbody_terms) {
        if(n == 1) J1 = 1;
        if(n == 2) J2 = 1;
        if(n == 3) J3 = 1;
    }
    Eigen::Tensor<Scalar, 4> MPO_nbody                                     = MPO();
    Eigen::Tensor<Scalar, 2> n                                             = Textra::TensorCast(0.5 * (id + sz));
    Eigen::Tensor<Scalar, 2> i                                             = Textra::TensorMap(id);
    for (auto && [r, J2r] : iter::enumerate(h5tb.param.J2_rand)){
        if(r == 0) continue;
        if(r > h5tb.param.J2_span) break;
        long rl = static_cast<long>(r);
        MPO_nbody.slice(Textra::array4{10, rl, 0, 0}, extent4).reshape(extent2) = J2 * J2r * n;
    }
    MPO_nbody.slice(Textra::array4{10, 9, 0, 0}, extent4).reshape(extent2) = J3 * h5tb.param.J3_rand * n;
    return MPO_nbody;
}

Eigen::Tensor<Scalar, 4> class_lbit::MPO_reduced_view() const {
    if(e_reduced == 0) { return MPO(); }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> class_lbit::MPO_reduced_view(double site_energy) const {
    if(site_energy == 0) { return MPO(); }
    Eigen::Tensor<Scalar, 4> temp                                        = MPO();
    long                     row                                         = temp.dimension(0) - 1;
    long                     col                                         = 0;
    Eigen::Tensor<Scalar, 2> n                                           = Textra::TensorCast(0.5 * (id + sz));
    Eigen::Tensor<Scalar, 2> i                                           = Textra::TensorMap(id);
    temp.slice(Textra::array4{row, col, 0, 0}, extent4).reshape(extent2) = h5tb.param.J1_rand * n - site_energy * i;
    return temp;
}

std::unique_ptr<class_mpo_site> class_lbit::clone() const { return std::make_unique<class_lbit>(*this); }

long class_lbit::get_spin_dimension() const { return h5tb.param.spin_dim; }

void class_lbit::set_averages([[maybe_unused]] std::vector<TableMap> lattice_parameters, bool infinite, bool reverse) {
    if(reverse) {
        std::reverse(lattice_parameters.begin(), lattice_parameters.end());
        for(size_t pos = 0; pos < lattice_parameters.size(); pos++) lattice_parameters[pos]["position"] = pos;
    }
    if(not infinite) {
        lattice_parameters.back()["J2_rand"]    = h5tb_lbit::J2Type{0};
        lattice_parameters.back()["J3_rand"]    = 0.0;
        lattice_parameters.end()[-2]["J3_rand"] = 0.0;
    }

    double J_sum = 0;
    for(auto &site_params : lattice_parameters) {
        auto J1_ = std::any_cast<double>(site_params["J1_rand"]);
        auto J3_ = std::any_cast<double>(site_params["J3_rand"]);
        auto J2_ = std::any_cast<h5tb_lbit::J2Type>(site_params["J2_rand"]);
        J_sum += J1_ + J3_;
        for(const auto &j2 : J2_) J_sum += j2;
    }
    if(parity_sep) psfactor = J_sum;
    set_parameters(lattice_parameters[get_position()]);
}

void class_lbit::save_hamiltonian(h5pp::File &file, const std::string &table_path) const {
    if(not file.linkExists(table_path)) file.createTable(h5tb_lbit::h5_type, table_path, "LBIT");
    file.appendTableRecords(h5tb, table_path);
    // Position 0 is also responsible for writing attributes
    if(position.value() != 0) return;

    file.writeAttribute(h5tb.param.J1_mean, "J1_mean", table_path);
    file.writeAttribute(h5tb.param.J2_mean, "J2_mean", table_path);
    file.writeAttribute(h5tb.param.J3_mean, "J3_mean", table_path);
    file.writeAttribute(h5tb.param.J1_wdth, "J1_wdth", table_path);
    file.writeAttribute(h5tb.param.J2_wdth, "J2_wdth", table_path);
    file.writeAttribute(h5tb.param.J3_wdth, "J3_wdth", table_path);
    file.writeAttribute(h5tb.param.J2_base, "J2_base", table_path);
    file.writeAttribute(h5tb.param.J2_span, "J2_span", table_path);
    file.writeAttribute(h5tb.param.f_mixer, "f_mixer", table_path);
    file.writeAttribute(h5tb.param.u_layer, "u_layer", table_path);
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
    if(std::abs(h5tb.param.J1_mean - J1_mean) > 1e-6)
        throw std::runtime_error(fmt::format("J1_mean {:.16f} != {:.16f} lbit::J1_mean", h5tb.param.J1_mean, J1_mean));
    if(std::abs(h5tb.param.J2_mean - J2_mean) > 1e-6)
        throw std::runtime_error(fmt::format("J2_mean {:.16f} != {:.16f} lbit::J2_mean", h5tb.param.J2_mean, J2_mean));
    if(std::abs(h5tb.param.J3_mean - J3_mean) > 1e-6)
        throw std::runtime_error(fmt::format("J3_mean {:.16f} != {:.16f} lbit::J3_mean", h5tb.param.J3_mean, J3_mean));
    if(std::abs(h5tb.param.J1_wdth - J1_wdth) > 1e-6)
        throw std::runtime_error(fmt::format("J1_wdth {:.16f} != {:.16f} lbit::J1_wdth", h5tb.param.J1_wdth, J1_wdth));
    if(std::abs(h5tb.param.J2_wdth - J2_wdth) > 1e-6)
        throw std::runtime_error(fmt::format("J2_wdth {:.16f} != {:.16f} lbit::J2_wdth", h5tb.param.J2_wdth, J2_wdth));
    if(std::abs(h5tb.param.J3_wdth - J3_wdth) > 1e-6)
        throw std::runtime_error(fmt::format("J3_wdth {:.16f} != {:.16f} lbit::J3_wdth", h5tb.param.J3_wdth, J3_wdth));
    if(std::abs(h5tb.param.J2_base - J2_base) > 1e-6)
        throw std::runtime_error(fmt::format("J2_base {:.16f} != {:.16f} lbit::J2_base", h5tb.param.J2_base, J2_base));
    if(h5tb.param.J2_span != J2_span )
        throw std::runtime_error(fmt::format("J2_span {} != {} lbit::J2_span", h5tb.param.J2_span, J2_span));
    if(std::abs(h5tb.param.f_mixer - f_mixer) > 1e-6)
        throw std::runtime_error(fmt::format("f_mixer {:.16f} != {:.16f} lbit::f_mixer", h5tb.param.f_mixer, f_mixer));

    if(h5tb.param.u_layer != u_layer) throw std::runtime_error(fmt::format("u_layer {:.16f} != {:.16f} lbit::u_layer", h5tb.param.u_layer, u_layer));

    // We can use the mpo's on file here to check everything is correct
    std::string mpo_path = fmt::format("{}/mpo/H_{}", model_prefix, get_position());
    if(file.linkExists(mpo_path)) {
        auto mpo_dset = file.readDataset<Eigen::Tensor<Scalar, 4>>(mpo_path);
        if(Textra::VectorCast(MPO()) != Textra::VectorCast(mpo_dset)) throw std::runtime_error("Built MPO does not match the MPO on file");
    }
}