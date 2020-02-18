//
// Created by david on 2018-07-04.
//

#include "class_tf_ising.h"
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <simulation/nmspc_settings.h>
#include <h5pp/h5pp.h>

using namespace qm::spinOneHalf;
using Scalar = std::complex<double>;

class_tf_ising::class_tf_ising(size_t position_) : class_model_base(position_) {
    log      = Logger::setLogger("tf-ising");
    pm.spin_dim = settings::model::tf_ising::d;
    pm.J_nn     = settings::model::tf_ising::J;
    pm.h_field  = settings::model::tf_ising::g;
    pm.h_sigma  = settings::model::tf_ising::w;
    std::strcpy(pm.distribution,
            settings::model::selfdual_tf_rf_ising::distribution.c_str());
    extent4 = {1, 1, pm.spin_dim, pm.spin_dim};
    extent2 = {pm.spin_dim, pm.spin_dim};

    qm::spinOneHalf::SX = qm::gen_manybody_spin(sx, 2);
    qm::spinOneHalf::SY = qm::gen_manybody_spin(sy, 2);
    qm::spinOneHalf::SZ = qm::gen_manybody_spin(sz, 2);
    qm::spinOneHalf::II = qm::gen_manybody_spin(Id, 2);

    h5table_define();
    all_mpo_parameters_have_been_set = true; // There are no full lattice parameters on this model so we set it true immediately!
    class_tf_ising::build_mpo();
}

double class_tf_ising::get_field() const { return pm.h_field + std::pow(pm.h_ptb + pm.h_rnd, 1 - beta); }
double class_tf_ising::get_coupling() const { return std::pow(pm.J_nn, 1 - alpha); }



void class_tf_ising::h5table_define() {
    // Create a type for the char array from the template H5T_C_S1
    // The template describes a string with a single char.
    // Set the size with H5Tset_size, or h5pp::hdf5::setStringSize(...)
    h5pp::hid::h5t h5t_custom_string = H5Tcopy(H5T_C_S1);
    H5Tset_size(h5t_custom_string, settings::model::model_type.size());

    // Optionally set the null terminator '\0'
    H5Tset_strpad(h5t_custom_string, H5T_STR_NULLTERM);


    // Define the compound type
    h5paramtype = H5Tcreate(H5T_COMPOUND, sizeof(p_tf_ising));
    H5Tinsert(h5paramtype, "J_nn", HOFFSET(p_tf_ising, J_nn), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "J_nnn", HOFFSET(p_tf_ising, J_nnn), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "h_field", HOFFSET(p_tf_ising, h_field), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "h_rnd", HOFFSET(p_tf_ising, h_rnd), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "h_ptb", HOFFSET(p_tf_ising, h_ptb), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "h_mean", HOFFSET(p_tf_ising, h_mean), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "h_sigma", HOFFSET(p_tf_ising, h_sigma), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "spin_dim", HOFFSET(p_tf_ising, spin_dim), H5T_NATIVE_ULONG);
    H5Tinsert(h5paramtype, "distribution", HOFFSET(p_tf_ising, distribution), h5t_custom_string);

}




void class_tf_ising::build_mpo()
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
    mpo_internal.resize(3, 3, pm.spin_dim, pm.spin_dim);
    mpo_internal.setZero();
    mpo_internal.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    mpo_internal.slice(Eigen::array<long, 4>{2, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 2, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
}

void class_tf_ising::randomize_hamiltonian() {
    pm.h_rnd = rn::uniform_double_box(-pm.h_sigma, pm.h_sigma);
    if(all_mpo_parameters_have_been_set or mpo_internal.size() > 3) {
        mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    }
}

void class_tf_ising::set_coupling_damping(double alpha_) {
    alpha                            = alpha_;
}
void class_tf_ising::set_field_damping(double beta_) {
    beta                             = beta_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    }
}

void class_tf_ising::set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) {
    switch(ptbMode) {
        case PerturbMode::ABSOLUTE: {
            pm.h_ptb = field_ptb;
            break;
        }
        case PerturbMode::PERCENTAGE: {
            pm.h_ptb = pm.h_rnd * field_ptb;
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_ABSOLUTE: {
            pm.h_ptb = rn::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_PERCENTAGE: {
            pm.h_ptb = pm.h_rnd * rn::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
    }
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
    }
    if(coupling_ptb == 0.0 and field_ptb == 0 and is_perturbed())
        throw std::runtime_error("MPO(" + std::to_string(get_position()) + ": Should have become unperturbed!");
}

bool class_tf_ising::is_perturbed() const { return pm.h_ptb != 0.0; }


Eigen::Tensor<Scalar, 4> class_tf_ising::MPO_reduced_view() const {
    if(e_reduced == 0) {
        return MPO();
    }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> class_tf_ising::MPO_reduced_view(double site_energy) const {
    if(site_energy == 0) {
        return MPO();
    }
    Eigen::Tensor<Scalar, 4> temp                                           = MPO();
    temp.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - site_energy * Id);
    return temp;
}

Eigen::MatrixXcd class_tf_ising::single_site_hamiltonian(int position, int sites, std::vector<Eigen::MatrixXcd> &SX,
                                                         std::vector<Eigen::MatrixXcd> &SY [[maybe_unused]], std::vector<Eigen::MatrixXcd> &SZ) const {
    int i = math::mod(position, sites);
    int j = math::mod(position + 1, sites);
    return -(pm.J_nn * SZ[i] * SZ[j] + pm.h_field * 0.5 * (SX[i] + SX[j]));
}

std::unique_ptr<class_model_base> class_tf_ising::clone() const { return std::make_unique<class_tf_ising>(*this); }


size_t class_tf_ising::get_spin_dimension() const { return pm.spin_dim; }

Eigen::Tensor<Scalar, 1> class_tf_ising::get_MPO_edge_left() const {
    Eigen::Tensor<Scalar, 1> ledge(3);
    ledge.setZero();
    ledge(2) = 1;
    return ledge;
}
Eigen::Tensor<Scalar, 1> class_tf_ising::get_MPO_edge_right() const {
    Eigen::Tensor<Scalar, 1> redge(3);
    redge.setZero();
    redge(0) = 1;
    return redge;
}

void class_tf_ising::set_parameters(const Parameters &parameters) {
    pm.J_nn         = get_val<double>(parameters, "J_nn");
    pm.J_nnn        = get_val<double>(parameters, "J_nnn");
    pm.h_field      = get_val<double>(parameters, "h_field");
    pm.h_rnd        = get_val<double>(parameters, "h_rnd");
    pm.h_mean       = get_val<double>(parameters, "h_mean");
    pm.h_sigma      = get_val<double>(parameters, "h_sigma");
    pm.spin_dim     = get_val<size_t>(parameters, "spin_dim");
    std::strcpy(pm.distribution,
                get_val<std::string>(parameters, "distribution").c_str());

    if(pm.J_nnn != 0.0) throw std::runtime_error("Use of [J_nnn] - Next-nearest neighbor coupling - is not implemented yet");
    all_mpo_parameters_have_been_set = true;
}

class_tf_ising::Parameters class_tf_ising::get_parameters() const {
    /* clang-format off */
    Parameters parameters;
    parameters.push_back({"J_nn", pm.J_nn});
    parameters.push_back({"J_nnn", pm.J_nnn});
    parameters.push_back({"h_field", pm.h_field});
    parameters.push_back({"h_rnd", pm.h_rnd});
    parameters.push_back({"h_ptb", pm.h_ptb});
    parameters.push_back({"h_mean", pm.h_mean});
    parameters.push_back({"h_sigma", pm.h_sigma});
    parameters.push_back({"spin_dim", pm.spin_dim});
    parameters.push_back({"distribution", std::string(pm.distribution)});

    return parameters;
    /* clang-format on */
}


void class_tf_ising::set_full_lattice_parameters([[maybe_unused]] std::vector<Parameters> lattice_parameters, bool reverse) {
    if(reverse) {
        std::reverse(lattice_parameters.begin(), lattice_parameters.end());
        for(size_t pos = 0; pos < lattice_parameters.size(); pos++)
        find_val(lattice_parameters[pos], "position") = pos;
    }
    find_val(lattice_parameters.back(), "J_nn")  = 0.0;
    find_val(lattice_parameters.back(), "J_nnn") = 0.0;
    // Recompute J_avg and h_avg from given J_rnd and h_rnd on all sites
    double J_sum = 0;
    double h_sum = 0;
    for(auto &site_params : lattice_parameters) {
        double J_nn_    = get_val<double>(site_params,"J_nn");
        double J_nnn_   = get_val<double>(site_params,"J_nnn");
        double h_field_ = get_val<double>(site_params,"h_field");
        double h_rnd_   = get_val<double>(site_params,"h_rnd");
        J_sum += J_nn_ + J_nnn_;
        h_sum += h_field_ + h_rnd_;
    }

    if(parity_sep) psfactor =  J_sum + h_sum;
    set_parameters(lattice_parameters[get_position()]);
}


void  class_tf_ising::write_parameters (h5pp::File & file, std::string_view table_name) const{
    if(not file.linkExists(table_name))
        file.createTable(h5paramtype,table_name, "Transverse-field Ising");

    file.appendTableEntries(pm,table_name);


}


void  class_tf_ising::read_parameters (h5pp::File & file, std::string_view table_name, size_t position ){
    file.readTableEntries(pm,table_name,position);
}
