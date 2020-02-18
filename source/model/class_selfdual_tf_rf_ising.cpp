//
// Created by david on 2018-07-06.
//

#include "class_selfdual_tf_rf_ising.h"
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <simulation/nmspc_settings.h>
#include <h5pp/h5pp.h>

using namespace qm::spinOneHalf;
using Scalar = std::complex<double>;

class_selfdual_tf_rf_ising::class_selfdual_tf_rf_ising(size_t position_) : class_model_base(position_) {
    log             = Logger::setLogger("selfdual-tf-rf-ising");
    pm.J_mean       = settings::model::selfdual_tf_rf_ising::J_mean;
    pm.h_mean       = settings::model::selfdual_tf_rf_ising::h_mean;
    pm.J_sigma      = settings::model::selfdual_tf_rf_ising::J_sigma;
    pm.h_sigma      = settings::model::selfdual_tf_rf_ising::h_sigma;
    pm.lambda       = settings::model::selfdual_tf_rf_ising::lambda;
    pm.delta        = pm.J_mean - pm.h_mean;
    pm.spin_dim     = settings::model::selfdual_tf_rf_ising::d;
    std::strcpy(pm.distribution,
            settings::model::selfdual_tf_rf_ising::distribution.c_str());
    parity_sep   = settings::model::selfdual_tf_rf_ising::parity_sep;

    extent4      = {1, 1, pm.spin_dim, pm.spin_dim};
    extent2      = {pm.spin_dim, pm.spin_dim};
    h5table_define();
    class_selfdual_tf_rf_ising::randomize_hamiltonian();
}



double class_selfdual_tf_rf_ising::get_coupling() const { return std::pow(pm.J_rnd + pm.J_ptb, 1 - alpha); }
double class_selfdual_tf_rf_ising::get_field() const { return std::pow(pm.h_rnd + pm.h_ptb, 1 - beta); }
double class_selfdual_tf_rf_ising::get_coupling(double J_rnd_, double J_ptb_, double alpha_) const { return std::pow(J_rnd_ + J_ptb_, 1 - alpha_); }
double class_selfdual_tf_rf_ising::get_field(double h_rnd_, double h_ptb_, double beta_) const { return std::pow(h_rnd_ + h_ptb_, 1 - beta_); }

void class_selfdual_tf_rf_ising::h5table_define() {
    // Create a type for the char array from the template H5T_C_S1
    // The template describes a string with a single char.
    // Set the size with H5Tset_size, or h5pp::hdf5::setStringSize(...)
    h5pp::hid::h5t h5t_custom_string = H5Tcopy(H5T_C_S1);
    H5Tset_size(h5t_custom_string, 10);

    // Optionally set the null terminator '\0'
    H5Tset_strpad(h5t_custom_string, H5T_STR_NULLTERM);


    // Define the compound type
    h5paramtype = H5Tcreate(H5T_COMPOUND, sizeof(p_selfdual_tf_rf_ising));
    H5Tinsert(h5paramtype, "J_rnd", HOFFSET(p_selfdual_tf_rf_ising, J_rnd), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "h_rnd", HOFFSET(p_selfdual_tf_rf_ising, h_rnd), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "J_ptb", HOFFSET(p_selfdual_tf_rf_ising, J_ptb), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "h_ptb", HOFFSET(p_selfdual_tf_rf_ising, h_ptb), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "J_avg", HOFFSET(p_selfdual_tf_rf_ising, J_avg), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "h_avg", HOFFSET(p_selfdual_tf_rf_ising, h_avg), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "J_mean", HOFFSET(p_selfdual_tf_rf_ising, J_mean), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "h_mean", HOFFSET(p_selfdual_tf_rf_ising, h_mean), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "J_sigma", HOFFSET(p_selfdual_tf_rf_ising, J_sigma), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "h_sigma", HOFFSET(p_selfdual_tf_rf_ising, h_sigma), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "lambda", HOFFSET(p_selfdual_tf_rf_ising, lambda), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "delta", HOFFSET(p_selfdual_tf_rf_ising, delta), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5paramtype, "spin_dim", HOFFSET(p_selfdual_tf_rf_ising, spin_dim), H5T_NATIVE_ULONG);
    H5Tinsert(h5paramtype, "distribution", HOFFSET(p_selfdual_tf_rf_ising, distribution), h5t_custom_string);

}


void class_selfdual_tf_rf_ising::set_parameters(const Parameters &parameters) {
    pm.J_rnd        = get_val<double>(parameters, "J_rnd");
    pm.J_ptb        = get_val<double>(parameters, "J_ptb");
    pm.h_rnd        = get_val<double>(parameters, "h_rnd");
    pm.h_ptb        = get_val<double>(parameters, "h_ptb");
    pm.J_avg        = get_val<double>(parameters, "J_avg");
    pm.h_avg        = get_val<double>(parameters, "h_avg");
    pm.J_mean       = get_val<double>(parameters, "J_mean");
    pm.h_mean       = get_val<double>(parameters, "h_mean");
    pm.J_sigma      = get_val<double>(parameters, "J_sigma");
    pm.h_sigma      = get_val<double>(parameters, "h_sigma");
    pm.lambda       = get_val<double>(parameters, "lambda");
    pm.delta        = get_val<double>(parameters, "delta");
    pm.spin_dim     = get_val<size_t>(parameters, "spin_dim");
    std::strcpy(pm.distribution,
                get_val<std::string>(parameters, "distribution").c_str());
    all_mpo_parameters_have_been_set = true;
    build_mpo();
}

class_selfdual_tf_rf_ising::Parameters class_selfdual_tf_rf_ising::get_parameters() const {
    /* clang-format off */
    Parameters parameters;
    parameters.push_back({"J_rnd", pm.J_rnd});
    parameters.push_back({"h_rnd", pm.h_rnd});
    parameters.push_back({"J_ptb", pm.J_ptb});
    parameters.push_back({"h_ptb", pm.h_ptb});
    parameters.push_back({"J_avg", pm.J_avg});
    parameters.push_back({"h_avg", pm.h_avg});
    parameters.push_back({"J_mean", pm.J_mean});
    parameters.push_back({"h_mean", pm.h_mean});
    parameters.push_back({"J_sigma", pm.J_sigma});
    parameters.push_back({"h_sigma", pm.h_sigma});
    parameters.push_back({"lambda", pm.lambda});
    parameters.push_back({"delta",  pm.delta});
    parameters.push_back({"spin_dim", pm.spin_dim});
    parameters.push_back({"distribution", std::string(pm.distribution)});
    return parameters;
    /* clang-format on */
}


void class_selfdual_tf_rf_ising::build_mpo()
/*! Builds the MPO hamiltonian as a rank 4 tensor. Notation following Schollwöck (2010)
 *
 * H = - Σ J_{i} sz_{i} sz_{i+1} +  h_{i} sx_{i} + l*(h_avg sx_i sx_{i+1} + J_avg sz_{i} sz_{i+2})
 *
 *  |     I                 0           0              0            0   |
 *  |     sz                0           0              0            0   |
 *  |     sx                0           0              0            0   |
 *  |     0                 I           0              0            0   |
 *  | -(h_rnd)*sx       -J_rnd*sz   -l*h_avg*sx     -l*J_avg*sz     I   |
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
        mpo_internal.resize(6, 6, pm.spin_dim, pm.spin_dim);
        mpo_internal.setZero();
        mpo_internal.slice(Eigen::array<long, 4>{5, 5, 0, 0}, extent4).reshape(extent2) =
            Textra::MatrixTensorMap(sx); // Multiply the psfactor on the edge! Not on each MPO!
    } else {
        mpo_internal.resize(5, 5, pm.spin_dim, pm.spin_dim);
        mpo_internal.setZero();
    }

    mpo_internal.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sx);
    mpo_internal.slice(Eigen::array<long, 4>{3, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx -  e_reduced * Id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(pm.lambda * pm.h_avg) * sx);
    mpo_internal.slice(Eigen::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(pm.lambda * pm.J_avg) * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 4, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);

    if(Textra::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw std::runtime_error("MPO at position " + std::to_string(get_position()) + " has been constructed with NAN's");
    }
}

void class_selfdual_tf_rf_ising::randomize_hamiltonian() {
    if(std::string(pm.distribution) == "normal") {
        pm.J_rnd = rn::normal(pm.J_mean, pm.J_sigma);
        pm.h_rnd = rn::normal(pm.h_mean, pm.h_sigma);
    } else if(std::string(pm.distribution) == "lognormal") {
        pm.J_rnd = rn::log_normal(pm.J_mean, pm.J_sigma);
        pm.h_rnd = rn::log_normal(pm.h_mean, pm.h_sigma);
    } else if(std::string(pm.distribution) == "uniform") {
        pm.J_rnd = rn::uniform_double_box(pm.J_mean - pm.J_sigma / 2.0, pm.J_mean + pm.J_sigma / 2.0);
        pm.h_rnd = rn::uniform_double_box(pm.h_mean - pm.h_sigma / 2.0, pm.h_mean + pm.J_sigma / 2.0);
    } else {
        throw std::runtime_error("Wrong distribution given. Expected one of <normal>, <lognormal>, <uniform>");
    }

    all_mpo_parameters_have_been_set = false;
}

void class_selfdual_tf_rf_ising::set_coupling_damping(double alpha_) {
    alpha = alpha_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    }
}
void class_selfdual_tf_rf_ising::set_field_damping(double beta_) {
    beta = beta_;
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx -  e_reduced * Id);
    }
}

void class_selfdual_tf_rf_ising::set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) {
    switch(ptbMode) {
        case PerturbMode::ABSOLUTE: {
            pm.J_ptb = coupling_ptb;
            pm.h_ptb = field_ptb;
            break;
        }
        case PerturbMode::PERCENTAGE: {
            pm.J_ptb = pm.J_rnd * coupling_ptb;
            pm.h_ptb = pm.h_rnd * field_ptb;
            break;
        }

        case PerturbMode::UNIFORM_RANDOM_ABSOLUTE: {
            pm.J_ptb = rn::uniform_double_box(-coupling_ptb, coupling_ptb);
            pm.h_ptb = rn::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_PERCENTAGE: {
            pm.J_ptb = pm.J_rnd * rn::uniform_double_box(-coupling_ptb, coupling_ptb);
            pm.h_ptb = pm.h_rnd * rn::uniform_double_box(-field_ptb, field_ptb);
            break;
        }
    }
    if(all_mpo_parameters_have_been_set) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx -  e_reduced * Id);
        mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    }
    if(coupling_ptb == 0.0 and field_ptb == 0 and is_perturbed())
        throw std::runtime_error("MPO(" + std::to_string(get_position()) + ": Should have become unperturbed!");
}

bool class_selfdual_tf_rf_ising::is_perturbed() const { return pm.J_ptb != 0.0 or pm.h_ptb != 0.0; }



Eigen::Tensor<Scalar, 4> class_selfdual_tf_rf_ising::MPO_reduced_view() const {
    if(e_reduced == 0) {
        return MPO();
    }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> class_selfdual_tf_rf_ising::MPO_reduced_view(double site_energy) const {
    if(site_energy == 0) {
        return MPO();
    }
    Eigen::Tensor<Scalar, 4> temp                                           = MPO();
    temp.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - site_energy * Id);
    temp.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    return temp;
}

Eigen::Tensor<Scalar, 1> class_selfdual_tf_rf_ising::get_MPO_edge_left() const {
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

Eigen::Tensor<Scalar, 1> class_selfdual_tf_rf_ising::get_MPO_edge_right() const {
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

Eigen::MatrixXcd class_selfdual_tf_rf_ising::single_site_hamiltonian(int position, int sites, std::vector<Eigen::MatrixXcd> &SX,
                                                                     std::vector<Eigen::MatrixXcd> &SY [[maybe_unused]],
                                                                     std::vector<Eigen::MatrixXcd> &SZ) const {
    int i = math::mod(position, sites);
    int j = math::mod(position + 1, sites);
    int k = math::mod(position + 2, sites);
    return -(pm.J_rnd * SZ[i] * SZ[j] + pm.h_rnd * 0.5 * (SX[i] + SX[j]) + pm.lambda * (pm.h_avg * SX[i] * SX[j] + pm.J_avg * SZ[i] * SZ[k]));
}

std::unique_ptr<class_model_base> class_selfdual_tf_rf_ising::clone() const { return std::make_unique<class_selfdual_tf_rf_ising>(*this); }

size_t class_selfdual_tf_rf_ising::get_spin_dimension() const { return pm.spin_dim; }

void class_selfdual_tf_rf_ising::set_full_lattice_parameters(std::vector<Parameters> all_parameters, bool reverse) {
    if(reverse) {
        std::reverse(all_parameters.begin(), all_parameters.end());
        for(size_t pos = 0; pos < all_parameters.size(); pos++) {
            find_val(all_parameters[pos], "position") = pos;
            if(pos < all_parameters.size() - 1) {
                find_val(all_parameters[pos], "J_rnd") = find_val(all_parameters[pos + 1], "J_rnd");
                find_val(all_parameters[pos], "J_ptb") = find_val(all_parameters[pos + 1], "J_ptb");
            } else {
                find_val(all_parameters[pos], "J_rnd") = 0.0;
                find_val(all_parameters[pos], "J_ptb") = 0.0;
            }
        }
    } else {
        find_val(all_parameters.back(), "J_rnd") = 0.0;
        find_val(all_parameters.back(), "J_ptb") = 0.0;
    }

    // Recompute J_avg and pm.h_avg from given pm.J_rnd and pm.h_rnd on all sites
    double J_avg_ = 0;
    double h_avg_ = 0;
    for(auto &site_params : all_parameters) {
        J_avg_ += get_val<double>(site_params, "J_rnd");
        h_avg_ += get_val<double>(site_params, "h_rnd");
    }
    J_avg_ /= (double)(all_parameters.size() - 1);
    h_avg_ /= (double)(all_parameters.size());
    for(auto &site_params : all_parameters) {
        find_val(site_params, "J_avg") = J_avg_;
        find_val(site_params, "h_avg") = h_avg_;
    }
    if(parity_sep) psfactor = all_parameters.size() * (J_avg_ + h_avg_) * (1.0 + pm.lambda);
    set_parameters(all_parameters[get_position()]);
}



void  class_selfdual_tf_rf_ising::write_parameters (h5pp::File & file, std::string_view table_name) const{
    if(not file.linkExists(table_name))
        file.createTable(h5paramtype,table_name, "Selfdual Ising");

    file.appendTableEntries(pm,table_name);


}
void  class_selfdual_tf_rf_ising::read_parameters (h5pp::File & file, std::string_view table_name, size_t position ){
    file.readTableEntries(pm,table_name,position);
}
