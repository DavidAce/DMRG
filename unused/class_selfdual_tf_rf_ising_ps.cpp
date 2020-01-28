//
// Created by david on 2018-07-06.
//

#include "class_selfdual_tf_rf_ising_ps.h"
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <simulation/nmspc_settings.h>


using Scalar = std::complex<double>;

class_selfdual_tf_rf_ising_ps::class_selfdual_tf_rf_ising_ps(size_t position_, std::string logName) : class_model_base(position_, logName) {
    spin_dim   = settings::model::selfdual_tf_rf_ising_ps::d; /*!< Spin dimension */
    J_log_mean = settings::model::selfdual_tf_rf_ising_ps::J_log_mean;
    h_log_mean = settings::model::selfdual_tf_rf_ising_ps::h_log_mean;
    J_sigma    = settings::model::selfdual_tf_rf_ising_ps::J_sigma;
    h_sigma    = settings::model::selfdual_tf_rf_ising_ps::h_sigma;
    lambda     = settings::model::selfdual_tf_rf_ising_ps::lambda;
    e_reduced  = 0;

    extent4 = {1, 1, spin_dim, spin_dim};
    extent2 = {spin_dim, spin_dim};
    J_rnd   = rn::log_normal(J_log_mean, J_sigma);
    h_rnd   = rn::log_normal(h_log_mean, h_sigma);
    delta   = J_log_mean - h_log_mean;
}

double class_selfdual_tf_rf_ising_ps::get_coupling() const { return std::pow(J_rnd + J_ptb, 1 - alpha); }
double class_selfdual_tf_rf_ising_ps::get_field() const { return std::pow(h_rnd + h_ptb, 1 - beta); }
double class_selfdual_tf_rf_ising_ps::get_coupling(double J_rnd_, double J_ptb_, double alpha_) const { return std::pow(J_rnd_ + J_ptb_, 1 - alpha_); }
double class_selfdual_tf_rf_ising_ps::get_field(double h_rnd_, double h_ptb_, double beta_) const { return std::pow(h_rnd_ + h_ptb_, 1 - beta_); }


void class_selfdual_tf_rf_ising_ps::set_hamiltonian(const Eigen::Tensor<Scalar, 4> &MPO_, std::vector<double> &parameters) {
    mpo_internal = MPO_;
    set_hamiltonian(parameters);
    auto mpo1 = Eigen::Map<const Eigen::VectorXcd>(MPO_.data(), MPO_.size());
    auto mpo2 = Eigen::Map<const Eigen::VectorXcd>(MPO().data(), MPO().size());
    assert(mpo1 == mpo2 and "MPO mismatch!");
    if(mpo1 != mpo2) throw std::runtime_error("MPO mismatch");
}

void class_selfdual_tf_rf_ising_ps::set_hamiltonian(const std::vector<double> &parameters) {
    auto temp = Eigen::Map<const Eigen::VectorXd>(parameters.data(), parameters.size());
    set_hamiltonian(temp);
}

void class_selfdual_tf_rf_ising_ps::set_hamiltonian(const Eigen::MatrixXd &all_parameters, int position) { set_hamiltonian(all_parameters.row(position)); }

void class_selfdual_tf_rf_ising_ps::set_hamiltonian(const Eigen::VectorXd &parameters) {
    if((int) parameters.size() != num_params) throw std::runtime_error("Wrong number of parameters given to initialize this model");
    assert((int) parameters.size() == num_params and "ERROR: wrong number of parameters given to initialize this model");
    position                         = parameters(0);
    J_rnd                            = parameters(1);
    h_rnd                            = parameters(2);
    J_ptb                            = parameters(3);
    h_ptb                            = parameters(4);
    J_log_mean                       = parameters(5);
    h_log_mean                       = parameters(6);
    J_avg                            = parameters(7);
    h_avg                            = parameters(8);
    J_sigma                          = parameters(9);
    h_sigma                          = parameters(10);
    lambda                           = parameters(11);
    delta                            = parameters(12);
    psfactor                         = parameters(13);
    e_reduced                        = parameters(14);
    spin_dim                         = parameters(15);
    all_mpo_parameters_have_been_set = true;
    build_mpo();
}

void class_selfdual_tf_rf_ising_ps::set_realization_averages(double J_avg_, double h_avg_) {
    J_avg                            = J_avg_;
    h_avg                            = h_avg_;
    all_mpo_parameters_have_been_set = true;
    build_mpo();
}

void class_selfdual_tf_rf_ising_ps::build_mpo()
/*! Builds the MPO hamiltonian as a rank 4 tensor. Notation following Schollwöck (2010)
 *
 * H = - Σ J_{i} sz_{i} sz_{i+1} +  h_{i} sx_{i} + l*(h_avg sx_i sx_{i+1} + J_avg sz_{i} sz_{i+2}) + ps * Π sx_{i}
 *
 *  |     I                 0           0              0            0           0       |
 *  |     sz                0           0              0            0           0       |
 *  |     sx                0           0              0            0           0       |
 *  |     0                 I           0              0            0           0       |
 *  | -(h_rnd)*sx       -J_rnd*sz   -l*h_avg*sx     -l*J_avg*sz     I           0       |
 *  |     0                 0           0              0            0        ps * sx    |
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
    using namespace qm::spinOneHalf;
    mpo_internal.resize(6, 6, spin_dim, spin_dim);
    mpo_internal.setZero();
    mpo_internal.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sx);
    mpo_internal.slice(Eigen::array<long, 4>{3, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(h_rnd + h_ptb) * sx - e_reduced * Id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(J_rnd + J_ptb) * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(lambda * h_avg) * sx);
    mpo_internal.slice(Eigen::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(lambda * J_avg) * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 4, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(Id);
    mpo_internal.slice(Eigen::array<long, 4>{5, 5, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(sx); //Multiply the psfactor on the edge! Not on each MPO!
}

void class_selfdual_tf_rf_ising_ps::randomize_hamiltonian() {
    using namespace qm::spinOneHalf;
    J_rnd = rn::log_normal(J_log_mean, J_sigma);
    h_rnd = rn::log_normal(h_log_mean, h_sigma);
    if(all_mpo_parameters_have_been_set or mpo_internal.size() > 5) {
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(h_rnd + h_ptb) * sx - e_reduced * Id);
        mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(J_rnd + J_ptb) * sz);
    }
}

void class_selfdual_tf_rf_ising_ps::set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) {
    switch(ptbMode) {
        case PerturbMode::ABSOLUTE: {
            J_ptb = coupling_ptb;
            h_ptb = field_ptb;
            break;
        }
        case PerturbMode::PERCENTAGE: {
            J_ptb = J_rnd * coupling_ptb;
            h_ptb = h_rnd * field_ptb;
        }

        case PerturbMode::UNIFORM_RANDOM_ABSOLUTE: {
            J_ptb = rn::uniform_double_box(-coupling_ptb, coupling_ptb);
            h_ptb = rn::uniform_double_box(-field_ptb, field_ptb);
        }
        case PerturbMode::UNIFORM_RANDOM_PERCENTAGE: {
            J_ptb = J_rnd * rn::uniform_double_box(-coupling_ptb, coupling_ptb);
            h_ptb = h_rnd * rn::uniform_double_box(-field_ptb, field_ptb);
        }
    }
    if(all_mpo_parameters_have_been_set or mpo_internal.size() == 5 * 5 * spin_dim * spin_dim) {
        using namespace qm::spinOneHalf;
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_field() * sx - e_reduced * Id);
        mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-get_coupling() * sz);
    }
    if(coupling_ptb == 0.0 and field_ptb == 0 and is_perturbed()) throw std::runtime_error("MPO(" + std::to_string(get_position()) + ": Should have become unperturbed!");
}


bool class_selfdual_tf_rf_ising_ps::is_perturbed() const { return J_ptb != 0.0 or h_ptb != 0.0; }

Eigen::Tensor<Scalar, 4> class_selfdual_tf_rf_ising_ps::MPO_reduced_view() const {
    if(e_reduced == 0) { return MPO(); }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> class_selfdual_tf_rf_ising_ps::MPO_reduced_view(double site_energy) const {
    using namespace qm::spinOneHalf;
    if(site_energy == 0) { return MPO(); }
    Eigen::Tensor<Scalar, 4> temp                                           = MPO();
    temp.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::MatrixTensorMap(-(h_rnd + h_ptb) * sx - site_energy * Id);
    return temp;
}


Eigen::Tensor<Scalar, 1> class_selfdual_tf_rf_ising_ps::get_MPO_edge_left() const {
    Eigen::Tensor<Scalar, 1> ledge(6);
    ledge.setZero();
    ledge(4) = 1;
    ledge(5) = psfactor;
    return ledge;
}
Eigen::Tensor<Scalar, 1> class_selfdual_tf_rf_ising_ps::get_MPO_edge_right() const {
    Eigen::Tensor<Scalar, 1> redge(6);
    redge.setZero();
    redge(0) = 1;
    redge(5) = 1;
    return redge;
}

Eigen::MatrixXcd class_selfdual_tf_rf_ising_ps::single_site_hamiltonian(int                            position,
                                                                        int                            sites,
                                                                        std::vector<Eigen::MatrixXcd> &SX,
                                                                        std::vector<Eigen::MatrixXcd> &SY [[maybe_unused]],
                                                                        std::vector<Eigen::MatrixXcd> &SZ) const {
    int i = math::mod(position, sites);
    int j = math::mod(position + 1, sites);
    int k = math::mod(position + 2, sites);

    return -(J_rnd * SZ[i] * SZ[j] + h_rnd * 0.5 * (SX[i] + SX[j]) + lambda * (h_avg * SX[i] * SX[j] + J_avg * SZ[i] * SZ[k])); // TODO: How do we add prod(SX) here?
}

std::unique_ptr<class_model_base> class_selfdual_tf_rf_ising_ps::clone() const { return std::make_unique<class_selfdual_tf_rf_ising_ps>(*this); }

size_t class_selfdual_tf_rf_ising_ps::get_spin_dimension() const { return spin_dim; }

std::map<std::string,double> class_selfdual_tf_rf_ising_ps::get_parameters() const {
    /* clang-format off */
    std::map<std::string,double> parameters;
    parameters.insert({"position", get_position()});
    parameters.insert({"J_rnd", J_rnd});
    parameters.insert({"h_rnd", h_rnd});
    parameters.insert({"J_ptb", J_ptb});
    parameters.insert({"h_ptb", h_ptb});
    parameters.insert({"J_avg", J_avg});
    parameters.insert({"h_avg", h_avg});
    parameters.insert({"J_log_mean", J_log_mean});
    parameters.insert({"h_log_mean", h_log_mean});
    parameters.insert({"J_sigma", J_sigma});
    parameters.insert({"h_sigma", h_sigma});
    parameters.insert({"lambda", lambda});
    parameters.insert({"delta", delta});
    parameters.insert({"alpha", alpha});
    parameters.insert({"beta", beta});
    parameters.insert({"e_reduced", e_reduced});
    parameters.insert({"spin_dim", spin_dim});
    parameters.insert({"psfactor", psfactor});
    return parameters;
    /* clang-format on */
}

void class_selfdual_tf_rf_ising_ps::set_full_lattice_parameters(std::vector<std::map<std::string, double>> chain_parameters, bool reverse) {
    if(reverse) {
        std::reverse(chain_parameters.begin(), chain_parameters.end());
        for(size_t pos = 0; pos < chain_parameters.size(); pos++) {
            chain_parameters[pos].at("position") = pos;
            if(pos < chain_parameters.size() - 1) {
                chain_parameters[pos].at("J_rnd") = chain_parameters[pos + 1].at("J_rnd");
                chain_parameters[pos].at("J_ptb") = chain_parameters[pos + 1].at("J_ptb");
            } else {
                chain_parameters[pos].at("J_rnd") = 0;
                chain_parameters[pos].at("J_ptb") = 0;
            }
        }
    } else {
        chain_parameters.back().at("J_rnd") = 0;
        chain_parameters.back().at("J_ptb") = 0;
    }

    // Recompute J_avg and h_avg from given J_rnd and h_rnd on all sites
    double J_avg_ = 0;
    double h_avg_ = 0;
    for(auto &param : chain_parameters) {
        J_avg_ += get_coupling(param.at("J_rnd"), param.at("J_ptb"), param.at("alpha"));
        h_avg_ += get_coupling(param.at("h_rnd"), param.at("h_ptb"), param.at("beta"));
    }
    J_avg_ /= (chain_parameters.size() - 1);
    h_avg_ /= (chain_parameters.size());
    for(auto & param: chain_parameters){
        param.at("J_avg") = J_avg_;
        param.at("h_avg") = h_avg_;
    }

    set_parameters(chain_parameters[get_position()]);
}


