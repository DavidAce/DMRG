//
// Created by david on 2018-07-06.
//

#include "class_selfdual_tf_rf_ising.h"
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_random_numbers.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <simulation/nmspc_settings.h>

using namespace qm::spinOneHalf;
using Scalar = std::complex<double>;


class_selfdual_tf_rf_ising::class_selfdual_tf_rf_ising(size_t position_, std::string logName): class_model_base(position_,logName){

    spin_dim            = settings::model::selfdual_tf_rf_ising::d;           /*!< Spin dimension */
    J_log_mean          = settings::model::selfdual_tf_rf_ising::J_log_mean;
    h_log_mean          = settings::model::selfdual_tf_rf_ising::h_log_mean;
    J_sigma             = settings::model::selfdual_tf_rf_ising::J_sigma;
    h_sigma             = settings::model::selfdual_tf_rf_ising::h_sigma;
    lambda              = settings::model::selfdual_tf_rf_ising::lambda;
    e_reduced           = 0;


    extent4     = {1, 1, spin_dim, spin_dim};
    extent2     = {spin_dim, spin_dim};
    J_rnd       = rn::log_normal(J_log_mean,J_sigma);
    h_rnd       = rn::log_normal(h_log_mean,h_sigma);
    delta       = J_log_mean - h_log_mean;

}


void   class_selfdual_tf_rf_ising::set_hamiltonian(const Eigen::Tensor<Scalar,4> & MPO_, std::vector<double> & parameters) {
    mpo_internal = MPO_;
    set_hamiltonian(parameters);
    auto mpo1 = Eigen::Map<const Eigen::VectorXcd>(MPO_ .data(),MPO_ .size());
    auto mpo2 = Eigen::Map<const Eigen::VectorXcd>(MPO().data(),MPO().size());
    assert(mpo1 == mpo2 and "MPO mismatch!");
    if(mpo1 != mpo2)throw std::runtime_error("MPO mismatch");
}

void   class_selfdual_tf_rf_ising::set_hamiltonian(const std::vector<double> & parameters) {
    auto temp = Eigen::Map<const Eigen::VectorXd>(parameters.data(),parameters.size());
    set_hamiltonian(temp);
}


void   class_selfdual_tf_rf_ising::set_hamiltonian(const Eigen::MatrixXd & all_parameters, int position) {
    set_hamiltonian (all_parameters.row(position));
}


void   class_selfdual_tf_rf_ising::set_hamiltonian(const Eigen::VectorXd & parameters) {
    if((int)parameters.size() != num_params ) throw std::runtime_error("Wrong number of parameters given to initialize this model");
    assert((int)parameters.size() == num_params and "ERROR: wrong number of parameters given to initialize this model");
    position       = parameters(0);
    J_rnd          = parameters(1);
    h_rnd          = parameters(2);
    J_log_mean     = parameters(3);
    h_log_mean     = parameters(4);
    J_avg          = parameters(5);
    h_avg          = parameters(6);
    J_sigma        = parameters(7);
    h_sigma        = parameters(8);
    lambda         = parameters(9);
    delta          = parameters(10);
    e_reduced      = parameters(11);
    spin_dim       = parameters(12);
    all_mpo_parameters_have_been_set = true;
    build_mpo();
}



void class_selfdual_tf_rf_ising::set_realization_averages(double J_avg_,double h_avg_){
    J_avg=J_avg_;
    h_avg=h_avg_;
    all_mpo_parameters_have_been_set = true;
    build_mpo();
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
    if (not all_mpo_parameters_have_been_set) throw std::runtime_error("Improperly built MPO: Full lattice parameters haven't been set yet.");
    mpo_internal.resize(5, 5, spin_dim, spin_dim);
    mpo_internal.setZero();
    mpo_internal.slice(Eigen::array<long, 4>{0, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(Id);
    mpo_internal.slice(Eigen::array<long, 4>{1, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(sz);
    mpo_internal.slice(Eigen::array<long, 4>{2, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(sx);
    mpo_internal.slice(Eigen::array<long, 4>{3, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(Id);
    mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-h_rnd * sx);
    mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-J_rnd * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 2, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-(lambda*h_avg) * sx);
    mpo_internal.slice(Eigen::array<long, 4>{4, 3, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-(lambda*J_avg) * sz);
    mpo_internal.slice(Eigen::array<long, 4>{4, 4, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(Id);
}

void class_selfdual_tf_rf_ising::randomize_hamiltonian(){
    J_rnd       = rn::log_normal(J_log_mean,J_sigma);
    h_rnd       = rn::log_normal(h_log_mean,h_sigma);
    if(all_mpo_parameters_have_been_set or mpo_internal.size()>5){
        mpo_internal.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-h_rnd * sx);
        mpo_internal.slice(Eigen::array<long, 4>{4, 1, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-J_rnd * sz);
    }
}


Eigen::Tensor<Scalar,4> class_selfdual_tf_rf_ising::MPO_reduced_view() const {
    if (e_reduced == 0){return MPO();}
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar,4> class_selfdual_tf_rf_ising::MPO_reduced_view(double site_energy) const {
    if (site_energy == 0){return MPO();}
    Eigen::Tensor<Scalar,4> temp  = MPO();
    temp.slice(Eigen::array<long, 4>{4, 0, 0, 0}, extent4).reshape(extent2) = Textra::Matrix_to_Tensor2(-h_rnd * sx - site_energy * Id);
    return temp;
}

Eigen::MatrixXcd class_selfdual_tf_rf_ising::single_site_hamiltonian(
        int position,
        int sites,
        std::vector<Eigen::MatrixXcd> &SX,
        std::vector<Eigen::MatrixXcd> &SY[[maybe_unused]],
        std::vector<Eigen::MatrixXcd> &SZ)
        const
{
    int i = math::mod(position,     sites);
    int j = math::mod(position + 1, sites);
    int k = math::mod(position + 2, sites);
    return - (J_rnd * SZ[i]*SZ[j] + h_rnd*0.5*(SX[i]+SX[j]) + lambda*(h_avg * SX[i]*SX[j] + J_avg*SZ[i]*SZ[k]));
}


std::unique_ptr<class_model_base> class_selfdual_tf_rf_ising::clone() const {return std::make_unique<class_selfdual_tf_rf_ising>(*this);}


size_t class_selfdual_tf_rf_ising::get_spin_dimension()                     const {return spin_dim;}
//double class_selfdual_tf_rf_ising::get_energy_reduced()                        const {return e_reduced;}
//double class_selfdual_tf_rf_ising::get_random_field()                          const {return h_rnd;}
//double class_selfdual_tf_rf_ising::get_randomness_strength()                   const {return w_rnd_strength;}


void class_selfdual_tf_rf_ising::print_parameter_names() const {
    std::cout
            << std::setprecision(10)
            << std::setw(16) << std::left << "MPO #"
            << std::setw(16) << std::left << "J_rnd"
            << std::setw(16) << std::left << "h_rnd"
            << std::setw(16) << std::left << "J_log_mean"
            << std::setw(16) << std::left << "h_log_mean"
            << std::setw(16) << std::left << "J_avg"
            << std::setw(16) << std::left << "h_avg"
            << std::setw(16) << std::left << "J_sigma"
            << std::setw(16) << std::left << "h_sigma"
            << std::setw(16) << std::left << "lambda"
            << std::setw(16) << std::left << "delta"
            << std::setw(16) << std::left << "e_reduced"
            << std::setw(16) << std::left << "spin_dim"
            << std::endl;
}

void class_selfdual_tf_rf_ising::print_parameter_values() const {
    std::cout
            << std::setprecision(10)
            << std::setw(16) << std::left << get_position()
            << std::setw(16) << std::left << J_rnd
            << std::setw(16) << std::left << h_rnd
            << std::setw(16) << std::left << J_log_mean
            << std::setw(16) << std::left << h_log_mean
            << std::setw(16) << std::left << J_avg
            << std::setw(16) << std::left << h_avg
            << std::setw(16) << std::left << J_sigma
            << std::setw(16) << std::left << h_sigma
            << std::setw(16) << std::left << lambda
            << std::setw(16) << std::left << delta
            << std::setw(16) << std::left << e_reduced
            << std::setw(16) << std::left << spin_dim
            << std::endl;
}

std::vector<std::string> class_selfdual_tf_rf_ising::get_parameter_names() const {
    return {"position",
            "J_rnd",
            "h_rnd",
            "J_log_mean",
            "h_log_mean",
            "J_avg",
            "h_avg",
            "J_sigma",
            "h_sigma",
            "lambda",
            "delta",
            "e_reduced",
            "spin_dim"
    };
}

//std::vector<double> class_selfdual_tf_rf_ising::get_random_parameter_values() const {
//
//    return {(double)get_position(),
//            J_rnd,
//            h_rnd,
//            J_log_mean,
//            h_log_mean,
//            J_avg,
//            h_avg,
//            J_sigma,
//            h_sigma,
//            lambda,
//            delta,
//            e_reduced,
//            (double)spin_dim
//    };
//}

std::vector<double> class_selfdual_tf_rf_ising::get_parameter_values() const {
    return {(double)get_position(),
            J_rnd,
            h_rnd,
            J_log_mean,
            h_log_mean,
            J_avg,
            h_avg,
            J_sigma,
            h_sigma,
            lambda,
            delta,
            e_reduced,
            (double)spin_dim
            };
}



void class_selfdual_tf_rf_ising::set_full_lattice_parameters(const std::vector<std::vector<double>> chain_parameters, bool reverse){
    // Calculate average J_rnd on the whole state
    all_mpo_parameters_have_been_set = true;
    std::list<double> J_rnd_list;
    std::list<double> h_rnd_list;
    if(reverse){
        J_rnd_list.emplace_front(0.0);
        for (auto &params : chain_parameters){
            J_rnd_list.emplace_front(params[1]);
            h_rnd_list.emplace_front(params[2]);
        }
        J_rnd_list.pop_front();
    }else{
        for (auto &params : chain_parameters){
            J_rnd_list.push_back(params[1]);
            h_rnd_list.push_back(params[2]);
        }
        J_rnd_list.back() = 0.0;
    }

    J_rnd = *std::next(J_rnd_list.begin(), get_position());
    h_rnd = *std::next(h_rnd_list.begin(), get_position());
    J_rnd_list.pop_back();   // Take average of all minus last
//    h_rnd_list.pop_back(); // Take average of all

    double J_rnd_avg = std::accumulate(J_rnd_list.begin(),J_rnd_list.end(),0.0)/J_rnd_list.size();
    double h_rnd_avg = std::accumulate(h_rnd_list.begin(),h_rnd_list.end(),0.0)/h_rnd_list.size();
    set_realization_averages(J_rnd_avg,h_rnd_avg);

}



//void   class_selfdual_tf_rf_ising::write_to_hdf5_table(){
//
//}
