//
// Created by david on 2018-07-06.
//

#pragma once

#include "class_model_base.h"
#include "class_selfdual_tf_rf_ising.h"
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include <iostream>
#include <simulation/nmspc_settings.h>



class class_selfdual_tf_rf_ising_ps : public class_selfdual_tf_rf_ising {
    using Scalar = class_selfdual_tf_rf_ising::Scalar;
    double psfactor   = 0;      /*!< Parity separating factor (specific for this model) */
    public:
    explicit class_selfdual_tf_rf_ising_ps(size_t position_, std::string logName = "SDUAL-ISING");

    // Functions that extend the base (no override)
    void set_realization_averages(double J_avg_, double h_avg_);

    // Functions that override the base
    std::unique_ptr<class_model_base> clone() const override;
    Eigen::Tensor<Scalar, 4>          MPO_reduced_view() const override;
    Eigen::Tensor<Scalar, 4>          MPO_reduced_view(double site_energy) const override;
    Eigen::Tensor<Scalar, 1>          get_MPO_edge_left() const override;
    Eigen::Tensor<Scalar, 1>          get_MPO_edge_right() const override;
    std::map<std::string, double>     get_parameters() const override;
    void                              set_parameters(const std::map<std::string, double> &parameters) override;
    void                              set_hamiltonian(const Eigen::Tensor<Scalar, 4> &MPO_, std::vector<double> &parameters) override;
    void                              set_hamiltonian(const std::vector<double> &parameters) override;
    void                              set_hamiltonian(const Eigen::MatrixXd &all_parameters, int position) override;
    void                              set_hamiltonian(const Eigen::VectorXd &parameters) override;
    void                              build_mpo() override;
    void                              randomize_hamiltonian() override;
    void                              set_full_lattice_parameters(std::vector<std::map<std::string, double>> chain_parameters, bool reverse = false) override;
};
