//
// Created by david on 2018-07-06.
//

#pragma once

#include "class_model_base.h"
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include <iostream>
#include <simulation/nmspc_settings.h>

template<typename log_type>
class class_h5table_buffer;
class class_selfdual_tf_rf_ising_table;

class class_selfdual_tf_rf_ising_normal : public class_model_base {
    using Scalar = std::complex<double>;

    private:
    double J_rnd      = 0;      /*!< Random normal distributed nearest neighbour coupling */
    double h_rnd      = 0;      /*!< Random normal distributed on-site field */
    double J_ptb      = 1;      /*!< Perturbation to the coupling, std::pow(J_rnd + J_ptb,1-alpha) */
    double h_ptb      = 1;      /*!< Perturbation to the coupling, std::pow(J_rnd + J_ptb,1-alpha) */
    double J_avg      = 0;      /*!< Average of J_rnd between all sites*/
    double h_avg      = 0;      /*!< Average of h_rnd on all sites */
    double J_mean     = 0;      /*!< Mean for the normal distrbution of J_rnd */
    double h_mean     = 0;      /*!< Mean for the normal distrbution of h_rnd */
    double J_sigma    = 0;      /*!< Standard deviation for the normal distribution of J_rnd */
    double h_sigma    = 0;      /*!< Standard deviation for the normal distribution of h_rnd */
    double lambda     = 0;      /*!< Factor involved in next-nearest neighbor interaction */
    double delta      = 0;      /*!< Difference J_mean - h_mean  */
    double alpha      = 0;      /*!< Damping factor [0,1] on random couplings, std::pow(J_rnd + J_ptb,1-alpha)  */
    double beta       = 0;      /*!< Disorder damping [0,1] on random couplings, std::pow(h_rnd + h_ptb,1-alpha)  */
    int    spin_dim   = 0;      /*!< Spin dimension */
    int    num_params = 15 + 2; /*!< Number of parameters for this model excluding this one. Don't forget to add position and e_reduced  */

    double get_coupling() const;
    double get_field() const;
    double get_coupling(double J_rnd_, double J_ptb_, double alpha_) const;
    double get_field(double h_rnd_, double h_ptb_, double beta_) const;

    public:
    explicit class_selfdual_tf_rf_ising_normal(size_t position_, std::string logName = "SDUAL-ISING");

    // Functions that extend the base (no override)
    void set_realization_averages(double J_avg_, double h_avg_);
    auto get_J_rnd() { return J_rnd; };
    auto get_h_rnd() { return h_rnd; };

    // Functions that override the base
    std::unique_ptr<class_model_base> clone() const override;
    Eigen::Tensor<Scalar, 4>          MPO_reduced_view() const override;
    Eigen::Tensor<Scalar, 4>          MPO_reduced_view(double site_energy) const override;
    Eigen::Tensor<Scalar, 1>          get_MPO_edge_left() const override;
    Eigen::Tensor<Scalar, 1>          get_MPO_edge_right() const override;
    size_t                            get_spin_dimension() const override;
    std::map<std::string, double>     get_parameters() const override;
    void                              set_parameters(const std::map<std::string, double> & parameters) override;
    void                              set_hamiltonian(const Eigen::Tensor<Scalar, 4> &MPO_, std::vector<double> &parameters) override;
    void                              set_hamiltonian(const std::vector<double> &parameters) override;
    void                              set_hamiltonian(const Eigen::MatrixXd &all_parameters, int position) override;
    void                              set_hamiltonian(const Eigen::VectorXd &parameters) override;
    void                              set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) override;
    void                              set_coupling_damping(double alpha) override;
    void                              set_field_damping(double beta) override;
    void                              build_mpo() override;
    void                              randomize_hamiltonian() override;
    bool                              is_perturbed() const override;
    void                              set_full_lattice_parameters(std::vector<std::map<std::string,double>> chain_parameters, bool reverse = false) override;
    Eigen::MatrixXcd
        single_site_hamiltonian(int position, int sites, std::vector<Eigen::MatrixXcd> &SX, std::vector<Eigen::MatrixXcd> &SY, std::vector<Eigen::MatrixXcd> &SZ) const override;

    //    void   write_to_hdf5_table()                                                override;
};
