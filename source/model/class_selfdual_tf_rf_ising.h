//
// Created by david on 2018-07-06.
//

#pragma once

#include "class_model_base.h"
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include <iostream>
#include <simulation/nmspc_settings.h>

class class_selfdual_tf_rf_ising : public class_model_base {
    using Scalar = std::complex<double>;

    private:
    size_t      spin_dim     = 0;           /*!< Spin dimension */
    std::string distribution = "lognormal"; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
    bool        parity_sep   = false;       /*!< Parity sector separation on/off */
    double      J_rnd        = 0;           /*!< Randomly distributed nearest neighbour coupling */
    double      h_rnd        = 0;           /*!< Randomly distributed on-site field */
    double      J_ptb        = 0;           /*!< Perturbation to the coupling, std::pow(J_rnd + J_ptb,1-alpha) */
    double      h_ptb        = 0;           /*!< Perturbation to the coupling, std::pow(J_rnd + J_ptb,1-alpha) */
    double      J_avg        = 0;           /*!< Average of J_rnd between all sites*/
    double      h_avg        = 0;           /*!< Average of h_rnd on all sites */
    double      J_mean       = 0;           /*!< Mean for the distrbution of J_rnd */
    double      h_mean       = 0;           /*!< Mean for the distrbution of h_rnd */
    double      J_sigma      = 0;           /*!< Standard deviation for the distribution of J_rnd */
    double      h_sigma      = 0;           /*!< Standard deviation for the distribution of h_rnd */
    double      lambda       = 0;           /*!< Factor involved in next-nearest neighbor interaction */
    double      delta        = 0;           /*!< Difference J_log_mean - h_log_mean  */
    double      alpha        = 0;           /*!< Damping factor [0,1] on random couplings, std::pow(J_rnd + J_ptb,1-alpha)  */
    double      beta         = 0;           /*!< Damping factor [0,1] on random fields, std::pow(h_rnd + h_ptb,1-alpha)  */
    double      psfactor     = 0;           /*!< Parity sector separation factor */

    double get_coupling() const;
    double get_field() const;
    double get_coupling(double J_rnd_, double J_ptb_, double alpha_) const;
    double get_field(double h_rnd_, double h_ptb_, double beta_) const;

    public:
    explicit class_selfdual_tf_rf_ising(size_t position_);

    // Functions that extend the base (no override)
    void set_realization_averages(double J_avg_, double h_avg_);
    // Functions that override the base
    std::unique_ptr<class_model_base> clone() const override;
    Eigen::Tensor<Scalar, 4>          MPO_reduced_view() const override;
    Eigen::Tensor<Scalar, 4>          MPO_reduced_view(double site_energy) const override;
    Eigen::Tensor<Scalar, 1>          get_MPO_edge_left() const override;
    Eigen::Tensor<Scalar, 1>          get_MPO_edge_right() const override;
    size_t                            get_spin_dimension() const override;
    Parameters                        get_parameters() const override;
    void                              set_parameters(const Parameters &parameters) override;
    void                              set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) override;
    void                              set_coupling_damping(double alpha) override;
    void                              set_field_damping(double beta) override;
    void                              build_mpo() override;
    void                              randomize_hamiltonian() override;
    bool                              is_perturbed() const override;
    bool                              is_damped() const override;
    void                              set_full_lattice_parameters(std::vector<Parameters> all_parameters, bool reverse = false) override;
    Eigen::MatrixXcd                  single_site_hamiltonian(int position, int sites, std::vector<Eigen::MatrixXcd> &SX, std::vector<Eigen::MatrixXcd> &SY,
                                                              std::vector<Eigen::MatrixXcd> &SZ) const override;
};
