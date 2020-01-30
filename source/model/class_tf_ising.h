//
// Created by david on 2018-07-04.
//

#pragma once

#include "class_model_base.h"
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include <iostream>

class class_tf_ising : public class_model_base {
    using Scalar = std::complex<double>;

    private:
    size_t      spin_dim     = 0;           /*!< Spin dimension */
    std::string distribution = "lognormal"; /*!< The random distribution of r_rnd_field. Choose between lognormal, normal or uniform */
    bool        parity_sep   = false;       /*!< Parity sector separation on/off */
    double      J_nn         = 0;           /*!< Nearest neighbor coupling */
    double      J_nnn        = 0;           /*!< Next-nearest neighbor coupling */
    double      h_field      = 0;           /*!< On-site magnetic field */
    double      h_rnd        = 0;           /*!< Random field value */
    double      h_ptb        = 0;           /*!< Perturbation */
    double      h_mean       = 0;           /*!< Mean of the distribution for the random field */
    double      h_sigma      = 0;           /*!< Randomness strength. In distribution this is N(r_mean,r_sigma) or U(r_mean-r_sigma,r_mean+r_sigma) */
    double      alpha        = 0;           /*!< Damping factor [0,1] on couplings fields */
    double      beta         = 0;           /*!< Damping factor [0,1] on random fields */
    double      psfactor     = 0;           /*!< Parity sector separation factor */

    [[nodiscard]] double get_field() const;
    [[nodiscard]] double get_coupling() const;

    public:
    class_tf_ising(size_t position_);
    std::unique_ptr<class_model_base> clone() const override;
    Eigen::Tensor<Scalar, 4>          MPO_reduced_view() const override;
    Eigen::Tensor<Scalar, 4>          MPO_reduced_view(double site_energy) const override;
    Eigen::Tensor<Scalar, 1>          get_MPO_edge_left() const override;
    Eigen::Tensor<Scalar, 1>          get_MPO_edge_right() const override;
    size_t                            get_spin_dimension() const override;
    Parameters                        get_parameters() const override;
    void                              set_parameters(const Parameters &parameters) override;
    void                              register_h5_parameters() override;
    void                              set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) override;
    void                              set_coupling_damping(double alpha) override;
    void                              set_field_damping(double beta) override;
    void                              build_mpo() override;
    void                              randomize_hamiltonian() override;
    bool                              is_perturbed() const override;
    bool                              is_damped() const override;
    void                              set_full_lattice_parameters(std::vector<Parameters> lattice_parameters, bool reverse = false) override;
    Eigen::MatrixXcd                  single_site_hamiltonian(int position, int sites, std::vector<Eigen::MatrixXcd> &SX, std::vector<Eigen::MatrixXcd> &SY,
                                                              std::vector<Eigen::MatrixXcd> &SZ) const override;
};
