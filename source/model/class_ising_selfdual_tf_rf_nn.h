//
// Created by david on 2018-07-06.
//

#pragma once

#include "class_model_base.h"
#include <general/nmspc_tensor_extra.h>
#include <h5pp/details/h5ppHid.h>

namespace h5pp::hid {
    class h5t;
}

class class_ising_selfdual_tf_rf_nn : public class_model_base {
    using Scalar = std::complex<double>;

    private:
    h5tb_ising_selfdual_tf_rf_nn::table pm;
    [[nodiscard]] double        get_coupling() const;
    [[nodiscard]] double        get_field() const;
    [[nodiscard]] double        get_coupling(double J_rnd_, double J_ptb_, double alpha_) const;
    [[nodiscard]] double        get_field(double h_rnd_, double h_ptb_, double beta_) const;

    public:
    explicit class_ising_selfdual_tf_rf_nn(std::string_view model_type,size_t position_);

    // Functions that extend the base (no override)
    void set_realization_averages(double J_avrg_, double h_avrg_);
    // Functions that override the base
    [[nodiscard]] std::unique_ptr<class_model_base> clone() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 4>          MPO_reduced_view() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 4>          MPO_reduced_view(double site_energy) const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 1>          get_MPO_edge_left() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 1>          get_MPO_edge_right() const override;
    [[nodiscard]] size_t                            get_spin_dimension() const override;
    [[nodiscard]] TableMap                          get_parameters() const override;
    [[nodiscard]] bool                              is_perturbed() const override;
    [[nodiscard]] Eigen::MatrixXcd single_site_hamiltonian(size_t position, size_t sites, std::vector<Eigen::MatrixXcd> &SX, std::vector<Eigen::MatrixXcd> &SY,
                                                           std::vector<Eigen::MatrixXcd> &SZ) const override;

    void set_parameters(TableMap &parameters) override;
    void set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) override;
    void set_coupling_damping(double alpha) override;
    void set_field_damping(double beta) override;
    void build_mpo() override;
    void randomize_hamiltonian() override;
    void set_averages(std::vector<TableMap> all_parameters, bool reverse = false) override;
    void write_hamiltonian(h5pp::File &file, const std::string &model_prefix) const override;
    void read_hamiltonian(const h5pp::File &file, const std::string &model_prefix) override;
};
