//
// Created by david on 2018-07-06.
//

#pragma once

#include "class_mpo_site.h"
#include <general/eigen_tensor_fwd_decl.h>
#include <h5pp/details/h5ppHid.h>

namespace h5pp::hid {
    class h5t;
}

class class_ising_sdual : public class_mpo_site {
    using Scalar = std::complex<double>;

    private:
    h5tb_ising_sdual     h5tb;
    [[nodiscard]] double get_coupling() const;
    [[nodiscard]] double get_field() const;
    [[nodiscard]] double get_coupling(double J_rnd_, double J_ptb_, double alpha_) const;
    [[nodiscard]] double get_field(double h_rnd_, double h_ptb_, double beta_) const;

    public:
    explicit class_ising_sdual(ModelType model_type_, size_t position_);

    // Functions that extend the base (no override)
    void set_realization_averages(double J_avrg_, double h_avrg_);
    // Functions that override the base
    [[nodiscard]] std::unique_ptr<class_mpo_site> clone() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 4>        MPO_nbody_view(const std::vector<size_t> &nbody_terms) const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 4>        MPO_reduced_view() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 4>        MPO_reduced_view(double site_energy) const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 1>        get_MPO_edge_left() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 1>        get_MPO_edge_right() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 1>        get_MPO2_edge_left() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 1>        get_MPO2_edge_right() const override;
    [[nodiscard]] long                            get_spin_dimension() const override;
    [[nodiscard]] TableMap                        get_parameters() const override;
    [[nodiscard]] bool                            is_perturbed() const override;

    void print_parameter_names() const override;
    void print_parameter_values() const override;
    void set_parameters(TableMap &parameters) override;
    void set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) override;
    void set_coupling_damping(double alpha) override;
    void set_field_damping(double beta) override;
    void build_mpo() override;
    void randomize_hamiltonian() override;
    void set_averages(std::vector<TableMap> all_parameters, bool infinite = false, bool reverse = false) override;
    void save_hamiltonian(h5pp::File &file, const std::string &hamiltonian_table_path) const override;
    void load_hamiltonian(const h5pp::File &file, const std::string &model_prefix) override;
};
