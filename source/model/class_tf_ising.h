//
// Created by david on 2018-07-04.
//

#pragma once

#include "class_model_base.h"
#include <general/nmspc_tensor_extra.h>
#include <h5pp/details/h5ppHid.h>
#include <iomanip>
#include <iostream>

class class_tf_ising : public class_model_base {
    using Scalar = std::complex<double>;

    private:
    h5tb_tf_ising::table pm;

    [[nodiscard]] double get_field() const;
    [[nodiscard]] double get_coupling() const;
    void                 h5table_define();

    public:
    class_tf_ising(size_t position_);
    std::unique_ptr<class_model_base> clone() const override;
    Eigen::Tensor<Scalar, 4>          MPO_reduced_view() const override;
    Eigen::Tensor<Scalar, 4>          MPO_reduced_view(double site_energy) const override;
    Eigen::Tensor<Scalar, 1>          get_MPO_edge_left() const override;
    Eigen::Tensor<Scalar, 1>          get_MPO_edge_right() const override;
    size_t                            get_spin_dimension() const override;
    TableMap                          get_parameters() const override;
    void                              set_parameters(TableMap &parameters) override;
    void                              set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) override;
    void                              set_coupling_damping(double alpha) override;
    void                              set_field_damping(double beta) override;
    void                              build_mpo() override;
    void                              randomize_hamiltonian() override;
    bool                              is_perturbed() const override;
    void                              set_averages(std::vector<TableMap> all_parameters, bool reverse = false) override;
    Eigen::MatrixXcd                  single_site_hamiltonian(size_t position, size_t sites, std::vector<Eigen::MatrixXcd> &SX, std::vector<Eigen::MatrixXcd> &SY,
                                                              std::vector<Eigen::MatrixXcd> &SZ) const override;

    void write_parameters(h5pp::File &file, std::string_view table_name) const override;
    void read_parameters(const h5pp::File &file, std::string_view table_name) override;
};
