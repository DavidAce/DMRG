#pragma once

#include "MpoSite.h"
#include <general/eigen_tensor_fwd_decl.h>
#include <h5pp/details/h5ppHid.h>

class LBit : public MpoSite {
    using Scalar = std::complex<double>;

    private:
    h5tb_lbit h5tb;
    //    [[nodiscard]] double get_field() const;
    //    [[nodiscard]] double get_coupling() const;

    public:
    LBit(ModelType model_type_, size_t position_);
    [[nodiscard]] std::unique_ptr<MpoSite> clone() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 4> MPO_nbody_view(std::optional<std::vector<size_t>> nbody,
                                                          std::optional<std::vector<size_t>> skip = std::nullopt) const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 4> MPO_reduced_view() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 4> MPO_reduced_view(double site_energy) const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 1> get_MPO_edge_left() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 1> get_MPO_edge_right() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 1> get_MPO2_edge_left() const override;
    [[nodiscard]] Eigen::Tensor<Scalar, 1> get_MPO2_edge_right() const override;
    [[nodiscard]] long                     get_spin_dimension() const override;
    [[nodiscard]] TableMap                 get_parameters() const override;
    [[nodiscard]] bool                     is_perturbed() const override;
    //    [[nodiscard]] Eigen::MatrixXcd single_site_hamiltonian(size_t position, size_t sites, std::vector<Eigen::MatrixXcd> &SX, std::vector<Eigen::MatrixXcd>
    //    &SY,
    //                                                           std::vector<Eigen::MatrixXcd> &SZ) const override;
    //
    void print_parameter_names() const override;
    void print_parameter_values() const override;
    void set_parameters(TableMap &parameters) override;
    void set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) override;
    void build_mpo() override;
    void randomize_hamiltonian() override;
    void set_averages(std::vector<TableMap> all_parameters, bool infinite = false, bool reverse = false) override;

    void save_hamiltonian(h5pp::File &file, std::string_view table_path) const override;
    void load_hamiltonian(const h5pp::File &file, std::string_view model_prefix) override;
};
