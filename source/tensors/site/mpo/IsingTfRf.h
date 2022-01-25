#pragma once

#include "MpoSite.h"
#include <h5pp/details/h5ppHid.h>
#include <math/tenx/fwd_decl.h>

class IsingTfRf : public MpoSite {
    private:
    h5tb_ising_tf_rf h5tb;

    [[nodiscard]] double get_field() const;
    [[nodiscard]] double get_coupling() const;

    public:
    IsingTfRf(ModelType model_type_, size_t position_);
    [[nodiscard]] std::unique_ptr<MpoSite> clone() const override;
    [[nodiscard]] Eigen::Tensor<cplx, 4>   MPO_nbody_view(std::optional<std::vector<size_t>> nbody,
                                                          std::optional<std::vector<size_t>> skip = std::nullopt) const override;
    [[nodiscard]] Eigen::Tensor<cplx, 4>   MPO_shifted_view() const override;
    [[nodiscard]] Eigen::Tensor<cplx, 4>   MPO_shifted_view(double site_energy) const override;
    [[nodiscard]] long                     get_spin_dimension() const override;
    [[nodiscard]] TableMap                 get_parameters() const override;
    [[nodiscard]] bool                     is_perturbed() const override;

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
