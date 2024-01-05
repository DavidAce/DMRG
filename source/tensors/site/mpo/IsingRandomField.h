#pragma once

#include "MpoSite.h"
#include <h5pp/details/h5ppHid.h>
#include "math/tenx/fwd_decl.h"

class IsingRandomField : public MpoSite {
    private:
    h5tb_ising_tf_rf h5tb;

    [[nodiscard]] double get_field() const;
    [[nodiscard]] double get_coupling() const;

    public:
    IsingRandomField(ModelType model_type_, size_t position_);
    [[nodiscard]] std::unique_ptr<MpoSite> clone() const final;
    [[nodiscard]] Eigen::Tensor<cplx, 4>   MPO_nbody_view(std::optional<std::vector<size_t>> nbody,
                                                          std::optional<std::vector<size_t>> skip = std::nullopt) const final;
    [[nodiscard]] Eigen::Tensor<cplx_t, 4> MPO_nbody_view_t(std::optional<std::vector<size_t>> nbody,
                                                            std::optional<std::vector<size_t>> skip = std::nullopt) const final;

    [[nodiscard]] Eigen::Tensor<cplx, 4> MPO_shifted_view() const final;
    [[nodiscard]] Eigen::Tensor<cplx, 4> MPO_shifted_view(double site_energy) const final;
    [[nodiscard]] long                   get_spin_dimension() const final;
    [[nodiscard]] TableMap               get_parameters() const final;
    [[nodiscard]] std::any               get_parameter(std::string_view name) const final;

    void print_parameter_names() const final;
    void print_parameter_values() const final;
    void set_parameters(TableMap &parameters) final;
    void build_mpo() final;
    void randomize_hamiltonian() final;
    void set_averages(std::vector<TableMap> all_parameters, bool infinite = false) final;

    void save_hamiltonian(h5pp::File &file, std::string_view table_path) const final;
    void load_hamiltonian(const h5pp::File &file, std::string_view model_prefix) final;
};
