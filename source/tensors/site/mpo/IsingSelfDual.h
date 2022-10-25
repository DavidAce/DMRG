#pragma once

#include "MpoSite.h"
#include <h5pp/details/h5ppHid.h>
#include <math/tenx/fwd_decl.h>

namespace h5pp::hid {
    class h5t;
}

class IsingSelfDual : public MpoSite {
    private:
    h5tb_ising_selfdual  h5tb;
    [[nodiscard]] double get_coupling() const;
    [[nodiscard]] double get_field() const;

    public:
    explicit IsingSelfDual(ModelType model_type_, size_t position_);

    // Functions that extend the base (no final)
    void set_realization_averages(double J_avrg_, double h_avrg_);
    // Functions that final the base
    [[nodiscard]] std::unique_ptr<MpoSite> clone() const final;
    [[nodiscard]] Eigen::Tensor<cplx, 4>   MPO_nbody_view(std::optional<std::vector<size_t>> nbody,
                                                          std::optional<std::vector<size_t>> skip = std::nullopt) const final;
    [[nodiscard]] Eigen::Tensor<cplx, 4>   MPO_shifted_view() const final;
    [[nodiscard]] Eigen::Tensor<cplx, 4>   MPO_shifted_view(double energy_shift_per_site) const final;
    [[nodiscard]] long                     get_spin_dimension() const final;
    [[nodiscard]] TableMap                 get_parameters() const final;
    [[nodiscard]] std::any                 get_parameter(const std::string &name) const final;

    void print_parameter_names() const final;
    void print_parameter_values() const final;
    void set_parameters(TableMap &parameters) final;
    void build_mpo() final;
    void randomize_hamiltonian() final;
    void set_averages(std::vector<TableMap> all_parameters, bool infinite = false) final;
    void save_hamiltonian(h5pp::File &file, std::string_view hamiltonian_table_path) const final;
    void load_hamiltonian(const h5pp::File &file, std::string_view model_prefix) final;
};
