#pragma once

#include "math/tenx/fwd_decl.h"
#include "MpoSite.h"
#include <h5pp/details/h5ppHid.h>

namespace h5pp::hid {
    class h5t;
}

class XXZ : public MpoSite {
    private:
    h5tb_xxz               h5tb;
    [[nodiscard]] double   get_coupling() const;
    [[nodiscard]] double   get_field() const;
    Eigen::Tensor<cx64, 4> get_mpo(cx64 energy_shift_per_site, std::optional<std::vector<size_t>> nbody = std::nullopt,
                                   std::optional<std::vector<size_t>> skip = std::nullopt) const final;

    public:
    explicit XXZ(ModelType model_type_, size_t position_);

    // Functions that extend the base (no override)
    void set_realization_averages(double J_avrg_, double h_avrg_);
    // Functions that override the base
    [[nodiscard]] std::unique_ptr<MpoSite> clone() const override;
    [[nodiscard]] long                     get_spin_dimension() const final;
    [[nodiscard]] TableMap                 get_parameters() const final;
    [[nodiscard]] std::any                 get_parameter(std::string_view name) const override;
    void                                   set_parameter(std::string_view name, std::any value) final;
    void                                   print_parameter_names() const final;
    void                                   print_parameter_values() const final;
    void                                   set_parameters(TableMap &parameters) final;
    void                                   randomize_hamiltonian() final;
    void                                   set_averages(std::vector<TableMap> all_parameters, bool infinite = false) final;
    void                                   save_hamiltonian(h5pp::File &file, std::string_view hamiltonian_table_path) const final;
    void                                   load_hamiltonian(const h5pp::File &file, std::string_view model_prefix) final;
};
