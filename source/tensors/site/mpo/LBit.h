#pragma once

#include "math/tenx/fwd_decl.h"
#include "MpoSite.h"
#include <h5pp/details/h5ppHid.h>

class LBit : public MpoSite {
    private:
    h5tb_lbit h5tb;
    Eigen::Tensor<cplx, 4> get_mpo(cplx energy_shift_per_site, std::optional<std::vector<size_t>> nbody = std::nullopt,
                                   std::optional<std::vector<size_t>> skip = std::nullopt) const final;
    Eigen::Tensor<cplx_t, 4> get_mpo_t(cplx_t energy_shift_per_site, std::optional<std::vector<size_t>> nbody = std::nullopt,
                                 std::optional<std::vector<size_t>> skip = std::nullopt) const final;
    public:
                                           LBit(ModelType model_type_, size_t position_);
    [[nodiscard]] std::unique_ptr<MpoSite> clone() const final;
    [[nodiscard]] long                     get_spin_dimension() const final;
    [[nodiscard]] TableMap                 get_parameters() const final;
    [[nodiscard]] std::any                 get_parameter(std::string_view name) const final;
    //    [[nodiscard]] Eigen::MatrixXcd single_site_hamiltonian(size_t position, size_t sites, std::vector<Eigen::MatrixXcd> &SX, std::vector<Eigen::MatrixXcd>
    //    &SY,
    //                                                           std::vector<Eigen::MatrixXcd> &SZ) const final;
    //
    void print_parameter_names() const final;
    void print_parameter_values() const final;
    void set_parameters(TableMap &parameters) final;
    void randomize_hamiltonian() final;
    void set_averages(std::vector<TableMap> all_parameters, bool infinite = false) final;

    void save_hamiltonian(h5pp::File &file, std::string_view table_path) const final;
    void load_hamiltonian(const h5pp::File &file, std::string_view model_prefix) final;
};
