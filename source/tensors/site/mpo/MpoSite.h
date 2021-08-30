#pragma once

#include "MpoParameters.h"
#include <any>
#include <config/enums.h>
#include <map>
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>

namespace h5pp {
    class File;
}

class MpoSite {
    public:
    using cplx     = std::complex<double>;
    using TableMap = std::map<std::string, std::any>;
    const ModelType model_type;
    bool            all_mpo_parameters_have_been_set = false;

    protected:
    mutable std::optional<std::size_t> unique_id;
    mutable std::optional<std::size_t> unique_id_sq;
    // Common parameters
    double e_reduced  = 0;     /*!< "Reduced" energy offset for this mpo (to make energy-reduced MPO views) */
    double psfactor   = 0;     /*!< Parity sector separation factor */
    bool   parity_sep = false; /*!< Parity sector separation on/off */

    std::array<long, 4>                   extent4{}; /*!< Extent of pauli matrices in a rank-4 tensor */
    std::array<long, 2>                   extent2{}; /*!< Extent of pauli matrices in a rank-2 tensor */
    std::optional<size_t>                 position;  /*!< Position on a finite chain */
    Eigen::Tensor<cplx, 4>                mpo_internal;
    std::optional<Eigen::Tensor<cplx, 4>> mpo_squared = std::nullopt;

    public:
    explicit MpoSite(ModelType model_type_, size_t position_);
    virtual ~MpoSite() = default;

    void                                        set_position(size_t new_pos);
    void                                        assert_validity() const;
    void                                        set_reduced_energy(double site_energy);
    void                                        build_mpo_squared();
    void                                        set_mpo_squared(const Eigen::Tensor<cplx, 4> &mpo_sq);
    void                                        clear_mpo_squared();
    [[nodiscard]] bool                          has_mpo_squared() const;
    [[nodiscard]] Eigen::Tensor<cplx, 4>        get_non_compressed_mpo_squared() const;
    [[nodiscard]] const Eigen::Tensor<cplx, 4> &MPO() const;
    [[nodiscard]] const Eigen::Tensor<cplx, 4> &MPO2() const;
    [[nodiscard]] Eigen::Tensor<cplx, 4>       &MPO2();
    [[nodiscard]] Eigen::Tensor<cplx, 4>        MPO2_nbody_view(const std::vector<size_t> &nbody_terms) const;
    [[nodiscard]] size_t                        get_position() const;
    [[nodiscard]] std::vector<std::string>      get_parameter_names() const;
    [[nodiscard]] std::vector<std::any>         get_parameter_values() const;
    [[nodiscard]] bool                          is_real() const;
    [[nodiscard]] bool                          has_nan() const;
    [[nodiscard]] bool                          is_reduced() const;
    [[nodiscard]] bool                          is_compressed_mpo_squared() const;
    [[nodiscard]] double                        get_reduced_energy() const;

    [[nodiscard]] virtual std::unique_ptr<MpoSite> clone() const                                                = 0;
    [[nodiscard]] virtual Eigen::Tensor<cplx, 4>   MPO_nbody_view(const std::vector<size_t> &nbody_terms) const = 0;
    [[nodiscard]] virtual Eigen::Tensor<cplx, 4>   MPO_reduced_view() const                                     = 0;
    [[nodiscard]] virtual Eigen::Tensor<cplx, 4>   MPO_reduced_view(double single_site_energy) const            = 0;
    [[nodiscard]] virtual Eigen::Tensor<cplx, 1>   get_MPO_edge_left() const                                    = 0;
    [[nodiscard]] virtual Eigen::Tensor<cplx, 1>   get_MPO_edge_right() const                                   = 0;
    [[nodiscard]] virtual Eigen::Tensor<cplx, 1>   get_MPO2_edge_left() const                                   = 0;
    [[nodiscard]] virtual Eigen::Tensor<cplx, 1>   get_MPO2_edge_right() const                                  = 0;
    [[nodiscard]] virtual long                     get_spin_dimension() const                                   = 0;
    [[nodiscard]] virtual TableMap                 get_parameters() const                                       = 0;
    [[nodiscard]] virtual bool                     is_perturbed() const                                         = 0;

    virtual void print_parameter_names() const                                                                   = 0;
    virtual void print_parameter_values() const                                                                  = 0;
    virtual void set_parameters(TableMap &parameters)                                                            = 0;
    virtual void set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode)                    = 0;
    virtual void build_mpo()                                                                                     = 0;
    virtual void randomize_hamiltonian()                                                                         = 0;
    virtual void set_averages(std::vector<TableMap> all_parameters, bool infinite = false, bool reverse = false) = 0;
    virtual void save_hamiltonian(h5pp::File &file, std::string_view model_prefix) const                         = 0;
    virtual void load_hamiltonian(const h5pp::File &file, std::string_view model_prefix)                         = 0;
    void         save_mpo(h5pp::File &file, std::string_view model_prefix) const;
    void         load_mpo(const h5pp::File &file, std::string_view model_prefix);
    std::size_t  get_unique_id() const;
    std::size_t  get_unique_id_sq() const;
};
