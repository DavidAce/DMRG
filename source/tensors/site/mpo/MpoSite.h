#pragma once

#include "config/enums.h"
#include "h5tb/h5tb.h"
#include <any>
#include <map>
#include <memory>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>

namespace h5pp {
    class File;
}

class MpoSite {
    public:
    using cplx     = std::complex<double>;
    using real     = double;
    using TableMap = std::map<std::string, std::any>;
    const ModelType model_type;
    bool            all_mpo_parameters_have_been_set = false;

    protected:
    mutable std::optional<std::size_t> unique_id;
    mutable std::optional<std::size_t> unique_id_sq;
    // Common parameters
    std::optional<size_t> position;                  /*!< Position on a finite chain */
    double                e_shift           = 0;     /*!< "Shifted" energy offset for this mpo (to make energy-shifted MPO views) */
    double                psfactor          = 0;     /*!< Parity sector separation factor */
    bool                  parity_sep        = false; /*!< Parity sector separation on/off */
    int                   parity_shift_sign = 0;     /*!< Sign of parity sector to shift */
    std::string           parity_shift_axus = {};    /*!< Unsigned axis {x,y,z} of spin parity sector to shift */

    std::array<long, 4>                   extent4{}; /*!< Extent of pauli matrices in a rank-4 tensor */
    std::array<long, 2>                   extent2{}; /*!< Extent of pauli matrices in a rank-2 tensor */
    Eigen::Tensor<cplx, 4>                mpo_internal;
    std::optional<Eigen::Tensor<cplx, 4>> mpo_squared = std::nullopt;

    public:
    explicit MpoSite(ModelType model_type_, size_t position_);
    virtual ~MpoSite() = default;

    void                                        set_position(size_t new_pos);
    void                                        assert_validity() const;
    void                                        set_energy_shift(double site_energy);
    void                                        set_psfactor(double psfactor);
    void                                        set_parity_shift_mpo_squared(int sign, std::string_view axis);
    std::pair<int, std::string_view>            get_parity_shift_mpo_squared() const;
    void                                        build_mpo_squared();
    void                                        set_mpo(const Eigen::Tensor<cplx, 4> &mpo_sq);
    void                                        set_mpo_squared(const Eigen::Tensor<cplx, 4> &mpo_sq);
    void                                        clear_mpo_squared();
    [[nodiscard]] bool                          has_mpo_squared() const;
    [[nodiscard]] Eigen::Tensor<cplx, 4>        get_non_compressed_mpo_squared() const;
    [[nodiscard]] const Eigen::Tensor<cplx, 4> &MPO() const;
    [[nodiscard]] const Eigen::Tensor<cplx, 4> &MPO2() const;
    [[nodiscard]] Eigen::Tensor<cplx, 4>        MPO2_nbody_view(std::optional<std::vector<size_t>> nbody,
                                                                std::optional<std::vector<size_t>> skip = std::nullopt) const;
    [[nodiscard]] size_t                        get_position() const;
    [[nodiscard]] std::vector<std::string>      get_parameter_names() const;
    [[nodiscard]] std::vector<std::any>         get_parameter_values() const;
    [[nodiscard]] bool                          is_real() const;
    [[nodiscard]] bool                          has_nan() const;
    [[nodiscard]] bool                          is_shifted() const;
    [[nodiscard]] bool                          is_compressed_mpo_squared() const;
    [[nodiscard]] double                        get_energy_shift() const;
    [[nodiscard]] Eigen::Tensor<cplx, 1>        get_MPO_edge_left() const;
    [[nodiscard]] Eigen::Tensor<cplx, 1>        get_MPO_edge_right() const;
    [[nodiscard]] Eigen::Tensor<cplx, 1>        get_MPO2_edge_left() const;
    [[nodiscard]] Eigen::Tensor<cplx, 1>        get_MPO2_edge_right() const;

    [[nodiscard]] virtual std::unique_ptr<MpoSite> clone() const                                                                = 0;
    [[nodiscard]] virtual Eigen::Tensor<cplx, 4>   MPO_nbody_view(std::optional<std::vector<size_t>> nbody,
                                                                  std::optional<std::vector<size_t>> skip = std::nullopt) const = 0;
    [[nodiscard]] virtual Eigen::Tensor<cplx, 4>   MPO_shifted_view() const                                                     = 0;
    [[nodiscard]] virtual Eigen::Tensor<cplx, 4>   MPO_shifted_view(double energy_shift_per_site) const                         = 0;
    [[nodiscard]] virtual long                     get_spin_dimension() const                                                   = 0;
    [[nodiscard]] virtual TableMap                 get_parameters() const                                                       = 0;
    [[nodiscard]] virtual std::any                 get_parameter(const std::string &name) const                                 = 0;

    virtual void print_parameter_names() const                                             = 0;
    virtual void print_parameter_values() const                                            = 0;
    virtual void set_parameters(TableMap &parameters)                                      = 0;
    virtual void build_mpo()                                                               = 0;
    virtual void randomize_hamiltonian()                                                   = 0;
    virtual void set_averages(std::vector<TableMap> all_parameters, bool infinite = false) = 0;
    virtual void save_hamiltonian(h5pp::File &file, std::string_view model_prefix) const   = 0;
    virtual void load_hamiltonian(const h5pp::File &file, std::string_view model_prefix)   = 0;
    void         save_mpo(h5pp::File &file, std::string_view model_prefix) const;
    void         load_mpo(const h5pp::File &file, std::string_view model_prefix);
    std::size_t  get_unique_id() const;
    std::size_t  get_unique_id_sq() const;
};
