#pragma once
#include <array>
#include <complex>
#include <config/enums.h>
#include <math/svd/settings.h>
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>
class MpoSite;
class TensorsFinite;

class ModelFinite {
    public:
    using Scalar = std::complex<double>;

    private:
    friend TensorsFinite;
    struct Cache {
        std::optional<std::vector<size_t>>      cached_sites          = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 4>> multisite_mpo         = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 2>> multisite_ham         = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 4>> multisite_mpo_squared = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 2>> multisite_ham_squared = std::nullopt;
    };
    mutable Cache                         cache;
    std::vector<Eigen::Tensor<Scalar, 4>> get_compressed_mpo_squared(std::optional<svd::settings> svd_settings = std::nullopt);
    void                                  randomize();
    bool                                  has_mpo_squared() const;
    void                                  reset_mpo_squared();
    void                                  clear_mpo_squared();
    void                                  compress_mpo_squared(std::optional<svd::settings> svd_settings = std::nullopt);
    void                                  set_reduced_energy(double total_energy);
    void                                  set_reduced_energy_per_site(double site_energy);
    void                                  perturb_hamiltonian(double coupling_ptb, double field_ptb, PerturbMode perturbMode);

    public:
    std::vector<std::unique_ptr<MpoSite>> MPO; /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */
    std::vector<size_t>                   active_sites;
    ModelType                             model_type = ModelType::ising_tf_rf;

    public:
    ModelFinite();
    ~ModelFinite();                                   // Read comment on implementation
    ModelFinite(ModelFinite &&other);                 // default move ctor
    ModelFinite &operator=(ModelFinite &&other);      // default move assign
    ModelFinite(const ModelFinite &other);            // copy ctor
    ModelFinite &operator=(const ModelFinite &other); // copy assign

    void           initialize(ModelType model_type_, size_t model_size);
    size_t         get_length() const;
    bool           is_real() const;
    bool           has_nan() const;
    void           assert_validity() const;
    const MpoSite &get_mpo(size_t pos) const;
    MpoSite       &get_mpo(size_t pos);

    // For reduced energy MPO's
    [[nodiscard]] bool   is_reduced() const;
    [[nodiscard]] bool   is_perturbed() const;
    [[nodiscard]] bool   is_compressed_mpo_squared() const;
    [[nodiscard]] double get_energy_reduced() const;
    [[nodiscard]] double get_energy_per_site_reduced() const;

    // For multisite
    std::array<long, 4>             active_dimensions() const;
    Eigen::Tensor<Scalar, 4>        get_multisite_mpo(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody = std::nullopt) const;
    Eigen::Tensor<Scalar, 2>        get_multisite_ham(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody = std::nullopt) const;
    const Eigen::Tensor<Scalar, 4> &get_multisite_mpo() const;
    const Eigen::Tensor<Scalar, 2> &get_multisite_ham() const;

    Eigen::Tensor<Scalar, 4> get_multisite_mpo_reduced_view(double energy_per_site) const;
    Eigen::Tensor<Scalar, 4> get_multisite_mpo_squared_reduced_view(double energy_per_site) const;

    std::array<long, 4>             active_dimensions_squared() const;
    Eigen::Tensor<Scalar, 4>        get_multisite_mpo_squared(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody = std::nullopt) const;
    Eigen::Tensor<Scalar, 2>        get_multisite_ham_squared(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody = std::nullopt) const;
    const Eigen::Tensor<Scalar, 4> &get_multisite_mpo_squared() const;
    const Eigen::Tensor<Scalar, 2> &get_multisite_ham_squared() const;
    void                            clear_cache(LogPolicy logPolicy = LogPolicy::QUIET) const;

    std::vector<size_t> get_active_ids() const;
};
