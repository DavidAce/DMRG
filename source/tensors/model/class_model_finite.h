#pragma once
#include <complex>
#include <config/enums.h>
#include <list>
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>

class class_mpo_site;
class class_tensors_finite;

class class_model_finite {
    public:
    using Scalar = std::complex<double>;

    private:
    friend class_tensors_finite;
    struct Cache {
        std::optional<std::vector<size_t>>      cached_sites          = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 4>> multisite_mpo         = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 2>> multisite_ham         = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 4>> multisite_mpo_squared = std::nullopt;
    };
    mutable Cache                         cache;
    std::vector<Eigen::Tensor<Scalar, 4>> get_compressed_mpo_squared(std::optional<SVDMode> svdMode = std::nullopt);
    void                                  randomize();
    void                                  reset_mpo_squared();
    void                                  rebuild_mpo_squared(std::optional<SVDMode> svdMode = std::nullopt);
    void                                  set_reduced_energy(double total_energy);
    void                                  set_reduced_energy_per_site(double site_energy);
    void                                  perturb_hamiltonian(double coupling_ptb, double field_ptb, PerturbMode perturbMode);
    void                                  damp_model_disorder(double coupling_damp, double field_damp);

    public:
    std::list<std::unique_ptr<class_mpo_site>> MPO; /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */
    std::vector<size_t>                        active_sites;
    ModelType                                  model_type = ModelType::ising_tf_rf;

    public:
    class_model_finite();
    ~class_model_finite();                                          // Read comment on implementation
    class_model_finite(class_model_finite &&other);                 // default move ctor
    class_model_finite &operator=(class_model_finite &&other);      // default move assign
    class_model_finite(const class_model_finite &other);            // copy ctor
    class_model_finite &operator=(const class_model_finite &other); // copy assign

    void                  initialize(ModelType model_type_, size_t model_size);
    size_t                get_length() const;
    bool                  is_real() const;
    bool                  has_nan() const;
    void                  assert_validity() const;
    const class_mpo_site &get_mpo(size_t pos) const;
    class_mpo_site &      get_mpo(size_t pos);

    // For reduced energy MPO's
    [[nodiscard]] bool   is_reduced() const;
    [[nodiscard]] bool   is_perturbed() const;
    [[nodiscard]] bool   is_damped() const;
    [[nodiscard]] double get_energy_reduced() const;
    [[nodiscard]] double get_energy_per_site_reduced() const;

    // For multisite
    Eigen::DSizes<long, 4>          active_dimensions() const;
    Eigen::Tensor<Scalar, 4>        get_multisite_mpo(const std::vector<size_t> &sites, const std::vector<size_t> &nbody = {}) const;
    Eigen::Tensor<Scalar, 2>        get_multisite_ham(const std::vector<size_t> &sites, const std::vector<size_t> & nbody_terms = {}) const;
    const Eigen::Tensor<Scalar, 4> &get_multisite_mpo() const;
    const Eigen::Tensor<Scalar, 2> &get_multisite_ham() const;

    Eigen::DSizes<long, 4>          active_dimensions_squared() const;
    Eigen::Tensor<Scalar, 4>        get_multisite_mpo_squared(const std::vector<size_t> &sites) const;
    const Eigen::Tensor<Scalar, 4> &get_multisite_mpo_squared() const;
    void                            clear_cache() const;
};
