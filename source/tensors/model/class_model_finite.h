#pragma once
#include <complex>
#include <config/enums.h>
#include <list>
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>
class class_mpo_site;

class class_model_finite {
    public:
    using Scalar = std::complex<double>;

    private:
    struct Cache {
        std::optional<Eigen::Tensor<Scalar, 4>> multisite_tensor = std::nullopt;
        std::optional<std::vector<size_t>>        cached_sites     = std::nullopt;
    };
    mutable Cache cache;

    public:
    std::list<std::unique_ptr<class_mpo_site>> MPO; /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */
    std::vector<size_t>                        active_sites;
    ModelType                                  model_type = ModelType::ising_tf_rf;

    public:
    class_model_finite();
    ~class_model_finite();                                              // Read comment on implementation
    class_model_finite(class_model_finite &&other);            // default move ctor
    class_model_finite &operator=(class_model_finite &&other); // default move assign
    class_model_finite(const class_model_finite &other);                // copy ctor
    class_model_finite &operator=(const class_model_finite &other);     // copy assign

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

    void randomize();
    void set_reduced_energy(double total_energy);
    void set_reduced_energy_per_site(double site_energy);
    void perturb_hamiltonian(double coupling_ptb, double field_ptb, PerturbMode perturbMode);
    void damp_hamiltonian(double coupling_damp, double field_damp);

    // For multisite
    Eigen::DSizes<long, 4>          active_dimensions() const;
    Eigen::Tensor<Scalar, 4>        get_multisite_tensor(const std::vector<size_t> &sites) const;
    const Eigen::Tensor<Scalar, 4> &get_multisite_tensor() const;

    void clear_cache() const;
};
