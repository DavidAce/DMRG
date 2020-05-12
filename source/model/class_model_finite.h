#pragma once
#include <complex>
#include <list>
#include <memory>
#include <simulation/enums.h>
#include <unsupported/Eigen/CXX11/Tensor>
class class_mpo_base;

class class_model_finite {
    public:
    using Scalar = std::complex<double>;

    class_model_finite() = default;
    ~class_model_finite();
    class_model_finite(const class_model_finite &other);
    class_model_finite &operator=(const class_model_finite &other);

    std::list<std::unique_ptr<class_mpo_base>> MPO; /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */

    size_t                get_length() const;
    bool                  isReal() const;
    bool                  hasNaN() const;
    void                  assertValidity() const;
    const class_mpo_base &get_MPO(size_t pos) const;
    class_mpo_base &      get_MPO(size_t pos);

    // For reduced energy MPO's
    bool   isReduced() const;
    double get_energy_reduced() const;
    double get_energy_per_site_reduced() const;
    void   set_reduced_energy(double total_energy);
    void   set_reduced_energy_per_site(double site_energy);
    void   perturb_hamiltonian(double coupling_ptb, double field_ptb, PerturbMode perturbMode);
    void   damp_hamiltonian(double coupling_damp, double field_damp);
    bool   is_perturbed() const;
    bool   is_damped() const;

    // For multisite
    std::list<size_t>               active_sites;
    const Eigen::Tensor<Scalar, 4> &get_multisite_mpo() const;

    private:
    struct Cache {
        std::optional<Eigen::Tensor<Scalar, 4>> multisite_mpo = std::nullopt;
        std::optional<std::list<size_t>>        cached_sites = std::nullopt;
    };

    public:
    mutable Cache cache;
    void          clear_cache() const;
};
