#pragma once
#include <complex>
#include <memory>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>

class class_mpo_base;

class class_model_infinite {
    public:
    using Scalar         = std::complex<double>;
    ModelType model_type = ModelType::ising_tf_rf;

    class_model_infinite() = default;
    ~class_model_infinite(); // Read comment on definition
    class_model_infinite(const class_model_infinite &other);
    class_model_infinite &operator=(const class_model_infinite &other);

    std::unique_ptr<class_mpo_base> HA; /*!< Left hamiltonian MPO */
    std::unique_ptr<class_mpo_base> HB; /*!< Right hamiltonian MPO */

    bool                            is_real() const;
    bool                            has_nan() const;
    void                            assert_validity() const;
    const Eigen::Tensor<Scalar, 4> &get_mpo() const;
    Eigen::DSizes<long, 4>          dimensions() const;

    [[nodiscard]] bool   is_reduced() const;
    [[nodiscard]] double get_energy_reduced() const;
    [[nodiscard]] double get_energy_per_site_reduced() const;

    void set_reduced_energy(double total_energy);
    void set_reduced_energy_per_site(double site_energy);

    void clear_cache();

    private:
    struct Cache {
        std::optional<Eigen::Tensor<Scalar, 4>> mpo = std::nullopt;
    };
    mutable Cache cache;
};
