#pragma once
#include <complex>
#include <config/enums.h>
#include <memory>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>

class class_mpo_site;

class class_model_infinite {
    public:
    using Scalar = std::complex<double>;

    private:
    struct Cache {
        std::optional<Eigen::Tensor<Scalar, 4>> twosite_mpo_AB = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 4>> twosite_mpo_BA = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 2>> twosite_ham_AB = std::nullopt;
        std::optional<Eigen::Tensor<Scalar, 2>> twosite_ham_BA = std::nullopt;
    };
    mutable Cache                   cache;
    std::unique_ptr<class_mpo_site> HA; /*!< Left hamiltonian MPO */
    std::unique_ptr<class_mpo_site> HB; /*!< Right hamiltonian MPO */

    public:
    ModelType model_type = ModelType::ising_tf_rf;

    class_model_infinite();
    ~class_model_infinite();                                            // Read comment on implementation
    class_model_infinite(class_model_infinite &&other);                 // default move ctor
    class_model_infinite &operator=(class_model_infinite &&other);      // default move assign
    class_model_infinite(const class_model_infinite &other);            // copy ctor
    class_model_infinite &operator=(const class_model_infinite &other); // copy assign

    void                                                        initialize(ModelType model_type_);
    void                                                        randomize();
    void                                                        reset_mpo_squared();
    void                                                        rebuild_mpo_squared(std::optional<SVDMode> svdMode = std::nullopt);
    std::vector<Eigen::Tensor<class_model_infinite::Scalar, 4>> get_compressed_mpo_squared(std::optional<SVDMode> svdMode = std::nullopt);

    bool is_real() const;
    bool has_nan() const;
    void assert_validity() const;

    [[nodiscard]] const class_mpo_site &get_mpo_siteA() const;
    [[nodiscard]] const class_mpo_site &get_mpo_siteB() const;
    [[nodiscard]] class_mpo_site &      get_mpo_siteA();
    [[nodiscard]] class_mpo_site &      get_mpo_siteB();
    const Eigen::Tensor<Scalar, 4> &    get_2site_mpo_AB() const;
    const Eigen::Tensor<Scalar, 4> &    get_2site_mpo_BA() const;
    const Eigen::Tensor<Scalar, 2> &    get_2site_ham_AB() const;
    const Eigen::Tensor<Scalar, 2> &    get_2site_ham_BA() const;

    Eigen::DSizes<long, 4> dimensions() const;

    [[nodiscard]] bool   is_reduced() const;
    [[nodiscard]] double get_energy_per_site_reduced() const;
    void                 set_reduced_energy_per_site(double site_energy);
    void                 clear_cache();
};
