#pragma once
#include <complex>
#include <config/enums.h>
#include <general/eigen_tensor_fwd_decl.h>
#include <memory>
#include <optional>
#include <vector>
#include "class_env_pair.h"

class class_env_ene;
class class_env_var;

class class_edges_finite {
    public:
    using Scalar = std::complex<double>;
    std::vector<size_t> active_sites;

    private:
    //    size_t iter      = 0;
    //    size_t step      = 0;
    //    int    direction = 1;

    std::vector<std::unique_ptr<class_env_ene>> eneL;
    std::vector<std::unique_ptr<class_env_ene>> eneR;
    std::vector<std::unique_ptr<class_env_var>> varL;
    std::vector<std::unique_ptr<class_env_var>> varR;

    public:
    class_edges_finite();
    ~class_edges_finite();                                              // Read comment on implementation
    class_edges_finite(class_edges_finite &&other) noexcept;            // default move ctor
    class_edges_finite &operator=(class_edges_finite &&other) noexcept; // default move assign
    class_edges_finite(const class_edges_finite &other);                // copy ctor
    class_edges_finite &operator=(const class_edges_finite &other);     // copy assign

    void initialize(size_t model_size);

    [[nodiscard]] const class_env_ene &get_eneL(size_t pos) const;
    [[nodiscard]] const class_env_ene &get_eneR(size_t pos) const;
    [[nodiscard]] const class_env_var &get_varL(size_t pos) const;
    [[nodiscard]] const class_env_var &get_varR(size_t pos) const;
    [[nodiscard]] class_env_ene &      get_eneL(size_t pos);
    [[nodiscard]] class_env_var &      get_varR(size_t pos);
    [[nodiscard]] class_env_var &      get_varL(size_t pos);
    [[nodiscard]] class_env_ene &      get_eneR(size_t pos);

    [[nodiscard]] size_t get_length() const;
    [[nodiscard]] bool   is_real() const;
    [[nodiscard]] bool   has_nan() const;
    void                 assert_validity() const;

    void eject_edges_inactive_ene(std::optional<std::vector<size_t>> sites = std::nullopt);
    void eject_edges_inactive_var(std::optional<std::vector<size_t>> sites = std::nullopt);
    void eject_edges_inactive(std::optional<std::vector<size_t>> sites = std::nullopt);

    void eject_edges_all_ene();
    void eject_edges_all_var();
    void eject_edges_all();

    [[nodiscard]] env_pair<const class_env_ene> get_ene(size_t posL, size_t posR) const;
    [[nodiscard]] env_pair<const class_env_var> get_var(size_t posL, size_t posR) const;
    [[nodiscard]] env_pair<class_env_ene>       get_ene(size_t posL, size_t posR);
    [[nodiscard]] env_pair<class_env_var>       get_var(size_t posL, size_t posR);

    [[nodiscard]] env_pair<const class_env_ene> get_ene(size_t pos) const;
    [[nodiscard]] env_pair<const class_env_var> get_var(size_t pos) const;
    [[nodiscard]] env_pair<class_env_ene>       get_ene(size_t pos);
    [[nodiscard]] env_pair<class_env_var>       get_var(size_t pos);

    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_ene_blk(size_t posL, size_t posR) const;
    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_var_blk(size_t posL, size_t posR) const;
    [[nodiscard]] env_pair<Eigen::Tensor<Scalar, 3>>       get_ene_blk(size_t posL, size_t posR);
    [[nodiscard]] env_pair<Eigen::Tensor<Scalar, 3>>       get_var_blk(size_t posL, size_t posR);

    [[nodiscard]] env_pair<const class_env_ene> get_multisite_ene(std::optional<std::vector<size_t>> sites = std::nullopt) const;
    [[nodiscard]] env_pair<const class_env_var> get_multisite_var(std::optional<std::vector<size_t>> sites = std::nullopt) const;
    [[nodiscard]] env_pair<class_env_ene>       get_multisite_ene(std::optional<std::vector<size_t>> sites = std::nullopt);
    [[nodiscard]] env_pair<class_env_var>       get_multisite_var(std::optional<std::vector<size_t>> sites = std::nullopt);

    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_multisite_ene_blk(std::optional<std::vector<size_t>> sites = std::nullopt) const;
    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_multisite_var_blk(std::optional<std::vector<size_t>> sites = std::nullopt) const;
    [[nodiscard]] env_pair<Eigen::Tensor<Scalar, 3>>       get_multisite_ene_blk(std::optional<std::vector<size_t>> sites = std::nullopt);
    [[nodiscard]] env_pair<Eigen::Tensor<Scalar, 3>>       get_multisite_var_blk(std::optional<std::vector<size_t>> sites = std::nullopt);

    [[nodiscard]] std::pair<std::vector<size_t>, std::vector<size_t>> get_active_ids() const;
};
