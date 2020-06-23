#pragma once
#include <complex>
#include <config/enums.h>
#include <list>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>
#include <memory>
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

    std::list<std::unique_ptr<class_env_ene>> eneL;
    std::list<std::unique_ptr<class_env_ene>> eneR;
    std::list<std::unique_ptr<class_env_var>> varL;
    std::list<std::unique_ptr<class_env_var>> varR;

    public:
    class_edges_finite();
    ~class_edges_finite();                                              // Read comment on implementation
    class_edges_finite(class_edges_finite &&other);                     // default move ctor
    class_edges_finite &operator=(class_edges_finite &&other);          // default move assign
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

    void eject_inactive_edges(std::optional<std::vector<size_t>> sites = std::nullopt);
    void eject_all_edges();

    // This is a reference wrapper for an edge pair
    template<typename env_type>
    struct env_pair {
        env_type &L;
        env_type &R;
        void      assert_validity() const;
        //        env_pair(env_type &L_, env_type &R_) : L(L_), R(R_) {}
    };

    [[nodiscard]] env_pair<const class_env_ene> get_ene(size_t posL, size_t posR) const;
    [[nodiscard]] env_pair<const class_env_var> get_var(size_t posL, size_t posR) const;
    [[nodiscard]] env_pair<class_env_ene>       get_ene(size_t posL, size_t posR);
    [[nodiscard]] env_pair<class_env_var>       get_var(size_t posL, size_t posR);

    [[nodiscard]] env_pair<const class_env_ene> get_ene(size_t pos) const;
    [[nodiscard]] env_pair<const class_env_var> get_var(size_t pos) const;
    [[nodiscard]] env_pair<class_env_ene>       get_ene(size_t pos);
    [[nodiscard]] env_pair<class_env_var>       get_var(size_t pos);

    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_ene_blk(size_t posL, size_t posR) const;
    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 4>> get_var_blk(size_t posL, size_t posR) const;
    [[nodiscard]] env_pair<Eigen::Tensor<Scalar, 3>>       get_ene_blk(size_t posL, size_t posR);
    [[nodiscard]] env_pair<Eigen::Tensor<Scalar, 4>>       get_var_blk(size_t posL, size_t posR);

    [[nodiscard]] env_pair<const class_env_ene> get_multisite_ene(std::optional<std::vector<size_t>> sites = std::nullopt) const;
    [[nodiscard]] env_pair<const class_env_var> get_multisite_var(std::optional<std::vector<size_t>> sites = std::nullopt) const;
    [[nodiscard]] env_pair<class_env_ene>       get_multisite_ene(std::optional<std::vector<size_t>> sites = std::nullopt);
    [[nodiscard]] env_pair<class_env_var>       get_multisite_var(std::optional<std::vector<size_t>> sites = std::nullopt);

    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_multisite_ene_blk(std::optional<std::vector<size_t>> sites = std::nullopt) const;
    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 4>> get_multisite_var_blk(std::optional<std::vector<size_t>> sites = std::nullopt) const;
    [[nodiscard]] env_pair<Eigen::Tensor<Scalar, 3>>       get_multisite_ene_blk(std::optional<std::vector<size_t>> sites = std::nullopt);
    [[nodiscard]] env_pair<Eigen::Tensor<Scalar, 4>>       get_multisite_var_blk(std::optional<std::vector<size_t>> sites = std::nullopt);
};
