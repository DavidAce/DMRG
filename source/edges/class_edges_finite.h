#pragma once
#include <complex>
#include <list>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>
class class_env_ene;
class class_env_var;

class class_edges_finite {
    public:
    using Scalar = std::complex<double>;

    private:
    size_t iter      = 0;
    size_t step      = 0;
    int    direction = 1;

    public:
    template<typename env_type>
    struct env_pair {
        env_type L, R;
        //        size_t get_position()const;
        void assertValidity() const;
        env_pair(const env_type &L_, const env_type &R_) : L(L_), R(R_) {}
    };

    std::list<env_pair<class_env_ene>> env_ene;
    std::list<env_pair<class_env_var>> env_var;
    std::list<size_t>                  active_sites;

    void                 assertValidity() const;
    [[nodiscard]] size_t get_length() const;
    //    [[nodiscard]] size_t get_position() const;
    [[nodiscard]] bool isReal() const;
    [[nodiscard]] bool hasNaN() const;

    template<typename env_type>
    using env_type_ref = std::reference_wrapper<const env_type>;
    template<typename env_type>
    using env_pair_ref = env_pair<env_type_ref<env_type>>;
    [[nodiscard]] env_pair_ref<class_env_ene> get_env_ene(size_t posL, size_t posR) const;
    [[nodiscard]] env_pair_ref<class_env_var> get_env_var(size_t posL, size_t posR) const;
    [[nodiscard]] env_pair_ref<class_env_ene> get_env_ene(std::optional <std::list<size_t>> sites = std::nullopt) const;
    [[nodiscard]] env_pair_ref<class_env_var> get_env_var(std::optional <std::list<size_t>> sites = std::nullopt) const;

    using ene_blk_ref = std::reference_wrapper<const Eigen::Tensor<Scalar,3>>;
    using var_blk_ref = std::reference_wrapper<const Eigen::Tensor<Scalar,4>>;
    [[nodiscard]] std::pair<ene_blk_ref,ene_blk_ref> get_multisite_ene(std::optional <std::list<size_t>> sites = std::nullopt) const;
    [[nodiscard]] std::pair<var_blk_ref,var_blk_ref> get_multisite_var(std::optional <std::list<size_t>> sites = std::nullopt) const;

};
