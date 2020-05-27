#pragma once
#include <tensors/edges/class_env_base.h>

/*! \brief Environment class with variance MPOs (i.e. double layer of energy MPOs) for environment blocks och type Left or Right corresponding to a single site.
 */

class class_env_var final : public class_env_base {
    private:
    void enlarge(const Eigen::Tensor<Scalar, 3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO) final;

    public:
    Eigen::Tensor<Scalar, 4> block; /*!< The environment block. */
    using class_env_base::class_env_base;
    explicit class_env_var(std::string side_, const class_mps_site &MPS, const class_mpo_site &MPO);
    [[nodiscard]] class_env_var enlarge(const class_mps_site &MPS, const class_mpo_site &MPO);

    void clear() final;
    bool has_block() const final;
    void assert_validity() const final;
    void set_edge_dims(const class_mps_site &MPS, const class_mpo_site &MPO) final;

    [[nodiscard]] bool is_real() const final;
    [[nodiscard]] bool has_nan() const final;
};
