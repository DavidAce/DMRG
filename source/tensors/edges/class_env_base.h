#pragma once

#include <complex>
#include <general/eigen_tensor_fwd_decl.h>
#include <memory>
#include <optional>

/*! \brief Base environment class for environment blocks och type Left or Right corresponding to a single site.
 */

class class_env_ene;
class class_env_var;
class class_mps_site;
class class_mpo_site;

class class_env_base {
    public:
    using Scalar = std::complex<double>;

    protected:
    bool edge_has_been_set = false;
    void build_block(Eigen::Tensor<Scalar, 3> &otherblock, const Eigen::Tensor<Scalar, 3> &mps, const Eigen::Tensor<Scalar, 4> &mpo);
    void enlarge(const Eigen::Tensor<Scalar, 3> &mps, const Eigen::Tensor<Scalar, 4> &mpo);
    void set_edge_dims(const Eigen::Tensor<Scalar, 3> &mps, const Eigen::Tensor<Scalar, 4> &mpo, const Eigen::Tensor<Scalar, 1> &edge);

    std::unique_ptr<Eigen::Tensor<Scalar, 3>> block;        /*!< The environment block. */
    size_t                                    sites    = 0; /*!< Number of particles that have been contracted into this environment. */
    std::optional<size_t>                     position = std::nullopt;
    std::string                               side;
    std::string                               tag;
    mutable std::optional<std::size_t>        unique_id;
    mutable std::optional<std::size_t>        unique_id_mps; // Unique identifiers of the neighboring site which are used to build this block
    mutable std::optional<std::size_t>        unique_id_mpo; // Unique identifiers of the neighboring site which are used to build this block
    mutable std::optional<std::size_t>        unique_id_env; // Unique identifiers of the neighboring site which are used to build this block

    public:
    class_env_base();
    ~class_env_base();                                          // Read comment on implementation
    class_env_base(class_env_base &&other) noexcept;            // default move ctor
    class_env_base &operator=(class_env_base &&other) noexcept; // default move assign
    class_env_base(const class_env_base &other);                // copy ctor
    class_env_base &operator=(const class_env_base &other);     // copy assign

    explicit class_env_base(std::string side_, size_t position_);
    explicit class_env_base(std::string side_, const class_mps_site &MPS, const class_mpo_site &MPO);

    void clear();

    void set_position(const size_t position_) { position = position_; }
    void assert_block() const;
    void assert_validity() const;

    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &get_block() const;
    [[nodiscard]] Eigen::Tensor<Scalar, 3> &      get_block();
    [[nodiscard]] bool                            has_block() const;
    [[nodiscard]] bool                            is_real() const;
    [[nodiscard]] bool                            has_nan() const;
    [[nodiscard]] size_t                          get_position() const;
    [[nodiscard]] size_t                          get_sites() const;

    virtual void set_edge_dims(const class_mps_site &MPS, const class_mpo_site &MPO) = 0;

    std::size_t get_unique_id() const;
    std::optional<std::size_t> get_unique_id_env() const;
    std::optional<std::size_t> get_unique_id_mps() const;
    std::optional<std::size_t> get_unique_id_mpo() const;
};
