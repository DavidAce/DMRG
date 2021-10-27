#pragma once

#include <complex>
#include <math/tenx/fwd_decl.h>
#include <memory>
#include <optional>

/*! \brief Base environment class for environment blocks och type Left or Right corresponding to a single site.
 */

class EnvEne;
class EnvVar;
class MpsSite;
class MpoSite;

class EnvBase {
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
    EnvBase();
    ~EnvBase();                                   // Read comment on implementation
    EnvBase(EnvBase &&other) noexcept;            // default move ctor
    EnvBase &operator=(EnvBase &&other) noexcept; // default move assign
    EnvBase(const EnvBase &other);                // copy ctor
    EnvBase &operator=(const EnvBase &other);     // copy assign

    explicit EnvBase(size_t position_, std::string side_, std::string tag_);
    explicit EnvBase(std::string side_, std::string tag_, const MpsSite &MPS, const MpoSite &MPO);

    void clear();

    void set_position(const size_t position_) { position = position_; }
    void assert_block() const;
    void assert_validity() const;
    void assert_unique_id(const EnvBase &env, const MpsSite &mps, const MpoSite &mpo) const;

    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &get_block() const;
    [[nodiscard]] Eigen::Tensor<Scalar, 3>       &get_block();
    [[nodiscard]] bool                            has_block() const;
    [[nodiscard]] bool                            is_real() const;
    [[nodiscard]] bool                            has_nan() const;
    [[nodiscard]] size_t                          get_position() const;
    [[nodiscard]] size_t                          get_sites() const;

    virtual void set_edge_dims(const MpsSite &MPS, const MpoSite &MPO) = 0;

    std::size_t                get_unique_id() const;
    std::optional<std::size_t> get_unique_id_env() const;
    std::optional<std::size_t> get_unique_id_mps() const;
    std::optional<std::size_t> get_unique_id_mpo() const;

    Eigen::Tensor<Scalar, 3> get_expansion_term(const MpsSite &mps, const MpoSite &mpo, double alpha = 0) const;
};
