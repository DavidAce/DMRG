#pragma once

#include <memory>

#include "general/nmspc_tensor_extra.h"
#include <complex>
#include <io/nmspc_logger.h>
#include <optional>

/*! \brief Base environment class for environment blocks och type Left or Right corresponding to a single site.
 */

class class_mps_site;
class class_mpo_base;

class class_env_base {
    public:
    using Scalar = std::complex<double>;

    protected:
    std::optional<size_t> position;
    bool                  edge_has_been_set                                                                 = false;
    virtual void          enlarge(const Eigen::Tensor<Scalar, 3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO) = 0;

    public:
    std::string side;
    size_t      sites = 0; /*!< Number of particles that have been contracted into this environment. */
    explicit class_env_base(std::string side_, size_t position_);
    explicit class_env_base(std::string side_, const class_mps_site &MPS, const class_mpo_base &MPO);

    void         set_position(const size_t position_) { position = position_; }
    size_t       get_position() const {
        if(position) {
            return position.value();
        } else {
            throw std::runtime_error(fmt::format("Position hasn't been set on env side {}", side));
        }
    }

    virtual void clear()                                                             = 0;
    virtual void assert_validity() const                                              = 0;
    virtual void set_edge_dims(const class_mps_site &MPS, const class_mpo_base &MPO) = 0;
    virtual bool is_real() const                                                      = 0;
    virtual bool has_nan() const                                                      = 0;

};
