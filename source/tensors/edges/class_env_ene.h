#pragma once
#include <tensors/edges/class_env_base.h>

/*! \brief Environment class with energy MPOs for environment blocks och type Left or Right corresponding to a single site.
 */

class class_env_ene final : public class_env_base {
    public:
    using class_env_base::class_env_base;
    using class_env_base::enlarge;
    using class_env_base::set_edge_dims;
    explicit class_env_ene(std::string side_, const class_mps_site &MPS, const class_mpo_site &MPO);
    [[nodiscard]] class_env_ene enlarge(const class_mps_site &MPS, const class_mpo_site &MPO);
    void set_edge_dims(const class_mps_site &MPS, const class_mpo_site &MPO) final;
};
