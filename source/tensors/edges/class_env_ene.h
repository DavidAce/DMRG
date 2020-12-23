#pragma once
#include <tensors/edges/class_env_base.h>

/*! \brief Environment class with energy MPOs for environment blocks och type Left or Right corresponding to a single site.
 */

class class_env_ene final : public class_env_base {
    public:
    using class_env_base::class_env_base;
    using class_env_base::enlarge;
    using class_env_base::set_edge_dims;
    explicit class_env_ene(std::string side_, const class_mps_site &mps, const class_mpo_site &mpo);
    [[nodiscard]] class_env_ene enlarge(const class_mps_site &mps, const class_mpo_site &mpo) const;
    void refresh(const class_env_ene & env, const class_mps_site &mps, const class_mpo_site &mpo);
    void set_edge_dims(const class_mps_site &mps, const class_mpo_site &mpo) final;
};
