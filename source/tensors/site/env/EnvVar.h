#pragma once
#include <tensors/site/env/EnvBase.h>

/*! \brief Environment class with variance MPOs (i.e. double layer of energy MPOs) for environment blocks och type Left or Right corresponding to a single site.
 */

class EnvVar final : public EnvBase {
    public:
    using EnvBase::enlarge;
    using EnvBase::EnvBase;
    using EnvBase::set_edge_dims;
    explicit EnvVar(std::string side_, const MpsSite &mps, const MpoSite &mpo);
    [[nodiscard]] EnvVar enlarge(const MpsSite &mps, const MpoSite &mpo) const;
    void                 refresh(const EnvVar &env, const MpsSite &mps, const MpoSite &mpo);
    void                 set_edge_dims(const MpsSite &MPS, const MpoSite &MPO) final;
};
