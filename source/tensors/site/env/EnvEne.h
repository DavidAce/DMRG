#pragma once
#include "tensors/site/env/EnvBase.h"

/*! \brief Environment class with energy MPOs for environment blocks och type Left or Right corresponding to a single site.
 */

class EnvEne final : public EnvBase {
    public:
    using EnvBase::enlarge;
    using EnvBase::EnvBase;
    using EnvBase::set_edge_dims;
    explicit EnvEne(std::string side_, const MpsSite &mps, const MpoSite &mpo);
    [[nodiscard]] EnvEne enlarge(const MpsSite &mps, const MpoSite &mpo) const;
    void                 refresh(const EnvEne &env, const MpsSite &mps, const MpoSite &mpo);
    void                 set_edge_dims(const MpsSite &mps, const MpoSite &mpo) final;
};
