#include "../opt_meta.h"
#include "debug/exceptions.h"
#include <config/enums.h>

namespace tools::finite::opt {
    OptMeta::OptMeta()
        : optMode(OptMode::VARIANCE), optSolver(OptSolver::LBFGS), optType(OptType::CPLX), optInit(OptInit::CURRENT_STATE), optWhen(OptWhen::ALWAYS),
          optRitz(OptRitz::SR), optExit(OptExit::NONE) {}

    OptMeta::OptMeta(OptRitz ritz) : OptMeta() { optRitz = ritz; }

    bool OptMeta::should_proceed(OptExit previous_exit) const {
        switch(optWhen) {
            case OptWhen::ALWAYS: return true;
            case OptWhen::NEVER: return false;
            default: {
                if(previous_exit == OptExit::SUCCESS)
                    return false;
                else
                    return have_common(optWhen, previous_exit);
            }
        }
    }
    void OptMeta::validate() const {
        if(optMode == OptMode::OVERLAP and optSolver == OptSolver::LBFGS)
            throw except::runtime_error("opt: mode [OVERLAP] and solver [LBFGS] are incompatible");
        if(optMode == OptMode::SUBSPACE and optSolver == OptSolver::LBFGS)
            throw except::runtime_error("opt: mode [ENERGY] and solver [LBFGS] are incompatible");
        if(optMode == OptMode::ENERGY and optSolver == OptSolver::LBFGS) throw except::runtime_error("opt: mode [ENERGY] and solver [LBFGS] are incompatible");
    }

}
