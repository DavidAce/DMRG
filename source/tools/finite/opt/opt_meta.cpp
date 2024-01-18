#include "../opt_meta.h"
#include "debug/exceptions.h"
#include "config/enums.h"

namespace tools::finite::opt {
    OptMeta::OptMeta()
        : optFunc(OptFunc::VARIANCE), optSolver(OptSolver::EIGS), optType(OptType::CPLX), optInit(OptInit::CURRENT_STATE), optWhen(OptWhen::ALWAYS),
          optRitz(OptRitz::SR), optExit(OptExit::NONE) {}

    OptMeta::OptMeta(OptRitz ritz, OptFunc mode) : OptMeta() {
        optRitz = ritz;
        optFunc = mode;
    }

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
//        if(optFunc == OptFunc::OVERLAP and optSolver == OptSolver::BFGS) throw except::runtime_error("opt: mode [OVERLAP] and solver [BFGS] are incompatible");
//        if(optFunc == OptFunc::SUBSPACE and optSolver == OptSolver::BFGS) throw except::runtime_error("opt: mode [ENERGY] and solver [BFGS] are incompatible");
//        if(optFunc == OptFunc::ENERGY and optSolver == OptSolver::BFGS) throw except::runtime_error("opt: mode [ENERGY] and solver [BFGS] are incompatible");
    }

}
