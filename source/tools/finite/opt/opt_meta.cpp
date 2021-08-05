#include "../opt_meta.h"
#include <config/enums.h>

namespace tools::finite::opt {
    OptMeta::OptMeta()
        : optMode(OptMode::VARIANCE), optSpace(OptSpace::DIRECT), optType(OptType::CPLX), optInit(OptInit::CURRENT_STATE), optWhen(OptWhen::ALWAYS),
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
                    return have_common(optWhen,previous_exit);
            }
        }
    }

}
