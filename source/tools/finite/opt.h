#pragma once

#include "math/tenx/fwd_decl.h"

class StateFinite;
class ModelFinite;
class EdgesFinite;
class TensorsFinite;
class AlgorithmStatus;
class ur;
namespace eig {
    class solver;
}

namespace tools::finite::opt {
    class opt_mps;
    struct OptMeta;
    extern opt_mps get_opt_initial_mps(const TensorsFinite &tensors, const OptMeta &meta);
    extern opt_mps get_updated_state(const TensorsFinite &tensors, const opt_mps &initial_state, const AlgorithmStatus &status, OptMeta &meta);
    extern opt_mps find_ground_state(const TensorsFinite &tensors, const opt_mps &initial_state, const AlgorithmStatus &status, OptMeta &meta);
}
