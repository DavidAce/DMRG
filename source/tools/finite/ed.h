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

namespace tools::finite::ed {
    extern StateFinite find_exact_state(const TensorsFinite &tensors, const AlgorithmStatus &status);
}