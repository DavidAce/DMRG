#pragma once

#include <general/eigen_tensor_fwd_decl.h>

class StateFinite;
class ModelFinite;
class EdgesFinite;
class TensorsFinite;
class AlgorithmStatus;
class ur;
namespace eig {
    class solver;
}
enum class OptSpace;
enum class OptType;
enum class OptMode;
enum class StateRitz;
namespace tools::finite::opt {
    class opt_mps;

    using real = double;
    using cplx = std::complex<double>;

    extern opt_mps find_excited_state(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status, OptMode optMode,
                                      OptSpace optSpace, OptType optType);
    extern opt_mps find_excited_state(const TensorsFinite &tensors, const AlgorithmStatus &status, OptMode optMode, OptSpace optSpace, OptType optType);
    extern opt_mps find_ground_state(const TensorsFinite &tensors, const AlgorithmStatus &status, StateRitz ritz);
}
