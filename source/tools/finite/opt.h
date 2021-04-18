//
// Created by david on 2019-03-18.
//

#pragma once

#include <general/eigen_tensor_fwd_decl.h>

class class_state_finite;
class class_model_finite;
class class_edges_finite;
class class_tensors_finite;
class class_algorithm_status;
class class_tic_toc;
namespace eig {class solver;}
enum class OptSpace;
enum class OptType;
enum class OptMode;
enum class StateRitz;
namespace tools::finite::opt {
    class opt_mps;

    using Scalar = std::complex<double>;
    using real = double;
    using cplx = std::complex<double>;

    extern opt_mps find_excited_state(const class_tensors_finite &tensors, const opt_mps &initial_mps, const class_algorithm_status &status,
                                        OptMode optMode, OptSpace optSpace, OptType optType);
    extern opt_mps find_excited_state(const class_tensors_finite &tensors, const class_algorithm_status &status, OptMode optMode, OptSpace optSpace,
                                        OptType optType);
    extern opt_mps find_ground_state(const class_tensors_finite &tensors, const class_algorithm_status &status, StateRitz ritz);
}
