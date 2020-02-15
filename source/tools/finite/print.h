#pragma once

class class_state_finite;
namespace tools::finite::print {
    extern void print_full_state    (const class_state_finite & state);
    extern void print_state         (const class_state_finite & state); /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
    extern void print_state_compact (const class_state_finite & state); /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
    extern void print_hamiltonians  (const class_state_finite & state);
}