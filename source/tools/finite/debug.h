#pragma once

class class_state_finite;
class class_model_finite;
class class_edges_finite;

namespace tools::finite::debug {
    /* clang-format off */
    extern void check_integrity             (const class_state_finite & state,
                                             const class_model_finite & model,
                                             const class_edges_finite & edges);
    extern void check_integrity_of_state    (const class_state_finite & state);
    extern void check_integrity_of_model    (const class_model_finite & model);
    extern void check_integrity_of_edges    (const class_edges_finite & edges);
    extern void check_normalization_routine (const class_state_finite & state);
    extern void print_parity_properties     (const class_state_finite & state);
    /* clang-format on */
}