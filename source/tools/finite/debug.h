#pragma once

class class_state_finite;
namespace tools::finite::debug {
    extern void check_integrity             (const class_state_finite & state);
    extern void check_integrity_of_mps      (const class_state_finite & state);
    extern void check_integrity_of_mpo      (const class_state_finite & state);
    extern void check_normalization_routine (const class_state_finite & state);
    extern void print_parity_properties     (const class_state_finite & state);
}