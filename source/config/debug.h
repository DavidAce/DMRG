#pragma once
namespace settings {

#ifdef NDEBUG
#pragma message "DISABLE constexpr debug IN RELEASE MODE"
    inline constexpr bool debug = true;
#else
    inline constexpr bool debug = true;
#endif
    inline constexpr bool debug_split = false;
    inline constexpr bool debug_merge = false;
    inline constexpr bool debug_gates = false;
    inline constexpr bool debug_moves = false;
}