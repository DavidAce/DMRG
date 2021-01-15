#pragma once
namespace settings {

#ifdef NDEBUG
    inline constexpr bool debug = false;
#else
    inline constexpr bool debug = true;
#endif
    inline constexpr bool debug_split = false;
    inline constexpr bool debug_merge = false;
    inline constexpr bool debug_gates = false;
    inline constexpr bool debug_moves = false;
}