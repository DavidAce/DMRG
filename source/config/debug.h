#pragma once
namespace settings {

#ifdef NDEBUG
    inline constexpr bool debug = false;
#else
    inline constexpr bool debug = true;
#endif
    // These are for very specific, very expensive debugging when developing those steps
    //    inline constexpr bool debug_split = false;
    //    inline constexpr bool debug_merge = false; // Deprecated. Change this in tools/finite/mps/mps.cpp
    //    inline constexpr bool debug_gates = false; // Deprecated. Change this in tools/finite/mps/mps.cpp
    //    inline constexpr bool debug_moves = false; // Deprecated. Change this in tools/finite/mps/mps.cpp
    //    inline constexpr bool debug_numen = false; // Deprecated. Change this in tools/finite/measure/entropy.cpp
}