#pragma once

#include <fmt/compile.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>

#if defined(FMT_HEADER_ONLY)
    #pragma message "{fmt} has been included as header-only library. This causes large compile-time overhead"
#endif

#if defined(FMT_FORMAT_H_) && !defined(FMT_USE_COMPLEX)
    #define FMT_USE_COMPLEX
    #include <complex>
    #include <type_traits>
template<typename T>
struct fmt::formatter<std::complex<T>, char, std::enable_if_t<std::is_arithmetic_v<T>>> : fmt::formatter<typename std::complex<T>::value_type> {
    template<typename FormatContext>
    auto format(const std::complex<T> &number, FormatContext &ctx) const {
        return fmt::format_to(ctx.out(), "{0}{1:+}i", number.real(), number.imag());
    }
};
#endif

#define FMT_EXTERN extern
#include "fmt.txx"
#undef FMT_EXTERN
