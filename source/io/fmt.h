#pragma once

#include <fmt/compile.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/std.h>
#include "fmt_f128_t.h"
#if defined(FMT_HEADER_ONLY)
    #pragma message "{fmt} has been included as header-only library. This causes large compile-time overhead"
#endif

#if defined(FMT_FORMAT_H_) && !defined(FMT_USE_COMPLEX)
    #define FMT_USE_COMPLEX 1
    #include <complex>
    #include <type_traits>

template<typename T, typename Char>
struct fmt::formatter<std::complex<T>, Char> : fmt::formatter<T, Char> {
    private:
    typedef fmt::formatter<T, Char>         base;
    fmt::detail::dynamic_format_specs<Char> specs_;

    public:
    template<typename FormatCtx>
    auto format(const std::complex<T> &x, FormatCtx &ctx) const -> decltype(ctx.out()) {
        base::format(x.real(), ctx);
        if(x.imag() >= 0 && specs_.sign != sign::plus) format_to(ctx.out(), "+");
        base::format(x.imag(), ctx);
        return format_to(ctx.out(), "i");
    }
};

#endif
#if defined(FMT_FORMAT_H_)
    #include <complex>
    #include <type_traits>

template<typename T>
struct fmt::formatter<std::reference_wrapper<T>> : fmt::formatter<T> {
    private:
    typedef fmt::formatter<T> base;

    public:
    template<typename FormatCtx>
    auto format(const std::reference_wrapper<T> &x, FormatCtx &ctx) const -> decltype(ctx.out()) {
        base::format(static_cast<T>(x), ctx);
        return ctx.out();
    }
};
#endif

#define FMT_EXTERN extern
#include "fmt.txx"
#undef FMT_EXTERN
