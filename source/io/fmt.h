#pragma once

#include <fmt/compile.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/std.h>
#if defined(FMT_HEADER_ONLY)
    #pragma message "{fmt} has been included as header-only library. This causes large compile-time overhead"
#endif

#if defined(USE_QUADMATH)
    #include <quadmath.h>
template<typename Char>
struct fmt::formatter<__float128, Char> : fmt::formatter<double, Char> {
    auto format(__float128 c, format_context &ctx) const {
        std::string dec;
        dec.resize(64);
        quadmath_snprintf(dec.data(), 64, "%.36Qg", c);
        return fmt::format_to(ctx.out(), "{}", dec);
    }
};

#endif
#if defined(FMT_FORMAT_H_) && !defined(FMT_USE_COMPLEX)
    #define FMT_USE_COMPLEX 1
    #include <complex>
    #include <type_traits>

template<typename T, typename Char>
struct fmt::formatter<std::complex<T>, Char> : fmt::formatter<T, Char> {
    typedef fmt::formatter<T, Char> base;
    enum style { expr, pair } style_ = expr;
    fmt::detail::dynamic_format_specs<Char> specs_;
    FMT_CONSTEXPR auto                      parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
        using handler_type                            = fmt::detail::dynamic_specs_handler<format_parse_context>;
        auto                                     type = fmt::detail::type_constant<T, Char>::value;
        fmt::detail::specs_checker<handler_type> handler(handler_type(specs_, ctx), type);
        auto                                     it = ctx.begin();
        if(it != ctx.end()) {
            switch(*it) {
                case '$':
                    style_ = style::expr;
                    ctx.advance_to(++it);
                    break;
                case ',':
                    style_ = style::pair;
                    ctx.advance_to(++it);
                    break;
                default: break;
            }
        }
        parse_format_specs(ctx.begin(), ctx.end(), handler);
        // todo: fixup alignment
        return base::parse(ctx);
    }
    template<typename FormatCtx>
    auto format(const std::complex<T> &x, FormatCtx &ctx) const -> decltype(ctx.out()) {
        format_to(ctx.out(), "(");
        if(style_ == style::pair) {
            base::format(x.real(), ctx);
            format_to(ctx.out(), ",");
            base::format(x.imag(), ctx);
            return format_to(ctx.out(), ")");
        }
        if(x.real() != 0.0 || x.imag() == 0.0) base::format(x.real(), ctx);
        if(x.imag() != 0.0) {
            if(x.real() >= 0 && x.imag() >= 0 && specs_.sign != sign::plus) format_to(ctx.out(), "+");
            base::format(x.imag(), ctx);
            format_to(ctx.out(), "i");
            if(std::is_same<typename std::decay<T>::type, float>::value) format_to(ctx.out(), "f");
            if(std::is_same<typename std::decay<T>::type, long double>::value) format_to(ctx.out(), "l");
        }
        return format_to(ctx.out(), ")");
    }
};

#endif

#define FMT_EXTERN extern
#include "fmt.txx"
#undef FMT_EXTERN
