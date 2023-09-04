#pragma once
#include "math/f128.h"
#include <fmt/core.h>

template<typename Char>
struct fmt::formatter<f128_t, Char> : fmt::formatter<f128_t::format_type, Char> {
    protected:
    typedef fmt::formatter<f128_t::format_type, Char> base;
    fmt::detail::dynamic_format_specs<Char>           specs_;

    public:
#if FMT_VERSION >= 100000

    // Parses format specifiers stopping either at the end of the range or at the
    // terminating '}'.
    template<typename ParseContext>
    FMT_CONSTEXPR auto parse(ParseContext &ctx) -> const Char * {
        auto type = detail::type_constant<double, Char>::value;
        auto end  = detail::parse_format_specs(ctx.begin(), ctx.end(), specs_, ctx, type);
        return end;
    }
#else
    // Use this parser to populate specs_, so that we know what width and precision to use in quadmath_snprintf later in format()
    template<typename ParseContext>
    FMT_CONSTEXPR auto parse(ParseContext &ctx) -> decltype(ctx.begin()) {
        auto begin = ctx.begin(), end = ctx.end();
        if(begin == end) return begin;
        using handler_type                       = detail::dynamic_specs_handler<format_parse_context>;
        auto                                type = detail::type_constant<f128_t::format_type, Char>::value;
        detail::specs_checker<handler_type> handler(handler_type(specs_, ctx), type);
        detail::parse_format_specs(ctx.begin(), ctx.end(), handler);
        return base::parse(ctx);
    }
#endif
    auto format(f128_t num, format_context &ctx) const -> decltype(ctx.out()) {
        // In quadmath, precision is the number of digits, and width is the total number of chars in the string.
        char prsnt = 'g';
        switch(specs_.type) {
                /* clang-format off */
            case fmt::presentation_type::exp_lower: {prsnt = 'e'; break  ;}
            case fmt::presentation_type::exp_upper: {prsnt = 'E'; break  ;}
            case fmt::presentation_type::fixed_lower: {prsnt = 'f'; break  ;}
            case fmt::presentation_type::fixed_upper: {prsnt = 'F'; break  ;}
            case fmt::presentation_type::general_lower: {prsnt = 'g'; break  ;}
            case fmt::presentation_type::general_upper: {prsnt = 'G'; break  ;}
            default: prsnt='g';
                /* clang-format on */
        }
        return fmt::format_to(ctx.out(), "{}", num.string(specs_.precision, specs_.width, prsnt, specs_.align == fmt::align_t ::left ? "-" : ""));
    }
};
