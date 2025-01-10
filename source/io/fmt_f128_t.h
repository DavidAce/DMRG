#pragma once
#include "math/float.h"
#include <fmt/core.h>

template<typename Char>
struct fmt::formatter<f128_t, Char> : fmt::formatter<f128_t::format_type, Char> {
    protected:
    typedef fmt::formatter<f128_t::format_type, Char> base;
    fmt::detail::dynamic_format_specs<Char>           specs_;

    public:
    // Parses format specifiers stopping either at the end of the range or at the
    // terminating '}'.
    template<typename ParseContext>
    FMT_CONSTEXPR auto parse(ParseContext &ctx) -> const Char * {
        // std::printf("parse\n");
        auto type = detail::type_constant<double, Char>::value;
        auto end  = detail::parse_format_specs(ctx.begin(), ctx.end(), specs_, ctx, type);
        return end;
    }
    auto format(f128_t num, format_context &ctx) const -> decltype(ctx.out()) {
        // In quadmath, precision is the number of digits, and width is the total number of chars in the string.
        // std::printf("format\n");
        char prsnt = 'g';
        switch(specs_.type) {
                /* clang-format off */
            case fmt::presentation_type::exp: {prsnt = 'e'; break  ;}
            case fmt::presentation_type::fixed: {prsnt = 'f'; break  ;}
            case fmt::presentation_type::general: {prsnt = 'g'; break  ;}
            default: prsnt='g';
                /* clang-format on */
        }
        return fmt::format_to(ctx.out(), "{}",
                              num.string(specs_.precision, specs_.width, prsnt,
                                         specs_.align == fmt::align_t ::left      ? "<"
                                         : ">", specs_.sign == fmt::sign_t::plus ? "+"
                                                                                  : ""));
    }
};
