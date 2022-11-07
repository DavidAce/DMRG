#pragma once
#if defined(FMT_HEADER_ONLY)
    #pragma message "{fmt} has been included as header-only library. This causes large compile-time overhead"
#endif

#if !defined(SPDLOG_FMT_EXTERNAL)
    #define SPDLOG_FMT_EXTERNAL
#endif
#if !defined(SPDLOG_COMPILED_LIB)
    #if !defined(SPDLOG_HEADER_ONLY)
        #define SPDLOG_HEADER_ONLY
    #endif
    #pragma message "SPDLOG_COMPILED_LIB is not defined. Define it and link libspdlog to speed up compilation"
#endif

#include <fmt/compile.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

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

namespace eig {
    inline auto log = spdlog::get("eig") == nullptr ? spdlog::stdout_color_mt("eig", spdlog::color_mode::always) : spdlog::get("eig");
    extern void setLevel(spdlog::level::level_enum level);
    extern void setLevel(size_t level);
    extern void setTimeStamp(std::string_view stamp = "[%Y-%m-%d %H:%M:%S][%n]%^[%=8l]%$ %v");
}