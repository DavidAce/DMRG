#pragma once

#if !defined(SPDLOG_COMPILED_LIB)
    #if !defined(SPDLOG_HEADER_ONLY)
        #define SPDLOG_HEADER_ONLY
    #endif
#endif

#if __has_include(<spdlog/fmt/fmt.h>)
    // Spdlog will include the bundled fmt unless SPDLOG_FMT_EXTERNAL is defined, in which case <fmt/core.h> gets included instead
    // If SPDLOG_HEADER_ONLY is defined this will cause FMT_HEADER_ONLY to also get defined
    #include <spdlog/fmt/fmt.h>
    #if defined(SPDLOG_FMT_EXTERNAL)
        #include <fmt/compile.h>
        #include <fmt/ostream.h>
        #include <fmt/ranges.h>
    #else
        #include <spdlog/fmt/bundled/compile.h>
        #include <spdlog/fmt/bundled/ostream.h>
        #include <spdlog/fmt/bundled/ranges.h>
    #endif
#elif __has_include(<fmt/core.h>) &&  __has_include(<fmt/format.h>) && __has_include(<fmt/ranges.h>) &&  __has_include(<fmt/ostream.h>)
    #if defined(SPDLOG_HEADER_ONLY)
    // Since spdlog is header-only, let's assume fmt is as well
    // We do this because we have no way of knowing if this is getting linked to libfmt
        #define FMT_HEADER_ONLY
    #endif
    #include <fmt/compile.h>
    #include <fmt/core.h>
    #include <fmt/format.h>
    #include <fmt/ostream.h>
    #include <fmt/ranges.h>
#endif
