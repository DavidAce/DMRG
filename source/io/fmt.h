#pragma once

#if defined(SPDLOG_FMT_EXTERNAL) && \
    __has_include(<fmt/core.h>) &&  __has_include(<fmt/format.h>) && __has_include(<fmt/ranges.h>) &&  __has_include(<fmt/ostream.h>)
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#elif !defined(SPDLOG_FMT_EXTERNAL) && \
    __has_include(<spdlog/fmt/fmt.h>) && __has_include(<spdlog/fmt/bundled/ranges.h>) &&  __has_include(<spdlog/fmt/bundled/ostream.h>)
#include <spdlog/fmt/bundled/ostream.h>
    #include <spdlog/fmt/bundled/ranges.h>
    #include <spdlog/fmt/fmt.h>
    #if defined(FMT_HEADER_ONLY)
        #pragma message "{fmt} has been included as header-only library. This causes large compile-time overhead"
    #endif
#elif __has_include(<fmt/core.h>) &&  __has_include(<fmt/format.h>) && __has_include(<fmt/ranges.h>) &&  __has_include(<fmt/ostream.h>)
// Check if there are already fmt headers installed independently from Spdlog
// Note that in this case the user hasn't enabled Spdlog for h5pp, so the build hasn't linked any compiled FMT libraries
// To avoid undefined references we coult opt in to the header-only mode of FMT.
// Note that this check should be skipped if using conan. Then, SPDLOG_FMT_EXTERNAL is defined
    #include <fmt/core.h>
    #include <fmt/format.h>
    #include <fmt/ostream.h>
    #include <fmt/ranges.h>
#else
    #error "Could not find fmt headers"
#endif
