#pragma once
#include "fmt.h"

#if __has_include(<spdlog/spdlog.h>) && __has_include(<spdlog/sinks/stdout_color_sinks.h>)
    #include <spdlog/sinks/stdout_color_sinks.h>
    #include <spdlog/spdlog.h>

    #if !defined(SPDLOG_COMPILED_LIB)
        #pragma message "SPDLOG_COMPILED_LIB is not defined. Define it and link libspdlog to get shorter compilation"
    #endif

    #if defined(FMT_HEADER_ONLY)
        #pragma message "{fmt} has been included as header-only library. This causes large compile-time overhead"
    #endif
#endif
