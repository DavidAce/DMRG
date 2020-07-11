#pragma once

#if __has_include(<spdlog/logger.h>)
#include <spdlog/logger.h>
#else
    #error "Could not find spdlog headers"
#endif

#if !defined(SPDLOG_FMT_EXTERNAL) && __has_include(<spdlog/fmt/fmt.h>) && __has_include(<spdlog/fmt/bundled/ranges.h>) &&  __has_include(<spdlog/fmt/bundled/ostream.h>)
#include <spdlog/fmt/bundled/ostream.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/fmt/fmt.h>
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


namespace Logger{
    extern void enableTimestamp (const std::shared_ptr<spdlog::logger> &log);
    extern void disableTimestamp(const std::shared_ptr<spdlog::logger> &log);
    extern size_t getLogLevel   (const std::shared_ptr<spdlog::logger> &log);
    template<typename levelType>
    extern void setLogLevel(const std::shared_ptr<spdlog::logger> &log, levelType levelZeroToFive);
    extern void setLogger  (std::shared_ptr<spdlog::logger> &log, const std::string &name, std::optional<size_t> levelZeroToFive = std::nullopt, std::optional<bool> timestamp = std::nullopt);
    extern std::shared_ptr<spdlog::logger>  setLogger(const std::string &name, std::optional<size_t> levelZeroToFive = std::nullopt, std::optional<bool> timestamp = std::nullopt);
    extern std::shared_ptr<spdlog::logger>  getLogger(const std::string & name );
}

