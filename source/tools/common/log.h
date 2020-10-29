#pragma once
#if __has_include(<spdlog/logger.h>)
    #include <spdlog/logger.h>
    #if !defined(SPDLOG_COMPILED_LIB)
        #pragma message "SPDLOG_COMPILED_LIB is not defined. Define it and link libspdlog to get shorter compilation"
    #endif
    #if defined(FMT_HEADER_ONLY)
    #pragma message "{fmt} has been included as header-only library. This causes large compile-time overhead"
    #endif
#else
#error "Could not find spdlog headers"
#endif

namespace tools {
    inline std::shared_ptr<spdlog::logger> log;
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

}