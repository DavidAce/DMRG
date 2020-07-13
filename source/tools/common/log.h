#pragma once
#if __has_include(<spdlog/logger.h>)
#include <spdlog/logger.h>
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