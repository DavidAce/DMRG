//
// Created by david on 2019-03-27.
//

#pragma once
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#if defined(SPDLOG_FMT_EXTERNAL)
#include <fmt/ranges.h>
#include <fmt/ostream.h>
#else
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/fmt/bundled/ostream.h>
#endif

namespace Logger{
    extern void enableTimestamp(std::shared_ptr<spdlog::logger> &log);
    extern void disableTimestamp(std::shared_ptr<spdlog::logger> &log);
    extern size_t getLogLevel(std::shared_ptr<spdlog::logger> &log);
    template<typename levelType>
    extern void setLogLevel(std::shared_ptr<spdlog::logger> &log, levelType levelZeroToFive);
    extern void setLogger(std::shared_ptr<spdlog::logger> &log,const std::string &name, std::optional<size_t> levelZeroToFive = std::nullopt, std::optional<bool> timestamp = std::nullopt);
    extern std::shared_ptr<spdlog::logger>  setLogger(const std::string &name, std::optional<size_t> levelZeroToFive = std::nullopt, std::optional<bool> timestamp = std::nullopt);
    extern std::shared_ptr<spdlog::logger>  getLogger(const std::string & name );
}

