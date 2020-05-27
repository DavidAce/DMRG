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
    extern void enableTimeStamp(std::shared_ptr<spdlog::logger> &log);
    extern void disableTimeStamp(std::shared_ptr<spdlog::logger> &log);
    extern void setLogLevel(std::shared_ptr<spdlog::logger> &log, size_t levelZeroToSix);
    extern size_t getLogLevel(std::shared_ptr<spdlog::logger> &log);
    extern void setLogger(std::shared_ptr<spdlog::logger> &log, const std::string & name, size_t levelZeroToSix = 2, bool timestamp = true);
    extern std::shared_ptr<spdlog::logger>  setLogger( const std::string & name, size_t levelZeroToSix = 2, bool timestamp = true);
    extern std::shared_ptr<spdlog::logger>  getLogger( const std::string & name );
}

