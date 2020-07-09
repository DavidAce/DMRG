//
// Created by david on 2019-10-16.
//

#include "nmspc_logger.h"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/stdout_sinks.h>

void Logger::enableTimestamp(std::shared_ptr <spdlog::logger> &log) {
    if (log != nullptr) {
        log->set_pattern("[%Y-%m-%d %H:%M:%S][%n]%^[%=8l]%$ %v");
        log->trace("Enabled timestamp");
    }
}

void Logger::disableTimestamp(std::shared_ptr<spdlog::logger> &log){
    if(log != nullptr){
        log->set_pattern("[%n]%^[%=8l]%$ %v");
        log->trace("Disabled timestamp");
    }
}


template<typename levelType>
void Logger::setLogLevel(std::shared_ptr<spdlog::logger> &log, levelType levelZeroToFive) {
    if constexpr(std::is_same_v<levelType, spdlog::level::level_enum>)
        log->set_level(levelZeroToFive);
    else if constexpr(std::is_integral_v<levelType>) {
        if(levelZeroToFive > 5) { throw std::runtime_error("Expected verbosity level integer in [0-5]. Got: " + std::to_string(levelZeroToFive)); }
        return setLogLevel(log,static_cast<spdlog::level::level_enum>(levelZeroToFive));
    }else if constexpr(std::is_same_v<levelType, std::optional<size_t>>) {
        if(levelZeroToFive)
            return setLogLevel(log,levelZeroToFive.value());
        else
            return;
    }else if constexpr(std::is_same_v<levelType, std::optional<spdlog::level::level_enum>>){
        if(levelZeroToFive) return setLogLevel(log,levelZeroToFive.value());
        else return;
    } else {
        throw std::runtime_error("Given wrong type for spdlog verbosity level");
    }
//        log->info("Log verbosity level: {}   | trace:0 | debug:1 | info:2 | warn:3 | error:4 | critical:5 |", static_cast<int>(log->level()));
    log->debug("Log verbosity level: {}", static_cast<int>(log->level()));
}

template void Logger::setLogLevel(std::shared_ptr<spdlog::logger> &log, size_t levelZeroToFive);
template void Logger::setLogLevel(std::shared_ptr<spdlog::logger> &log, std::optional<size_t> levelZeroToFive);
template void Logger::setLogLevel(std::shared_ptr<spdlog::logger> &log, std::optional<spdlog::level::level_enum> levelZeroToFive);



size_t Logger::getLogLevel(std::shared_ptr<spdlog::logger> &log) {
    if(log != nullptr)
        return static_cast<size_t>(log->level());
    else
        return 2;
}


void Logger::setLogger(std::shared_ptr<spdlog::logger> &log, const std::string &name, std::optional<size_t> levelZeroToFive , std::optional<bool> timestamp) {
    if(spdlog::get(name) == nullptr)
        log = spdlog::stdout_color_mt(name);
    else
        log = spdlog::get(name);
    log->set_pattern("[%n]%^[%=8l]%$ %v"); // Disabled timestamp is the default
    setLogLevel(log,levelZeroToFive);
    if(timestamp and timestamp.value())
        enableTimestamp(log);
}

std::shared_ptr<spdlog::logger> Logger::setLogger(const std::string &name, std::optional<size_t> levelZeroToFive, std::optional<bool> timestamp) {
    std::shared_ptr<spdlog::logger> log;
    if(spdlog::get(name) == nullptr) {
        log = spdlog::stdout_color_mt(name);
    } else
        log = spdlog::get(name);
    log->set_pattern("[%n]%^[%=8l]%$ %v"); // Disabled timestamp is the default
    Logger::setLogLevel(log,levelZeroToFive);
    if(timestamp and timestamp.value())
        Logger::enableTimestamp(log);
    return log;
}


std::shared_ptr<spdlog::logger> Logger::getLogger( const std::string & name ){
    if(spdlog::get(name) == nullptr)
        throw std::runtime_error(fmt::format("Logger with name [{}] does not exist",name));
    else
        return spdlog::get(name);
}
