//
// Created by david on 2019-10-16.
//

#include "logger.h"
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/stdout_sinks.h>
#include <spdlog/spdlog.h>

void tools::logger::enableTimeStamp(std::shared_ptr<spdlog::logger> &log) {
    if(log != nullptr) { log->set_pattern("[%Y-%m-%d %H:%M:%S][%n]%^[%=8l]%$ %v"); }
}

void tools::logger::disableTimeStamp(std::shared_ptr<spdlog::logger> &log) {
    if(log != nullptr) { log->set_pattern("[%n]%^[%=8l]%$ %v"); }
}

void tools::logger::setLogLevel(std::shared_ptr<spdlog::logger> &log, size_t levelZeroToSix) {
    if(levelZeroToSix > 6) { throw std::logic_error("Expected verbosity level integer in [0-6]. Got: " + std::to_string(levelZeroToSix)); }
    auto lvlEnum = static_cast<spdlog::level::level_enum>(levelZeroToSix);

    // Set console settings
    log->set_level(lvlEnum);
}

size_t tools::logger::getLogLevel(std::shared_ptr<spdlog::logger> &log) { return static_cast<size_t>(log->level()); }

void tools::logger::setLogger(std::shared_ptr<spdlog::logger> &log, const std::string &name, size_t levelZeroToSix, bool timestamp) {
    if(spdlog::get(name) == nullptr) {
        log = spdlog::stdout_color_mt(name, spdlog::color_mode::always);
        //            log = spdlog::stdout_color_mt(name);

        if(timestamp) {
            enableTimeStamp(log);
        } else {
            disableTimeStamp(log);
        }
        setLogLevel(log, levelZeroToSix);
    } else {
        log = spdlog::get(name);
    }
}

std::shared_ptr<spdlog::logger> tools::logger::setLogger(const std::string &name, size_t levelZeroToSix, bool timestamp) {
    if(spdlog::get(name) == nullptr) {
        auto log_ = spdlog::stdout_color_mt(name, spdlog::color_mode::always);

        //            auto log = spdlog::stdout_color_mt(name);

        if(timestamp) {
            enableTimeStamp(log_);
        } else {
            disableTimeStamp(log_);
        }
        setLogLevel(log_, levelZeroToSix);
        return log_;
    } else {
        return spdlog::get(name);
    }
}

std::shared_ptr<spdlog::logger> tools::logger::getLogger(const std::string &name) {
    if(spdlog::get(name) == nullptr)
        throw std::runtime_error(fmt::format("Logger with name [{}] does not exist", name));
    else
        return spdlog::get(name);
}
