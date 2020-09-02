#include <tools/common/log.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

void tools::Logger::enableTimestamp(const std::shared_ptr <spdlog::logger> &other_log) {
    if (other_log != nullptr) {
        other_log->set_pattern("[%Y-%m-%d %H:%M:%S][%n]%^[%=8l]%$ %v");
        other_log->trace("Enabled timestamp");
    }
}

void tools::Logger::disableTimestamp(const std::shared_ptr<spdlog::logger> &other_log){
    if(other_log != nullptr){
        other_log->set_pattern("[%n]%^[%=8l]%$ %v");
        other_log->trace("Disabled timestamp");
    }
}


template<typename levelType>
void tools::Logger::setLogLevel(const std::shared_ptr<spdlog::logger> &other_log, levelType levelZeroToFive) {
    if constexpr(std::is_same_v<levelType, spdlog::level::level_enum>)
        other_log->set_level(levelZeroToFive);
    else if constexpr(std::is_integral_v<levelType>) {
        if(levelZeroToFive > 5) { throw std::runtime_error("Expected verbosity level integer in [0-5]. Got: " + std::to_string(levelZeroToFive)); }
        return tools::Logger::setLogLevel(other_log,static_cast<spdlog::level::level_enum>(levelZeroToFive));
    }else if constexpr(std::is_same_v<levelType, std::optional<size_t>>) {
        if(levelZeroToFive)
            return setLogLevel(other_log,levelZeroToFive.value());
        else
            return;
    }else if constexpr(std::is_same_v<levelType, std::optional<spdlog::level::level_enum>>){
        if(levelZeroToFive) return setLogLevel(other_log,levelZeroToFive.value());
        else return;
    } else {
        throw std::runtime_error("Given wrong type for spdlog verbosity level");
    }
//    other_log->info("Log verbosity level: {}   | trace:0 | debug:1 | info:2 | warn:3 | error:4 | critical:5 |", static_cast<int>(other_log->level()));
//    other_log->debug("Log verbosity level: {}", static_cast<int>(other_log->level()));
}

template void tools::Logger::setLogLevel(const std::shared_ptr<spdlog::logger> &other_log, size_t levelZeroToFive);
template void tools::Logger::setLogLevel(const std::shared_ptr<spdlog::logger> &other_log, std::optional<size_t> levelZeroToFive);
template void tools::Logger::setLogLevel(const std::shared_ptr<spdlog::logger> &other_log, std::optional<spdlog::level::level_enum> levelZeroToFive);



size_t tools::Logger::getLogLevel(const std::shared_ptr<spdlog::logger> &other_log) {
    if(other_log != nullptr)
        return static_cast<size_t>(other_log->level());
    else
        return 2;
}


void tools::Logger::setLogger(std::shared_ptr<spdlog::logger> &other_log, const std::string &name, std::optional<size_t> levelZeroToFive , std::optional<bool> timestamp) {
    if(spdlog::get(name) == nullptr)
        other_log = spdlog::stdout_color_st(name,spdlog::color_mode::always);
    else
        other_log = spdlog::get(name);
    other_log->set_pattern("[%n]%^[%=8l]%$ %v"); // Disabled timestamp is the default
    setLogLevel(other_log,levelZeroToFive);
    if(timestamp and timestamp.value())
        enableTimestamp(other_log);
}

std::shared_ptr<spdlog::logger> tools::Logger::setLogger(const std::string &name, std::optional<size_t> levelZeroToFive, std::optional<bool> timestamp) {
    std::shared_ptr<spdlog::logger> other_log;
    if(spdlog::get(name) == nullptr) {
        other_log = spdlog::stdout_color_mt(name);
    } else {
        other_log = spdlog::get(name);
    }
    tools::Logger::setLogLevel(other_log,levelZeroToFive);
    if(not timestamp or (timestamp and timestamp.value()))
        tools::Logger::enableTimestamp(other_log);
    else if(timestamp and not timestamp.value())
        tools::Logger::disableTimestamp(other_log);
    return other_log;
}


std::shared_ptr<spdlog::logger> tools::Logger::getLogger( const std::string & name ){
    if(spdlog::get(name) == nullptr)
        throw std::runtime_error(fmt::format("Logger with name [{}] does not exist",name));
    else
        return spdlog::get(name);
}
