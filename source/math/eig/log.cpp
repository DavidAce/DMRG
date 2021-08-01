#include "log.h"

void eig::setLevel(spdlog::level::level_enum level) { eig::log->set_level(level); }
void eig::setLevel(size_t level) { eig::log->set_level(static_cast<spdlog::level::level_enum>(level)); }
void eig::setTimeStamp(std::string_view stamp) { eig::log->set_pattern(std::string(stamp)); }
