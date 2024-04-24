#pragma once
#include <string>
namespace tools {
    std::string fmtBytes(bool on, size_t bytes, size_t base = 1024, size_t decimals = 2);
}