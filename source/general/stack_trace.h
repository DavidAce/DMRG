#pragma once
#include <string>

namespace debug{
    void print_stack_trace();
    void throw_stack_trace(const std::string & msg);
}