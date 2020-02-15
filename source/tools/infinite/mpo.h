#pragma once
#include <string>

class class_state_infinite;

namespace tools::infinite::mpo {
    extern void initialize(class_state_infinite &state, const std::string &model_type_str);
    extern void randomize(class_state_infinite &state);
}