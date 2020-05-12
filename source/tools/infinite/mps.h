#pragma once
#include <string>

class class_state_infinite;
namespace tools::infinite::mps {
    extern void                 initialize(class_state_infinite &state, ModelType model_type);
    extern class_state_infinite set_random_state(const class_state_infinite &state, [[maybe_unused]] const std::string &parity);
}