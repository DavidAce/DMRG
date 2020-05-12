#pragma once
#include <simulation/enums.h>
#include <string>

class class_state_infinite;

namespace tools::infinite::mpo {
    extern void initialize(class_state_infinite &state, ModelType model_type);
    extern void randomize(class_state_infinite &state);
}