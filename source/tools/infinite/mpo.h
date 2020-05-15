#pragma once
#include <config/enums.h>
#include <string>

class class_model_infinite;

namespace tools::infinite::mpo {
    extern void initialize(class_model_infinite &state, ModelType model_type);
    extern void randomize(class_model_infinite &model);
}