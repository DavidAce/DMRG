#pragma once
#include <string>

class class_state_infinite;
namespace tools::infinite::mps {
    extern void random_product_state(const class_state_infinite &state, [[maybe_unused]] const std::string &sector, [[maybe_unused]] long bitfield, bool use_eigenspinors );
}