//
// Created by david on 2019-06-25.
//

#include <tools/nmspc_tools.h>
#include <state/class_infinite_state.h>


class_infinite_state tools::infinite::mps::set_random_state(const class_infinite_state & state, std::string parity, int seed_state){
    throw std::runtime_error("You need to implement set random state for infinite state");
    return state;
}
