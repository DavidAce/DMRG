//
// Created by david on 2019-02-20.
//

#include <iostream>
#include <tensors/model/class_mpo_base.h>
#include <tensors/state/class_mps_2site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/infinite/print.h>

void tools::infinite::print::print_hamiltonians(const class_state_infinite & state){
    state.HA->print_parameter_names();
    state.HA->print_parameter_values();
    state.HB->print_parameter_values();

}

void tools::infinite::print::print_state(const class_state_infinite & state){
    using namespace Textra;
    std::cout << std::setprecision(10);
    std::cout << "State length              : "    << state.get_length()   << std::endl;
    std::cout << "State position            : "    << state.get_position() << std::endl;
}