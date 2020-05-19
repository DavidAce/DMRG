//
// Created by david on 2019-02-20.
//

#include "print.h"
#include <tensors/state/class_state_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/model/class_mpo_site.h>
#include <tools/common/log.h>

void tools::infinite::print::print_hamiltonians(const class_model_infinite &model) {
    model.get_mpo_siteA().print_parameter_names();
    model.get_mpo_siteA().print_parameter_values();
    model.get_mpo_siteB().print_parameter_values();
}

void tools::infinite::print::print_state(const class_state_infinite &state) { tools::log->warn("Print state not implemented yet"); }