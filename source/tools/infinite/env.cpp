//
// Created by david on 2019-10-24.
//

#include <tensors/state/class_environment.h>
#include <tensors/state/class_mps_2site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/infinite/env.h>
//#include <tensors/state/class_mps_site.h>
//#include <tensors/model/class_model_base.h>

void tools::infinite::env::initialize(class_state_infinite &state) {
    state.Lblock  = std::make_unique<class_environment>("L", *state.MPS->MPS_A, *state.HA);
    state.Rblock  = std::make_unique<class_environment>("R", *state.MPS->MPS_B, *state.HB);
    state.Lblock2 = std::make_unique<class_environment_var>("L", *state.MPS->MPS_A, *state.HA);
    state.Rblock2 = std::make_unique<class_environment_var>("R", *state.MPS->MPS_B, *state.HB);
}
