//
// Created by david on 2019-10-24.
//

#include <tools/nmspc_tools.h>
#include <state/class_state_infinite.h>
#include <state/class_environment.h>
#include <state/class_mps_2site.h>
#include <state/class_mps_site.h>
#include <model/class_model_base.h>
void tools::infinite::env::initialize  (class_state_infinite & state){
    state.Lblock = std::make_unique<class_environment>("L",0);
    state.Rblock = std::make_unique<class_environment>("R",1);
    state.Lblock2 = std::make_unique<class_environment_var>("L",0);
    state.Rblock2 = std::make_unique<class_environment_var>("R",1);
}
