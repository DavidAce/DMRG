//
// Created by david on 2019-10-24.
//

#include <tools/infinite/mpo.h>
#include <tools/common/log.h>
#include <state/class_state_infinite.h>
#include <model/class_model_base.h>
#include <model/class_model_factory.h>

void tools::infinite::mpo::initialize(class_state_infinite & state, const std::string  & model_type){
    tools::log->trace("Initializing mpo");
    //Generate MPO
    state.HA = class_model_factory::create_mpo(0,model_type);
    state.HB = class_model_factory::create_mpo(1,model_type);
}

void tools::infinite::mpo::randomize(class_state_infinite &state) {
    tools::log->trace("Setting random fields in MPO's");
    state.HA->randomize_hamiltonian();
    state.HB->randomize_hamiltonian();
    std::vector<class_model_base::TableMap> all_params;

    all_params.push_back(state.HA->get_parameters());
    all_params.push_back(state.HB->get_parameters());

    state.HA->set_averages(all_params);
    state.HB->set_averages(all_params);

}