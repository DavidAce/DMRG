//
// Created by david on 2019-10-24.
//

#include <state/class_state_infinite.h>
#include <tools/nmspc_tools.h>
#include <model/class_model_base.h>
#include <model/class_model_factory.h>
#include <math/nmspc_random.h>
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
    std::vector<class_model_base::Parameters> all_params;

    all_params.push_back(state.HA->get_parameters());
    all_params.push_back(state.HB->get_parameters());

    state.HA->set_full_lattice_parameters(all_params);
    state.HB->set_full_lattice_parameters(all_params);

}