//
// Created by david on 2019-10-24.
//

#include <state/class_state_infinite.h>
#include <tools/nmspc_tools.h>
#include <model/class_model_base.h>
#include <model/class_model_factory.h>
#include <general/nmspc_random_numbers.h>
void tools::infinite::mpo::initialize(class_state_infinite & state, std::string model_type){
    tools::log->trace("Initializing mpo");
    //Generate MPO
    state.HA = class_model_factory::create_mpo(0,model_type);
    state.HB = class_model_factory::create_mpo(1,model_type);
}

void tools::infinite::mpo::randomize(class_state_infinite &state, int seed_model) {
    tools::log->trace("Setting random fields in MPO's");
    if (seed_model >= 0){
        rn::seed(seed_model);
    }

    state.HA->randomize_hamiltonian();
    state.HB->randomize_hamiltonian();

    std::vector<std::vector<double>> all_params;

    all_params.push_back(state.HA->get_parameter_values());
    all_params.push_back(state.HB->get_parameter_values());

    state.HA->set_full_lattice_parameters(all_params);
    state.HB->set_full_lattice_parameters(all_params);

}