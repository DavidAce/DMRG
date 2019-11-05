//
// Created by david on 2019-06-25.
//
#include <state/class_state_finite.h>
#include <tools/nmspc_tools.h>
#include <model/class_model_base.h>
#include <model/class_model_factory.h>
#include <general/nmspc_random_numbers.h>

void tools::finite::mpo::initialize(class_state_finite & state, const size_t length, std::string model_type){
    tools::log->trace("Initializing mpo");
    //Generate MPO
    size_t pos = 0;
    state.MPO_L.emplace_back(class_model_factory::create_mpo(pos++,model_type));
    while(true){
        state.MPO_R.emplace_back(class_model_factory::create_mpo(pos++,model_type));
        if(state.MPO_L.size() + state.MPO_R.size() >= length){break;}
    }
}


void tools::finite::mpo::randomize(class_state_finite &state, int seed_model) {
    tools::log->trace("Setting random fields in MPO's");
    if (seed_model >= 0){
        rn::seed(seed_model);
    }
    std::vector<std::vector<double>> all_params;



    for (auto &mpo : state.MPO_L){
        mpo->randomize_hamiltonian();
        all_params.push_back(mpo->get_parameter_values());
    }
    for (auto &mpo : state.MPO_R){
        mpo->randomize_hamiltonian();
        all_params.push_back(mpo->get_parameter_values());
    }

    //

    for (auto &mpo : state.MPO_L){
        mpo->set_full_lattice_parameters(all_params,false);
    }
    for (auto &mpo : state.MPO_R){
        mpo->set_full_lattice_parameters(all_params,false);
    }
}



void tools::finite::mpo::reduce_mpo_energy(class_state_finite &state){
    state.unset_measurements();
    state.clear_cache();
    double energy_per_site   = tools::finite::measure::energy_per_site(state);
    tools::log->trace("Reducing MPO energy by: {}",energy_per_site);
    state.set_reduced_energy(energy_per_site);
}