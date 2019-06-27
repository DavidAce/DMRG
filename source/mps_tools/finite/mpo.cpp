//
// Created by david on 2019-06-25.
//
#include <mps_state/class_finite_chain_state.h>
#include <mps_tools/nmspc_mps_tools.h>
#include <model/class_hamiltonian_base.h>
#include <model/class_hamiltonian_factory.h>


void mpstools::finite::mpo::initialize(class_finite_chain_state & state, const size_t length, std::string model_type){
    log->info("Initializing mpo");
    //Generate MPO
    size_t pos = 0;
    state.MPO_L.emplace_back(class_hamiltonian_factory::create_mpo(pos++,model_type));
    while(true){
        state.MPO_R.emplace_back(class_hamiltonian_factory::create_mpo(pos++,model_type));
        if(state.MPO_L.size() + state.MPO_R.size() >= length){break;}
    }
}


void mpstools::finite::mpo::randomize(class_finite_chain_state &state) {
    log->info("Setting random fields in state");
    std::vector<std::vector<double>> all_params;
    for (auto &mpo : state.MPO_L){
        mpo->randomize_hamiltonian();
        all_params.push_back(mpo->get_parameter_values());
    }
    for (auto &mpo : state.MPO_R){
        mpo->randomize_hamiltonian();
        all_params.push_back(mpo->get_parameter_values());
    }

    for (auto &mpo : state.MPO_L){
        mpo->set_full_lattice_parameters(all_params);
    }
    for (auto &mpo : state.MPO_R){
        mpo->set_full_lattice_parameters(all_params);
    }
}