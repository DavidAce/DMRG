//
// Created by david on 2019-06-25.
//
#include <mps_state/class_finite_chain_state.h>
#include <mps_tools/nmspc_mps_tools.h>
#include <model/class_hamiltonian_base.h>
#include <model/class_hamiltonian_factory.h>


void mpstools::finite::mpo::initialize_mpo(class_finite_chain_state & state, std::string model_type, const size_t length){
    //Generate MPO
    while(true){
        state.MPO_L.emplace_back(class_hamiltonian_factory::create_mpo(model_type));
        if(state.MPO_L.size() + state.MPO_R.size() >= length){break;}
        state.MPO_R.emplace_front(class_hamiltonian_factory::create_mpo(model_type));
        if(state.MPO_L.size() + state.MPO_R.size() >= length){break;}
    }
    mpstools::finite::mpo::randomize_mpo(state);
}


void mpstools::finite::mpo::randomize_mpo(class_finite_chain_state &state) {
    mpstools::log->info("Setting random fields in state");
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