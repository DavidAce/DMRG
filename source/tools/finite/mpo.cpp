//
// Created by david on 2019-06-25.
//
#include <state/class_finite_state.h>
#include <tools/nmspc_tools.h>
#include <model/class_model_base.h>
#include <model/class_model_factory.h>


void tools::finite::mpo::initialize(class_finite_state & state, const size_t length, std::string model_type){
    log->info("Initializing mpo");
    //Generate MPO
    size_t pos = 0;
    state.MPO_L.emplace_back(class_model_factory::create_mpo(pos++,model_type));
    while(true){
        state.MPO_R.emplace_back(class_model_factory::create_mpo(pos++,model_type));
        if(state.MPO_L.size() + state.MPO_R.size() >= length){break;}
    }
}


void tools::finite::mpo::randomize(class_finite_state &state) {
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

    //

    for (auto &mpo : state.MPO_L){
        mpo->set_full_lattice_parameters(all_params,true);
    }
    for (auto &mpo : state.MPO_R){
        mpo->set_full_lattice_parameters(all_params,true);
    }
}
//
//MPO #           J_rnd           h_rnd           J_log_mean      h_log_mean      J_avg           h_avg           J_sigma         h_sigma         lambda          delta           e_reduced       spin_dim
//0               1.78029073      1.857672865     1               1               6.786844177     3.981962353     1               1               0               0               0               2
//1               0.8564592888    4.862651436     1               1               6.786844177     3.981962353     1               1               0               0               0               2
//2               67.10952028     6.259807039     1               1               6.786844177     3.981962353     1               1               0               0               0               2
//3               3.970975987     1.890555937     1               1               6.786844177     3.981962353     1               1               0               0               0               2
//4               2.175348662     11.12155811     1               1               6.786844177     3.981962353     1               1               0               0               0               2
//5               1.44676615      0.5403342528    1               1               6.786844177     3.981962353     1               1               0               0               0               2
//6               2.214192912     1.421169278     1               1               6.786844177     3.981962353     1               1               0               0               0               2
//7               1.036456611     5.217373846     1               1               6.786844177     3.981962353     1               1               0               0               0               2
//8               3.318821342     4.510961248     1               1               6.786844177     3.981962353     1               1               0               0               0               2
//9               7.14590453      1.21947362      1               1               6.786844177     3.981962353     1               1               0               0               0               2
//10              1.277359344     0.9150799398    1               1               6.786844177     3.981962353     1               1               0               0               0               2
//11              1.374393995     0.390866565     1               1               6.786844177     3.981962353     1               1               0               0               0               2
//12              2.822960709     0.3557954012    1               1               6.786844177     3.981962353     1               1               0               0               0               2
//13              2.549087719     12.74857849     1               1               6.786844177     3.981962353     1               1               0               0               0               2
//14              2.724124404     6.417557257     1               1               6.786844177     3.981962353     1               1               0               0               0               2
//15              0               1.403051852     1               1               6.786844177     3.981962353     1               1               0               0               0               2
