//
// Created by david on 2019-03-09.
//

#include <tools/nmspc_tools.h>
#include <state/class_state_infinite.h>
#include <state/class_mps_site.h>
#include <state/class_mps_2site.h>
#include <state/class_environment.h>
#include <model/class_model_base.h>
#include <h5pp/h5pp.h>

void tools::infinite::io::h5dset::write_all_state(const class_state_infinite &state, h5pp::File &h5ppFile,
                                          std::string sim_name){
    if(state.has_been_written){return;}


    if (settings::output::storage_level >= StorageLevel::NONE) {
    }

    if (settings::output::storage_level >= StorageLevel::LIGHT) {
        write_hamiltonian_params(state,h5ppFile,sim_name);
    }

    if (settings::output::storage_level >= StorageLevel::NORMAL) {
        write_2site_mps(state,h5ppFile,sim_name);
        write_2site_mpo(state,h5ppFile,sim_name);
        write_2site_env(state,h5ppFile,sim_name);
        write_2site_env2(state,h5ppFile,sim_name);
    }

    if (settings::output::storage_level >= StorageLevel::FULL) {
    }

    state.has_been_written = true;
}


void tools::infinite::io::h5dset::write_2site_mps (const class_state_infinite &state, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(state.MPS->MPS_A->get_M() , sim_name + "/state/2site/MPS_A");
    h5ppFile.writeDataset(state.MPS->MPS_A->get_LC(), sim_name + "/state/2site/L_C");
    h5ppFile.writeDataset(state.MPS->MPS_B->get_M(), sim_name + "/state/2site/MPS_B");
}

void tools::infinite::io::h5dset::write_2site_mpo (const class_state_infinite &state, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(state.HA->MPO(), sim_name + "/state/2site/MPO_A");
    h5ppFile.writeDataset(state.HB->MPO(), sim_name + "/state/2site/MPO_B");

}

void tools::infinite::io::h5dset::write_2site_env (const class_state_infinite &state, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(state.Lblock->block, sim_name + "/state/2site/ENV_L");
    h5ppFile.writeDataset(state.Rblock->block, sim_name + "/state/2site/ENV_R");
}

void tools::infinite::io::h5dset::write_2site_env2 (const class_state_infinite &state, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(state.Lblock2->block, sim_name + "/state/2site/ENV2_L");
    h5ppFile.writeDataset(state.Rblock2->block, sim_name + "/state/2site/ENV2_R");
}

void tools::infinite::io::h5dset::write_hamiltonian_params(const class_state_infinite &state, h5pp::File & h5ppFile, std::string sim_name){

    h5ppFile.writeDataset(settings::model::model_type, sim_name + "/model/model_type");
    auto paramsA = state.HA->get_parameters();
    auto paramsB = state.HB->get_parameters();
        // Write MPO properties as attributes
    std::string dataset_name = sim_name + "/model/hamiltonian_A";
    h5ppFile.writeDataset("A", dataset_name);
    for(auto &params : state.HA->get_parameters()) {
        if(params.second.type() == typeid(double)) h5ppFile.writeAttribute(std::any_cast<double>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(size_t)) h5ppFile.writeAttribute(std::any_cast<size_t>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(int)) h5ppFile.writeAttribute(std::any_cast<int>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(bool)) h5ppFile.writeAttribute(std::any_cast<bool>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(std::string)) h5ppFile.writeAttribute(std::any_cast<std::string>(params.second), params.first, dataset_name);
    }
    dataset_name = sim_name + "/model/hamiltonian_B";
    h5ppFile.writeDataset("B", dataset_name);
    for(auto &params : state.HB->get_parameters()) {
        if(params.second.type() == typeid(double)) h5ppFile.writeAttribute(std::any_cast<double>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(size_t)) h5ppFile.writeAttribute(std::any_cast<size_t>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(int)) h5ppFile.writeAttribute(std::any_cast<int>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(bool)) h5ppFile.writeAttribute(std::any_cast<bool>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(std::string)) h5ppFile.writeAttribute(std::any_cast<std::string>(params.second), params.first, dataset_name);
    }

}

void tools::infinite::io::h5dset::write_all_measurements  (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name){
    state.do_all_measurements();
    h5ppFile.writeDataset(state.measurements.length.value()                      , sim_name + "/measurements/2site/length");
    h5ppFile.writeDataset(state.measurements.bond_dimension.value()              , sim_name + "/measurements/2site/bond_dimension");
    h5ppFile.writeDataset(state.measurements.norm.value()                        , sim_name + "/measurements/2site/norm");
    h5ppFile.writeDataset(state.measurements.truncation_error.value()            , sim_name + "/measurements/2site/truncation_error");
    h5ppFile.writeDataset(state.measurements.energy_mpo.value()                  , sim_name + "/measurements/2site/energy");
    h5ppFile.writeDataset(state.measurements.energy_per_site_mpo.value()         , sim_name + "/measurements/2site/energy_per_site");
    h5ppFile.writeDataset(state.measurements.energy_per_site_ham.value()         , sim_name + "/measurements/2site/energy_per_site_mom");
    h5ppFile.writeDataset(state.measurements.energy_per_site_mom.value()         , sim_name + "/measurements/2site/energy_per_site_mom");
    h5ppFile.writeDataset(state.measurements.energy_variance_per_site_mpo.value(), sim_name + "/measurements/2site/energy_variance_per_site");
    h5ppFile.writeDataset(state.measurements.energy_variance_per_site_ham.value(), sim_name + "/measurements/2site/energy_variance_per_site_ham");
    h5ppFile.writeDataset(state.measurements.energy_variance_per_site_mom.value(), sim_name + "/measurements/2site/energy_variance_per_site_mom");
    h5ppFile.writeDataset(state.measurements.current_entanglement_entropy.value(), sim_name + "/measurements/2site/entanglement_entropy_midchain");
}

