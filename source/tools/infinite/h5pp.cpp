//
// Created by david on 2019-03-09.
//

#include <tools/nmspc_tools.h>
#include <state/class_infinite_state.h>
#include <state/class_vidal_site.h>
#include <state/class_mps_2site.h>
#include <state/class_environment.h>
#include <model/class_model_base.h>
#include <h5pp/h5pp.h>

void tools::infinite::io::write_all_state(const class_infinite_state &state, h5pp::File &h5ppFile,
                                          std::string sim_name){
    if(state.has_been_written){return;}


    if (settings::hdf5::storage_level >= StorageLevel::NONE) {
    }

    if (settings::hdf5::storage_level >= StorageLevel::LIGHT) {
        write_hamiltonian_params(state,h5ppFile,sim_name);
    }

    if (settings::hdf5::storage_level >= StorageLevel::NORMAL) {
        write_2site_mps(state,h5ppFile,sim_name);
        write_2site_mpo(state,h5ppFile,sim_name);
        write_2site_env(state,h5ppFile,sim_name);
        write_2site_env2(state,h5ppFile,sim_name);
    }

    if (settings::hdf5::storage_level >= StorageLevel::FULL) {
    }

    state.has_been_written = true;
}


void tools::infinite::io::write_2site_mps (const class_infinite_state &state, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(state.MPS->MPS_A->get_A(), sim_name + "/state/2site/MPS_A");
    h5ppFile.writeDataset(state.MPS->LC                                 , sim_name + "/state/2site/L_C");
    h5ppFile.writeDataset(state.MPS->MPS_B->get_B(), sim_name + "/state/2site/MPS_B");
}

void tools::infinite::io::write_2site_mpo (const class_infinite_state &state, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(state.HA->MPO(), sim_name + "/state/2site/MPO_A");
    h5ppFile.writeDataset(state.HB->MPO(), sim_name + "/state/2site/MPO_B");

}

void tools::infinite::io::write_2site_env (const class_infinite_state &state, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(state.Lblock->block, sim_name + "/state/2site/ENV_L");
    h5ppFile.writeDataset(state.Rblock->block, sim_name + "/state/2site/ENV_R");
}

void tools::infinite::io::write_2site_env2 (const class_infinite_state &state, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(state.Lblock2->block, sim_name + "/state/2site/ENV2_L");
    h5ppFile.writeDataset(state.Rblock2->block, sim_name + "/state/2site/ENV2_R");
}

void tools::infinite::io::write_hamiltonian_params(const class_infinite_state &state, h5pp::File & h5ppFile, std::string sim_name){
    // Write down the Hamiltonian metadata as a table
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    Eigen::MatrixXd hamiltonian_props;

    auto propsA = state.HA->get_parameter_values();
    auto propsB = state.HB->get_parameter_values();
    Eigen::ArrayXd  temp_rowA  = Eigen::Map<Eigen::ArrayXd> (propsA.data(),propsA.size());
    Eigen::ArrayXd  temp_rowB  = Eigen::Map<Eigen::ArrayXd> (propsB.data(),propsB.size());
    hamiltonian_props.resize(2, temp_rowA.size());
    hamiltonian_props.row(0) = temp_rowA.transpose();
    hamiltonian_props.row(1) = temp_rowB.transpose();
    h5ppFile.writeDataset(hamiltonian_props ,sim_name + "/model/2site/Hamiltonian");
    int col = 0;
    for (auto &name : state.HA->get_parameter_names()){
        std::string attr_value = name;
        std::string attr_name  = "FIELD_" + std::to_string(col) + "_NAME";
        h5ppFile.writeAttributeToLink(attr_value, attr_name,sim_name + "/model/2site/Hamiltonian" );
        col++;
    }
}

void tools::infinite::io::write_all_measurements  (const class_infinite_state & state, h5pp::File & h5ppFile, std::string sim_name){
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

