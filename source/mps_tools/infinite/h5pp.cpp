//
// Created by david on 2019-03-09.
//

#include <mps_tools/nmspc_mps_tools.h>
#include <mps_state/class_superblock.h>
#include <mps_state/class_mps_2site.h>
#include <mps_state/class_environment.h>
#include <model/class_hamiltonian_base.h>
#include <h5pp/h5pp.h>

void mpstools::infinite::io::write_all_superblock(const class_superblock &superblock, h5pp::File & h5ppFile,
                                                     std::string sim_name){
    if(superblock.has_been_written){return;}


    if (settings::hdf5::storage_level >= StorageLevel::NONE) {
    }

    if (settings::hdf5::storage_level >= StorageLevel::LIGHT) {
        write_hamiltonian_params(superblock,h5ppFile,sim_name);
    }

    if (settings::hdf5::storage_level >= StorageLevel::NORMAL) {
        write_2site_mps(superblock,h5ppFile,sim_name);
        write_2site_mpo(superblock,h5ppFile,sim_name);
        write_2site_env(superblock,h5ppFile,sim_name);
        write_2site_env2(superblock,h5ppFile,sim_name);
    }

    if (settings::hdf5::storage_level >= StorageLevel::FULL) {
    }


        superblock.has_been_written = true;
}


void mpstools::infinite::io::write_2site_mps (const class_superblock &superblock, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(superblock.MPS->MPS_A->get_A(), sim_name + "/state/2site/MPS_A");
    h5ppFile.writeDataset(superblock.MPS->LC                                 , sim_name + "/state/2site/L_C");
    h5ppFile.writeDataset(superblock.MPS->MPS_B->get_B(), sim_name + "/state/2site/MPS_B");
}

void mpstools::infinite::io::write_2site_mpo (const class_superblock &superblock, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(superblock.HA->MPO(), sim_name + "/state/2site/MPO_A");
    h5ppFile.writeDataset(superblock.HB->MPO(), sim_name + "/state/2site/MPO_B");

}

void mpstools::infinite::io::write_2site_env (const class_superblock &superblock, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(superblock.Lblock->block, sim_name + "/state/2site/ENV_L");
    h5ppFile.writeDataset(superblock.Rblock->block, sim_name + "/state/2site/ENV_R");
}

void mpstools::infinite::io::write_2site_env2 (const class_superblock &superblock, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(superblock.Lblock2->block, sim_name + "/state/2site/ENV2_L");
    h5ppFile.writeDataset(superblock.Rblock2->block, sim_name + "/state/2site/ENV2_R");
}

void mpstools::infinite::io::write_hamiltonian_params(const class_superblock &superblock, h5pp::File & h5ppFile, std::string sim_name){
    // Write down the Hamiltonian metadata as a table
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    Eigen::MatrixXd hamiltonian_props;

    auto propsA = superblock.HA->get_parameter_values();
    auto propsB = superblock.HB->get_parameter_values();
    Eigen::ArrayXd  temp_rowA  = Eigen::Map<Eigen::ArrayXd> (propsA.data(),propsA.size());
    Eigen::ArrayXd  temp_rowB  = Eigen::Map<Eigen::ArrayXd> (propsB.data(),propsB.size());
    hamiltonian_props.resize(2, temp_rowA.size());
    hamiltonian_props.row(0) = temp_rowA.transpose();
    hamiltonian_props.row(1) = temp_rowB.transpose();
    h5ppFile.writeDataset(hamiltonian_props ,sim_name + "/model/2site/Hamiltonian");
    int col = 0;
    for (auto &name : superblock.HA->get_parameter_names()){
        std::string attr_value = name;
        std::string attr_name  = "FIELD_" + std::to_string(col) + "_NAME";
        h5ppFile.writeAttributeToLink(attr_value, attr_name,sim_name + "/model/2site/Hamiltonian" );
        col++;
    }
}

void mpstools::infinite::io::write_all_measurements  (const class_superblock & superblock, h5pp::File & h5ppFile, std::string sim_name){
    superblock.do_all_measurements();
    h5ppFile.writeDataset(superblock.measurements.length.value()                      , sim_name + "/measurements/2site/length");
    h5ppFile.writeDataset(superblock.measurements.bond_dimension.value()              , sim_name + "/measurements/2site/bond_dimension");
    h5ppFile.writeDataset(superblock.measurements.norm.value()                        , sim_name + "/measurements/2site/norm");
    h5ppFile.writeDataset(superblock.measurements.truncation_error.value()            , sim_name + "/measurements/2site/truncation_error");
    h5ppFile.writeDataset(superblock.measurements.energy_mpo.value()                  , sim_name + "/measurements/2site/energy");
    h5ppFile.writeDataset(superblock.measurements.energy_per_site_mpo.value()         , sim_name + "/measurements/2site/energy_per_site");
    h5ppFile.writeDataset(superblock.measurements.energy_per_site_ham.value()         , sim_name + "/measurements/2site/energy_per_site_mom");
    h5ppFile.writeDataset(superblock.measurements.energy_per_site_mom.value()         , sim_name + "/measurements/2site/energy_per_site_mom");
    h5ppFile.writeDataset(superblock.measurements.energy_variance_per_site_mpo.value(), sim_name + "/measurements/2site/energy_variance_per_site");
    h5ppFile.writeDataset(superblock.measurements.energy_variance_per_site_ham.value(), sim_name + "/measurements/2site/energy_variance_per_site_ham");
    h5ppFile.writeDataset(superblock.measurements.energy_variance_per_site_mom.value(), sim_name + "/measurements/2site/energy_variance_per_site_mom");
    h5ppFile.writeDataset(superblock.measurements.current_entanglement_entropy.value(), sim_name + "/measurements/2site/middle_entanglement_entropy");
}

void mpstools::infinite::io::load_from_hdf5    ([[maybe_unused]] const h5pp::File & h5ppFile, [[maybe_unused]] class_superblock & superblock,[[maybe_unused]] class_simulation_state &sim_state,[[maybe_unused]] std::string sim_name){

}
