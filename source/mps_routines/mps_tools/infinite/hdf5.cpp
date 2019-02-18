//
// Created by david on 2019-02-14.
//

#include <mps_routines/nmspc_mps_tools.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mps_2site.h>
#include <mps_routines/class_environment.h>
#include <model/class_hamiltonian_base.h>
#include <io/class_hdf5_file.h>

void MPS_Tools::Infinite::Hdf5::write_superblock_state (class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name){
    if(superblock.has_been_written){return;}
    write_2site_mps(superblock,hdf5,sim_name);
    write_2site_mpo(superblock,hdf5,sim_name);
    write_2site_env(superblock,hdf5,sim_name);
    write_2site_env2(superblock,hdf5,sim_name);
    write_hamiltonian_params(superblock,hdf5,sim_name);
    superblock.has_been_written = true;
}


void MPS_Tools::Infinite::Hdf5::write_2site_mps (class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name){
    hdf5.write_dataset(superblock.MPS->MPS_A->get_A(), sim_name + "/state/2site/MPS_A");
    hdf5.write_dataset(superblock.MPS->LC                                 , sim_name + "/state/2site/L_C");
    hdf5.write_dataset(superblock.MPS->MPS_B->get_B(), sim_name + "/state/2site/MPS_B");
}

void MPS_Tools::Infinite::Hdf5::write_2site_mpo (class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name){
    hdf5.write_dataset(superblock.HA->MPO, sim_name + "/state/2site/MPO_A");
    hdf5.write_dataset(superblock.HB->MPO, sim_name + "/state/2site/MPO_B");

}

void MPS_Tools::Infinite::Hdf5::write_2site_env (class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name){
    hdf5.write_dataset(superblock.Lblock->block, sim_name + "/state/2site/ENV_L");
    hdf5.write_dataset(superblock.Rblock->block, sim_name + "/state/2site/ENV_R");
}

void MPS_Tools::Infinite::Hdf5::write_2site_env2 (class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name){
    hdf5.write_dataset(superblock.Lblock2->block, sim_name + "/state/2site/ENV2_L");
    hdf5.write_dataset(superblock.Rblock2->block, sim_name + "/state/2site/ENV2_R");
}

void MPS_Tools::Infinite::Hdf5::write_hamiltonian_params(class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name){
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
    hdf5.write_dataset(hamiltonian_props ,sim_name + "/model/Hamiltonian");
    int col = 0;
    for (auto &name : superblock.HA->get_parameter_names()){
        std::string attr_value = name;
        std::string attr_name  = "FIELD_" + std::to_string(col) + "_NAME";
        hdf5.write_attribute_to_dataset(sim_name + "/model/Hamiltonian", attr_value, attr_name );
        col++;
    }
}

void MPS_Tools::Infinite::Hdf5::write_all_measurements  (class_superblock &superblock, class_hdf5_file & hdf5, std::string sim_name){
    superblock.do_all_measurements();
    hdf5.write_dataset(superblock.measurements.length                      , sim_name + "/measurements/length");
    hdf5.write_dataset(superblock.measurements.bond_dimension              , sim_name + "/measurements/bond_dimensions");
    hdf5.write_dataset(superblock.measurements.norm                        , sim_name + "/measurements/norm");
    hdf5.write_dataset(superblock.measurements.truncation_error            , sim_name + "/measurements/truncation_error");
    hdf5.write_dataset(superblock.measurements.energy_mpo                  , sim_name + "/measurements/energy_mpo");
    hdf5.write_dataset(superblock.measurements.energy_per_site_mpo         , sim_name + "/measurements/energy_per_site_mpo");
    hdf5.write_dataset(superblock.measurements.energy_per_site_ham         , sim_name + "/measurements/energy_per_site_mom");
    hdf5.write_dataset(superblock.measurements.energy_per_site_mom         , sim_name + "/measurements/energy_per_site_mom");
    hdf5.write_dataset(superblock.measurements.energy_variance_per_site_mpo, sim_name + "/measurements/energy_variance_per_site_mpo");
    hdf5.write_dataset(superblock.measurements.energy_variance_per_site_ham, sim_name + "/measurements/energy_variance_per_site_ham");
    hdf5.write_dataset(superblock.measurements.energy_variance_per_site_mom, sim_name + "/measurements/energy_variance_per_site_mom");
    hdf5.write_dataset(superblock.measurements.current_entanglement_entropy, sim_name + "/measurements/entanglement_entropy");
}

void MPS_Tools::Infinite::Hdf5::load_from_hdf5    (class_superblock &superblock,class_simulation_state &sim_state,class_hdf5_file & hdf5, std::string sim_name){

}
