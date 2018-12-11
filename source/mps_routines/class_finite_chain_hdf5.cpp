//
// Created by david on 2018-12-06.
//

#include "class_finite_chain.h"


void class_finite_chain::set_hdf5_file(std::shared_ptr<class_hdf5_file> hdf5_){
    hdf5 = std::move(hdf5_);
    hdf5_file_is_set = true;
}

void class_finite_chain::write_all_to_hdf5() {

    write_2site_mps_to_hdf5();
    write_2site_mpo_to_hdf5();
    write_2site_env_to_hdf5();
    write_2site_env2_to_hdf5();

    write_bond_matrices_to_hdf5();
    write_hamiltonian_params_to_hdf5();
    write_entanglement_to_hdf5();

    if (settings::hdf5::full_storage){
        write_full_mps_to_hdf5();
        write_full_mpo_to_hdf5();
    }

}


void class_finite_chain::write_bond_matrices_to_hdf5()
/*! Writes down all the "Lambda" bond matrices (singular value matrices), so we can obtain the entanglement spectrum easily.
 */
 {
    unsigned long counter = 0;
    for (auto &mps : get_MPS_L()){
        hdf5->write_dataset(mps.get_L()                     ,sim_name + "/chain/bonds/L_" + std::to_string(counter++));
    }
    hdf5->write_dataset(get_MPS_C(),sim_name + "/chain/bonds/L_" + std::to_string(counter++) + "_C");
    for (auto &mps : get_MPS_R()){
        hdf5->write_dataset(mps.get_L()                      ,sim_name + "/chain/bonds/L_" + std::to_string(counter++));
    }
}

void class_finite_chain::write_2site_mps_to_hdf5()
/*! Writes down the local two-site MPS in "A - LC - B" notation.
 */
{
    // Write MPS in A-B notation.
    // Also write down all the center Lambdas (singular values) , so that we can obtain the entanglement spectrum easily.'
    // Remember to write tensors in row-major storage order because that's what hdf5 uses.
    hdf5->write_dataset(Textra::to_RowMajor(get_MPS_L().back().get_A()), sim_name + "/2site/MPS_A");
    hdf5->write_dataset(get_MPS_C(), sim_name + "/2site/L_C");
    hdf5->write_dataset(Textra::to_RowMajor(get_MPS_R().front().get_B()),sim_name + "/2site/MPS_B");
}

void class_finite_chain::write_2site_mpo_to_hdf5() {
    // Write two current the MPO's
    // Remember to write tensors in row-major storage order because that's what hdf5 uses.
    auto &mpoL = get_MPO_L().back();
    auto &mpoR = get_MPO_R().front();
    hdf5->write_dataset(Textra::to_RowMajor(mpoL->MPO), sim_name + "/2site/MPO_L");
    hdf5->write_dataset(Textra::to_RowMajor(mpoR->MPO), sim_name + "/2site/MPO_R");
    auto valuesL = mpoL->get_parameter_values();
    auto valuesR = mpoR->get_parameter_values();
    auto namesL  = mpoL->get_parameter_names();
    auto namesR  = mpoR->get_parameter_names();

    for (size_t i = 0; i < std::min(valuesL.size(), namesL.size()); i++){
        hdf5->write_attribute_to_dataset(sim_name + "/2site/MPO_L", valuesL[i], namesL[i]);
    }
    for (size_t i = 0; i < std::min(valuesR.size(), namesR.size()); i++){
        hdf5->write_attribute_to_dataset(sim_name + "/2site/MPO_R", valuesR[i], namesR[i]);
    }
}

void class_finite_chain::write_2site_env_to_hdf5(){
    // Write the environment blocks
    // Remember to write tensors in row-major storage order because that's what hdf5 uses.
    hdf5->write_dataset(Textra::to_RowMajor(get_ENV_L().back().block), sim_name + "/2site/ENV_L");
    hdf5->write_attribute_to_dataset(sim_name + "/2site/ENV_L", get_ENV_L().back().size, "sites");
    hdf5->write_dataset(Textra::to_RowMajor(get_ENV_R().front().block), sim_name + "/2site/ENV_R");
    hdf5->write_attribute_to_dataset(sim_name + "/2site/ENV_R", get_ENV_R().front().size, "sites");
}


void class_finite_chain::write_2site_env2_to_hdf5(){
    // Write the environment squared blocks
    // Remember to write tensors in row-major storage order because that's what hdf5 uses.
    hdf5->write_dataset(Textra::to_RowMajor(get_ENV2_L().back().block), sim_name + "/2site/ENV2_L");
    hdf5->write_attribute_to_dataset(sim_name + "/2site/ENV2_L", get_ENV2_L().back().size, "sites");
    hdf5->write_dataset(Textra::to_RowMajor(get_ENV2_R().front().block), sim_name + "/2site/ENV2_R");
    hdf5->write_attribute_to_dataset(sim_name + "/2site/ENV2_R", get_ENV2_R().front().size, "sites");
}



void class_finite_chain::write_full_mps_to_hdf5()
/*! Writes down the full MPS in "A-A-A...A-LC-B...B-B-B" notation.
 */
{
    if(full_mps_has_been_written){return;}
    unsigned long counter = 0;
    // Write MPS in A-B notation for simplicity.
    // Also write down all the Lambdas (singular value) , so that we can obtain the entanglement spectrum easily.'
    // Remember to write tensors in row-major storage order because that's what hdf5 uses.
    for (auto &mps : get_MPS_L()){
        hdf5->write_dataset(Textra::to_RowMajor(mps.get_A()),sim_name + "/chain/mps/A_" + std::to_string(counter));
    }
    hdf5->write_dataset(get_MPS_C(), sim_name + "/chain/mps/L_C");
    for (auto &mps : get_MPS_R()){
        hdf5->write_dataset(Textra::to_RowMajor(mps.get_B()) ,sim_name + "/chain/mps/B_" + std::to_string(counter));
    }
    full_mps_has_been_written = true;
}




void class_finite_chain::write_full_mpo_to_hdf5() {
    // Write all the MPO's
    // Remember to write tensors in row-major storage order because that's what hdf5 uses.
    if(full_mpo_has_been_written){return;}
    unsigned long counter = 0;
    for(auto &mpo : get_MPO_L()){
        hdf5->write_dataset(Textra::to_RowMajor(mpo->MPO), sim_name + "/chain/mpo/H_" + std::to_string(counter));
        //Write MPO properties as attributes
        auto values = mpo->get_parameter_values();
        auto names  = mpo->get_parameter_names();
        for (size_t i = 0; i < std::min(values.size(), names.size()); i++){
            hdf5->write_attribute_to_dataset(sim_name + "/chain/mpo/H_" + std::to_string(counter), values[i], names[i]);
        }
        counter++;
    }
    for(auto &mpo : get_MPO_R()){
        hdf5->write_dataset(Textra::to_RowMajor(mpo->MPO), sim_name + "/chain/mpo/H_" + std::to_string(counter));
        //Write MPO properties as attributes
        auto values = mpo->get_parameter_values();
        auto names  = mpo->get_parameter_names();
        for (size_t i = 0; i < std::min(values.size(), names.size()); i++){
            hdf5->write_attribute_to_dataset(sim_name + "/chain/mpo/H_" + std::to_string(counter), values[i], names[i]);
        }
        counter++;
    }
    full_mpo_has_been_written = true;
}



void class_finite_chain::write_hamiltonian_params_to_hdf5(){
    // Write down the Hamiltonian metadata as a table
    // Remember to write tensors in row-major storage order because that's what hdf5 uses.
    Eigen::MatrixXd hamiltonian_props;
    for(auto &mpo : get_MPO_L()){
        auto props = mpo->get_parameter_values();
        Eigen::ArrayXd  temp_row  = Eigen::Map<Eigen::ArrayXd> (props.data(),props.size());
        hamiltonian_props.conservativeResize(hamiltonian_props.rows()+1, temp_row.size());
        hamiltonian_props.bottomRows(1) = temp_row.transpose();
    }
    for(auto &mpo : get_MPO_R()){
        auto props = mpo->get_parameter_values();
        Eigen::ArrayXd  temp_row  = Eigen::Map<Eigen::ArrayXd> (props.data(),props.size());
        hamiltonian_props.conservativeResize(hamiltonian_props.rows()+1, temp_row.size());
        hamiltonian_props.bottomRows(1) = temp_row.transpose();
    }
    hdf5->write_dataset(Textra::to_RowMajor(hamiltonian_props.transpose()) ,sim_name + "/Hamiltonian");
    int col = 0;
    for (auto &name : get_MPO_L().front()->get_parameter_names()){
        std::string attr_value = name;
        std::string attr_name  = "FIELD_" + std::to_string(col) + "_NAME";
        hdf5->write_attribute_to_dataset(sim_name + "/Hamiltonian", attr_value, attr_name );
        col++;
    }
}



void class_finite_chain::write_entanglement_to_hdf5(){
    // Write relevant quantities
    std::vector<double> entanglement_entropies;
    for (auto &mps : get_MPS_L()){
        Eigen::Tensor<Scalar,0> SA  = -mps.get_L().square()
                .contract(mps.get_L().square().log().eval(), Textra::idx({0},{0}));
        entanglement_entropies.push_back(std::real(SA(0)));
    }
    Eigen::Tensor<Scalar,0> SA  = -get_MPS_C().square()
            .contract(get_MPS_C().square().log().eval(), Textra::idx({0},{0}));
    entanglement_entropies.push_back(std::real(SA(0)));
    for (auto &mps : get_MPS_R()){
        Eigen::Tensor<Scalar,0> SA  = -mps.get_L().square()
                .contract(mps.get_L().square().log().eval(), Textra::idx({0},{0}));
        entanglement_entropies.push_back(std::real(SA(0)));
    }
    hdf5->write_dataset(entanglement_entropies ,sim_name + "/entanglement_entropy");

}