//
// Created by david on 2019-01-30.
//

//
// Created by david on 2018-12-06.
//

#include <mps_routines/nmspc_mps_tools.h>
#include <mps_routines/class_finite_chain_state.h>
#include <mps_routines/class_superblock.h>
#include <IO/class_hdf5_file.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <spdlog/spdlog.h>
#include <stdexcept>

using Scalar    = std::complex<double>;



void MPS_Tools::Finite::Hdf5::write_all_state(class_finite_chain_state &state, class_hdf5_file &hdf5,
                                              std::string sim_name) {
    if(state.has_been_written()){return;}
    write_2site_mps(state,hdf5,sim_name);
    write_2site_mpo(state,hdf5,sim_name);
    write_2site_env(state,hdf5,sim_name);
    write_2site_env2(state,hdf5,sim_name);
    write_bond_matrices(state,hdf5,sim_name);
    write_hamiltonian_params(state,hdf5,sim_name);

    if (settings::hdf5::full_storage){
        write_full_mps(state,hdf5,sim_name);
        write_full_mpo(state,hdf5,sim_name);
    }
    state.set_written_true();
}




void MPS_Tools::Finite::Hdf5::write_2site_mps(class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name)
/*! Writes down the local two-site MPS in "A - LC - B" notation.
 */
{
    // Write MPS in A-B notation.
    // Also write down all the center Lambdas (singular values) , so that we can obtain the entanglement spectrum easily.'
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    hdf5.write_dataset(Textra::to_RowMajor(state.get_MPS_L().back().get_A()), sim_name + "/state/2site/MPS_A");
    hdf5.write_dataset(state.get_MPS_C(), sim_name + "/state/2site/L_C");
    hdf5.write_dataset(Textra::to_RowMajor(state.get_MPS_R().front().get_B()),sim_name + "/state/2site/MPS_B");
}

void MPS_Tools::Finite::Hdf5::write_2site_mpo(class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name) {
    // Write two current the MPO's
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    auto &mpoL = state.get_MPO_L().back();
    auto &mpoR = state.get_MPO_R().front();
    hdf5.write_dataset(Textra::to_RowMajor(mpoL->MPO), sim_name + "/state/2site/MPO_L");
    hdf5.write_dataset(Textra::to_RowMajor(mpoR->MPO), sim_name + "/state/2site/MPO_R");
    auto valuesL = mpoL->get_parameter_values();
    auto valuesR = mpoR->get_parameter_values();
    auto namesL  = mpoL->get_parameter_names();
    auto namesR  = mpoR->get_parameter_names();

    for (size_t i = 0; i < std::min(valuesL.size(), namesL.size()); i++){
        hdf5.write_attribute_to_dataset(sim_name + "/state/2site/MPO_L", valuesL[i], namesL[i]);
    }
    for (size_t i = 0; i < std::min(valuesR.size(), namesR.size()); i++){
        hdf5.write_attribute_to_dataset(sim_name + "/state/2site/MPO_R", valuesR[i], namesR[i]);
    }
}

void MPS_Tools::Finite::Hdf5::write_2site_env(class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name){
    // Write the environment blocks
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    hdf5.write_dataset(Textra::to_RowMajor(state.get_ENV_L().back().block), sim_name + "/state/2site/ENV_L");
    hdf5.write_attribute_to_dataset(sim_name + "/state/2site/ENV_L", state.get_ENV_L().back().size, "sites");
    hdf5.write_dataset(Textra::to_RowMajor(state.get_ENV_R().front().block), sim_name + "/state/2site/ENV_R");
    hdf5.write_attribute_to_dataset(sim_name + "/state/2site/ENV_R", state.get_ENV_R().front().size, "sites");
}


void MPS_Tools::Finite::Hdf5::write_2site_env2(class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name){
    // Write the environment squared blocks
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    hdf5.write_dataset(Textra::to_RowMajor(state.get_ENV2_L().back().block), sim_name + "/state/2site/ENV2_L");
    hdf5.write_attribute_to_dataset(sim_name + "/state/2site/ENV2_L", state.get_ENV2_L().back().size, "sites");
    hdf5.write_dataset(Textra::to_RowMajor(state.get_ENV2_R().front().block), sim_name + "/state/2site/ENV2_R");
    hdf5.write_attribute_to_dataset(sim_name + "/state/2site/ENV2_R", state.get_ENV2_R().front().size, "sites");
}

void MPS_Tools::Finite::Hdf5::write_bond_matrices(class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name)
/*! Writes down all the "Lambda" bond matrices (singular value matrices), so we can obtain the entanglement spectrum easily.
 */
{
    unsigned long counter = 0;
    for (auto &mps : state.get_MPS_L()){
        hdf5.write_dataset(mps.get_L()                     ,sim_name + "/state/full/bonds/L_" + std::to_string(counter++));
    }
    hdf5.write_dataset(state.get_MPS_C(),sim_name + "/state/full/bonds/L_C");
    for (auto &mps : state.get_MPS_R()){
        hdf5.write_dataset(mps.get_L()                      ,sim_name + "/state/full/bonds/L_" + std::to_string(counter++));
    }
}

void MPS_Tools::Finite::Hdf5::write_full_mps(class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name)
/*!
 * Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation.
 *
 */
{
    if(state.mps_have_been_written_to_hdf5){return;}
    unsigned long counter = 0;
    for (auto &mps : state.get_MPS_L()){
        hdf5.write_dataset(Textra::to_RowMajor(mps.get_G()),sim_name + "/state/full/mps/G_" + std::to_string(counter));
        hdf5.write_dataset(mps.get_L()                     ,sim_name + "/state/full/mps/L_" + std::to_string(counter++));
    }

    hdf5.write_dataset(state.get_MPS_C(), sim_name + "/state/full/mps/L_C");
    for (auto &mps : state.get_MPS_R()){
        hdf5.write_dataset(Textra::to_RowMajor(mps.get_G()) ,sim_name + "/state/full/mps/G_" + std::to_string(counter));
        hdf5.write_dataset(mps.get_L()                      ,sim_name + "/state/full/mps/L_" + std::to_string(counter++));
    }
    state.mps_have_been_written_to_hdf5 = true;
}




void MPS_Tools::Finite::Hdf5::write_full_mpo(class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name) {
    // Write all the MPO's
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    if(state.mpo_have_been_written_to_hdf5){return;}
    unsigned long counter = 0;
    for(auto &mpo : state.get_MPO_L()){
        hdf5.write_dataset(Textra::to_RowMajor(mpo->MPO), sim_name + "/state/full/mpo/H_" + std::to_string(counter));
        //Write MPO properties as attributes
        auto values = mpo->get_parameter_values();
        auto names  = mpo->get_parameter_names();
        for (size_t i = 0; i < std::min(values.size(), names.size()); i++){
            hdf5.write_attribute_to_dataset(sim_name + "/state/full/mpo/H_" + std::to_string(counter), values[i], names[i]);
        }
        counter++;
    }
    for(auto &mpo : state.get_MPO_R()){
        hdf5.write_dataset(Textra::to_RowMajor(mpo->MPO), sim_name + "/state/full/mpo/H_" + std::to_string(counter));
        //Write MPO properties as attributes
        auto values = mpo->get_parameter_values();
        auto names  = mpo->get_parameter_names();
        for (size_t i = 0; i < std::min(values.size(), names.size()); i++){
            hdf5.write_attribute_to_dataset(sim_name + "/state/full/mpo/H_" + std::to_string(counter), values[i], names[i]);
        }
        counter++;
    }
    state.mpo_have_been_written_to_hdf5 = true;
}



void MPS_Tools::Finite::Hdf5::write_hamiltonian_params(class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name){
    // Write down the Hamiltonian metadata as a table
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    Eigen::MatrixXd hamiltonian_props;
    for(auto &mpo : state.get_MPO_L()){
        auto props = mpo->get_parameter_values();
        Eigen::ArrayXd  temp_row  = Eigen::Map<Eigen::ArrayXd> (props.data(),props.size());
        hamiltonian_props.conservativeResize(hamiltonian_props.rows()+1, temp_row.size());
        hamiltonian_props.bottomRows(1) = temp_row.transpose();
    }
    for(auto &mpo : state.get_MPO_R()){
        auto props = mpo->get_parameter_values();
        Eigen::ArrayXd  temp_row  = Eigen::Map<Eigen::ArrayXd> (props.data(),props.size());
        hamiltonian_props.conservativeResize(hamiltonian_props.rows()+1, temp_row.size());
        hamiltonian_props.bottomRows(1) = temp_row.transpose();
    }
    hdf5.write_dataset(Textra::to_RowMajor(hamiltonian_props.transpose()) ,sim_name + "/model/Hamiltonian");
    int col = 0;
    for (auto &name : state.get_MPO_L().front()->get_parameter_names()){
        std::string attr_value = name;
        std::string attr_name  = "FIELD_" + std::to_string(col) + "_NAME";
        hdf5.write_attribute_to_dataset(sim_name + "/model/Hamiltonian", attr_value, attr_name );
        col++;
    }
}


void MPS_Tools::Finite::Hdf5::write_all_measurements(class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name){
    state.do_all_measurements();
    hdf5.write_dataset(state.measurements.length                      , sim_name + "/measurements/length");
    hdf5.write_dataset(state.measurements.norm                        , sim_name + "/measurements/norm");
    hdf5.write_dataset(state.measurements.bond_dimensions             , sim_name + "/measurements/bond_dimensions");
    hdf5.write_dataset(state.measurements.energy_per_site_mpo         , sim_name + "/measurements/energy_per_site_mpo");
    hdf5.write_dataset(state.measurements.energy_variance_per_site_mpo, sim_name + "/measurements/energy_variance_per_site_mpo");
    hdf5.write_dataset(state.measurements.entanglement_entropies      , sim_name + "/measurements/entanglement_entropies");
    hdf5.write_dataset(state.measurements.spin_components             , sim_name + "/measurements/spin_components");
    hdf5.write_dataset(state.measurements.spin_component_sx           , sim_name + "/measurements/spin_component_sx");
    hdf5.write_dataset(state.measurements.spin_component_sy           , sim_name + "/measurements/spin_component_sy");
    hdf5.write_dataset(state.measurements.spin_component_sz           , sim_name + "/measurements/spin_component_sz");
}



void MPS_Tools::Finite::Hdf5::load_state_from_hdf5(class_finite_chain_state & state, class_superblock & superblock, class_hdf5_file & hdf5, std::string sim_name){
    // Load into state
    state.clear();
    superblock.clear();

    if (not hdf5.link_exists(sim_name + "/state/full/mps")){
        spdlog::debug("Link does not exist: {}"  ,sim_name + "/state/full/mps");
        throw std::logic_error("Tried to load non-existing MPS");
    }

    auto mps_list = hdf5.print_contents_of_group(sim_name + "/state/full/mps");
    auto mpo_list = hdf5.print_contents_of_group(sim_name + "/state/full/mpo");


    for (auto &mps : mps_list){
        std::cout << mps << std::endl;
    }
    std::cerr << "HAVENT IMPLEMENTED LOADING" << std::endl;

    exit(0);
//
//
//    int chain_length =

}

