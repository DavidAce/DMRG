//
// Created by david on 2019-03-09.
//


#include <mps_routines/nmspc_mps_tools.h>
#include <mps_routines/class_finite_chain_state.h>
#include <mps_routines/class_superblock.h>
#include <model/class_hamiltonian_factory.h>
#include <algorithms/class_simulation_state.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <spdlog/spdlog.h>
#include <stdexcept>
#include <h5pp/h5pp.h>

using Scalar    = std::complex<double>;



void MPS_Tools::Finite::H5pp::write_all_state(class_finite_chain_state &state, h5pp::File & h5ppFile,
                                              std::string sim_name) {
    if(state.has_been_written()){return;}

    if (settings::hdf5::storage_level >= StorageLevel::LIGHT){
        write_bond_matrices(state,h5ppFile,sim_name);
        write_hamiltonian_params(state,h5ppFile,sim_name);
    }
    if (settings::hdf5::storage_level >= StorageLevel::NORMAL){
        write_2site_mps(state,h5ppFile,sim_name);
        write_2site_mpo(state,h5ppFile,sim_name);
        write_2site_env(state,h5ppFile,sim_name);
        write_2site_env2(state,h5ppFile,sim_name);
    }


    if (settings::hdf5::storage_level >= StorageLevel::FULL){
        write_full_mps(state,h5ppFile,sim_name);
        write_full_mpo(state,h5ppFile,sim_name);
    }
    state.set_written_true();
}




void MPS_Tools::Finite::H5pp::write_2site_mps(class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name)
/*! Writes down the local two-site MPS in "A - LC - B" notation.
 */
{
    // Write MPS in A-B notation.
    // Also write down all the center Lambdas (singular values) , so that we can obtain the entanglement spectrum easily.'
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    h5ppFile.writeDataset(state.get_MPS_L().back().get_A(), sim_name + "/state/2site/MPS_A");
    h5ppFile.writeDataset(state.get_MPS_C(), sim_name + "/state/2site/L_C");
    h5ppFile.writeDataset(state.get_MPS_R().front().get_B(),sim_name + "/state/2site/MPS_B");
}

void MPS_Tools::Finite::H5pp::write_2site_mpo(class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name) {
    // Write two current the MPO's
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    auto &mpoL = state.get_MPO_L().back();
    auto &mpoR = state.get_MPO_R().front();
    h5ppFile.writeDataset(mpoL->MPO, sim_name + "/state/2site/MPO_L");
    h5ppFile.writeDataset(mpoR->MPO, sim_name + "/state/2site/MPO_R");
    auto valuesL = mpoL->get_parameter_values();
    auto valuesR = mpoR->get_parameter_values();
    auto namesL  = mpoL->get_parameter_names();
    auto namesR  = mpoR->get_parameter_names();

    for (size_t i = 0; i < std::min(valuesL.size(), namesL.size()); i++){
        h5ppFile.writeAttributeToLink(valuesL[i], namesL[i],sim_name + "/state/2site/MPO_L");
    }
    for (size_t i = 0; i < std::min(valuesR.size(), namesR.size()); i++){
        h5ppFile.writeAttributeToLink(valuesR[i], namesR[i],sim_name + "/state/2site/MPO_R");
    }
}

void MPS_Tools::Finite::H5pp::write_2site_env(class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name){
    // Write the environment blocks
    h5ppFile.writeDataset(state.get_ENV_L().back().block, sim_name + "/state/2site/ENV_L");
    h5ppFile.writeAttributeToLink(state.get_ENV_L().back().size, "sites",sim_name + "/state/2site/ENV_L");
    h5ppFile.writeDataset(state.get_ENV_R().front().block, sim_name + "/state/2site/ENV_R");
    h5ppFile.writeAttributeToLink(state.get_ENV_R().front().size, "sites",sim_name + "/state/2site/ENV_R");
}


void MPS_Tools::Finite::H5pp::write_2site_env2(class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name){
    // Write the environment squared blocks
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    h5ppFile.writeDataset(state.get_ENV2_L().back().block, sim_name + "/state/2site/ENV2_L");
    h5ppFile.writeAttributeToLink(state.get_ENV2_L().back().size, "sites",sim_name + "/state/2site/ENV2_L");
    h5ppFile.writeDataset(state.get_ENV2_R().front().block, sim_name + "/state/2site/ENV2_R");
    h5ppFile.writeAttributeToLink(state.get_ENV2_R().front().size, "sites",sim_name + "/state/2site/ENV2_R");
}

void MPS_Tools::Finite::H5pp::write_bond_matrices(class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name)
/*! Writes down all the "Lambda" bond matrices (singular value matrices), so we can obtain the entanglement spectrum easily.
 */
{
    unsigned long counter = 0;
    for (auto &mps : state.get_MPS_L()){
        h5ppFile.writeDataset(mps.get_L()                     ,sim_name + "/state/full/bonds/L_" + std::to_string(counter++));
    }
    h5ppFile.writeDataset(state.get_MPS_C(),sim_name + "/state/full/bonds/L_C");
    for (auto &mps : state.get_MPS_R()){
        h5ppFile.writeDataset(mps.get_L()                      ,sim_name + "/state/full/bonds/L_" + std::to_string(counter++));
    }
    h5ppFile.writeDataset(state.get_position() ,sim_name + "/state/full/position");
    h5ppFile.writeDataset(state.get_length()   ,sim_name + "/state/full/sites");
}

void MPS_Tools::Finite::H5pp::write_full_mps(class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name)
/*!
 * Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation.
 *
 */
{
    if(state.mps_have_been_written_to_hdf5){return;}
    unsigned long counter = 0;
    for (auto &mps : state.get_MPS_L()){
        h5ppFile.writeDataset(mps.get_G(),sim_name + "/state/full/mps/G_" + std::to_string(counter));
        h5ppFile.writeDataset(mps.get_L()                     ,sim_name + "/state/full/mps/L_" + std::to_string(counter++));
    }

    h5ppFile.writeDataset(state.get_MPS_C(), sim_name + "/state/full/mps/L_C");
    for (auto &mps : state.get_MPS_R()){
        h5ppFile.writeDataset(mps.get_G() ,sim_name + "/state/full/mps/G_" + std::to_string(counter));
        h5ppFile.writeDataset(mps.get_L()                      ,sim_name + "/state/full/mps/L_" + std::to_string(counter++));
    }
    h5ppFile.writeDataset(state.get_position() ,sim_name + "/state/full/position");
    h5ppFile.writeDataset(state.get_length()   ,sim_name + "/state/full/sites");
    state.mps_have_been_written_to_hdf5 = true;
}




void MPS_Tools::Finite::H5pp::write_full_mpo(class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name) {
    // Write all the MPO's
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    if(state.mpo_have_been_written_to_hdf5){return;}
    unsigned long counter = 0;
    for(auto &mpo : state.get_MPO_L()){
        h5ppFile.writeDataset(mpo->MPO, sim_name + "/state/full/mpo/H_" + std::to_string(counter));
        //Write MPO properties as attributes
        auto values = mpo->get_parameter_values();
        auto names  = mpo->get_parameter_names();
        for (size_t i = 0; i < std::min(values.size(), names.size()); i++){
            h5ppFile.writeAttributeToLink(values[i], names[i],sim_name + "/state/full/mpo/H_" + std::to_string(counter));
        }
        counter++;
    }
    for(auto &mpo : state.get_MPO_R()){
        h5ppFile.writeDataset(mpo->MPO, sim_name + "/state/full/mpo/H_" + std::to_string(counter));
        //Write MPO properties as attributes
        auto values = mpo->get_parameter_values();
        auto names  = mpo->get_parameter_names();
        for (size_t i = 0; i < std::min(values.size(), names.size()); i++){
            h5ppFile.writeAttributeToLink(values[i], names[i],sim_name + "/state/full/mpo/H_" + std::to_string(counter));
        }
        counter++;
    }
    h5ppFile.writeDataset(state.get_position() ,sim_name + "/state/full/position");
    h5ppFile.writeDataset(state.get_length()   ,sim_name + "/state/full/sites");
    h5ppFile.writeDataset(settings::model::model_type,sim_name + "/model/model_type");
    state.mpo_have_been_written_to_hdf5 = true;
}

void MPS_Tools::Finite::H5pp::write_hamiltonian_params(const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name){
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
    h5ppFile.writeDataset(hamiltonian_props,sim_name + "/model/full/Hamiltonian");
    int col = 0;
    for (auto &name : state.get_MPO_L().front()->get_parameter_names()){
        std::string attr_value = name;
        std::string attr_name  = "FIELD_" + std::to_string(col) + "_NAME";
        h5ppFile.writeAttributeToLink(attr_value, attr_name,sim_name + "/model/full/Hamiltonian");
        col++;
    }
    h5ppFile.writeDataset(settings::model::model_type,sim_name + "/model/model_type");

}

void MPS_Tools::Finite::H5pp::write_all_measurements(class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name){
    h5ppFile.writeDataset(state.measurements.length                      , sim_name + "/measurements/full/length");
    h5ppFile.writeDataset(state.measurements.norm                        , sim_name + "/measurements/full/norm");
    h5ppFile.writeDataset(state.measurements.bond_dimensions             , sim_name + "/measurements/full/bond_dimensions");
    h5ppFile.writeDataset(state.measurements.energy_per_site_mpo         , sim_name + "/measurements/full/energy_per_site_mpo");
    h5ppFile.writeDataset(state.measurements.energy_variance_per_site_mpo, sim_name + "/measurements/full/energy_variance_per_site_mpo");
    h5ppFile.writeDataset(state.measurements.entanglement_entropies      , sim_name + "/measurements/full/entanglement_entropies");
    h5ppFile.writeDataset(state.measurements.spin_components             , sim_name + "/measurements/full/spin_components");
    h5ppFile.writeDataset(state.measurements.spin_component_sx           , sim_name + "/measurements/full/spin_component_sx");
    h5ppFile.writeDataset(state.measurements.spin_component_sy           , sim_name + "/measurements/full/spin_component_sy");
    h5ppFile.writeDataset(state.measurements.spin_component_sz           , sim_name + "/measurements/full/spin_component_sz");
    bool OK = true;
    h5ppFile.writeDataset(OK, sim_name + "/simOK");
}


void MPS_Tools::Finite::H5pp::write_all_parity_projections(const class_finite_chain_state & state, const class_superblock &superblock, h5pp::File & h5ppFile, std::string sim_name){
    double var_sxup = write_parity_projected_analysis(state,superblock,h5ppFile,sim_name, "projections/sx_up", qm::spinOneHalf::sx, 1);
    double var_sxdn = write_parity_projected_analysis(state,superblock,h5ppFile,sim_name, "projections/sx_dn", qm::spinOneHalf::sx, -1);

//    write_parity_projected_analysis(state,superblock,hdf5,sim_name, "projections/sy_up", qm::spinOneHalf::sy, 1);
//    write_parity_projected_analysis(state,superblock,hdf5,sim_name, "projections/sy_dn", qm::spinOneHalf::sy, -1);
//    write_parity_projected_analysis(state,superblock,hdf5,sim_name, "projections/sz_up", qm::spinOneHalf::sz, 1);
//    write_parity_projected_analysis(state,superblock,hdf5,sim_name, "projections/sz_dn", qm::spinOneHalf::sz, -1);
    std::string best = var_sxup < var_sxdn ? "sx_up" : "sx_dn";
//    hdf5.write_symbolic_link("/" +sim_name + "/projections/" + best, "/" + sim_name + "/projections/best");
}



double MPS_Tools::Finite::H5pp::write_parity_projected_analysis(const class_finite_chain_state & state, const class_superblock &superblock, h5pp::File & h5ppFile, std::string sim_name, std::string projection_name,const Eigen::MatrixXcd paulimatrix, const int sign){
    if (std::abs(sign) != 1) throw std::runtime_error("Expected 'sign' +1 or -1. Got: " + std::to_string(sign));
    auto state_projected = MPS_Tools::Finite::Ops::get_parity_projected_state(state,paulimatrix, sign);
    state_projected.set_measured_false();
    state_projected.do_all_measurements();
    MPS_Tools::Finite::H5pp::write_all_state(state_projected,h5ppFile, sim_name + "/" + projection_name);
    MPS_Tools::Finite::H5pp::write_all_measurements(state_projected,h5ppFile, sim_name + "/" + projection_name);

    auto superblock_projected = superblock;
    MPS_Tools::Finite::Ops::rebuild_superblock(state_projected,superblock_projected);
    superblock_projected.set_measured_false();
    superblock_projected.do_all_measurements();
    MPS_Tools::Infinite::H5pp::write_all_superblock(superblock_projected,h5ppFile, sim_name + "/" + projection_name);
    MPS_Tools::Infinite::H5pp::write_all_measurements(superblock_projected,h5ppFile, sim_name + "/" + projection_name);

    return superblock.measurements.energy_variance_mpo;

}



void MPS_Tools::Finite::H5pp::load_from_hdf5(class_finite_chain_state & state, class_superblock & superblock, class_simulation_state &sim_state, h5pp::File & h5ppFile, std::string sim_name){
    // Load into state
    try{
        load_sim_state_from_hdf5(sim_state,h5ppFile,sim_name);
        load_state_from_hdf5(state,h5ppFile,sim_name);
        state.set_sweeps(sim_state.iteration);
        MPS_Tools::Finite::Ops::rebuild_superblock(state,superblock);
        MPS_Tools::Finite::Debug::check_integrity(state,superblock,sim_state);
    }catch(std::exception &ex){
        throw std::runtime_error("Failed to load from hdf5: " + std::string(ex.what()));
    }
}

void MPS_Tools::Finite::H5pp::load_state_from_hdf5(class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name){
    int    position = 0;
    size_t sites   = 0;
    Eigen::Tensor<Scalar,3> G;
    Eigen::Tensor<Scalar,1> L;
    Eigen::Tensor<Scalar,4> H;
    Eigen::MatrixXd Hamiltonian_params;
    std::string model_type;
    try{
        h5ppFile.readDataset(position             , sim_name + "/state/full/position");
        h5ppFile.readDataset(sites                , sim_name + "/state/full/sites");
        h5ppFile.readDataset(Hamiltonian_params   , sim_name + "/model/full/Hamiltonian");
        h5ppFile.readDataset(model_type           , sim_name + "/model/model_type");
    }catch (std::exception &ex){
        throw std::runtime_error("Couldn't read necessary model parameters: " + std::string(ex.what()));
    }
    state.clear();
    state.set_max_sites(sites);

    try {
        for(size_t i = 0; i < sites; i++){
            h5ppFile.readDataset(G, sim_name + "/state/full/mps/G_" + std::to_string(i));
            h5ppFile.readDataset(L, sim_name + "/state/full/mps/L_" + std::to_string(i));
            h5ppFile.readDataset(H, sim_name + "/state/full/mpo/H_" + std::to_string(i));
            if(i <= position ) {
                if(not state.get_MPS_L().empty() and state.get_MPS_L().back().get_chiR() != G.dimension(1)){
                    throw std::runtime_error("Mismatch in adjacent MPS dimensions");
                }
                state.get_MPS_L().emplace_back(G,L,i);
                state.get_MPO_L().emplace_back(class_hamiltonian_factory::create_mpo(model_type,Hamiltonian_params.row(i)));
            }
            else{
                if(not state.get_MPS_R().empty() and state.get_MPS_R().back().get_chiR() != G.dimension(1)){
                    throw std::runtime_error("Mismatch in adjacent MPS dimensions");
                }
                state.get_MPS_R().emplace_back(G,L,i);
                state.get_MPO_R().emplace_back(class_hamiltonian_factory::create_mpo(model_type,Hamiltonian_params.row(i)));
            }
        }
        h5ppFile.readDataset(state.get_MPS_C()    , sim_name + "/state/full/mps/L_C");
        if (state.get_MPS_L().size() + state.get_MPS_R().size() != (size_t)sites){
            throw std::runtime_error("Number of sites loaded does not match the number of sites advertised by the hdf5 file");
        }
        if (position != state.get_position()){
            throw std::runtime_error("Position loaded does not match the position read from the hdf5 file");
        }

    }catch (std::exception &ex){
        throw std::runtime_error("Could not read MPS/MPO tensors from file: " + std::string(ex.what()));
    }
    MPS_Tools::Finite::Ops::rebuild_environments(state);
}

void MPS_Tools::Finite::H5pp::load_sim_state_from_hdf5 (class_simulation_state &sim_state, h5pp::File & h5ppFile, std::string sim_name){
    sim_state.clear();
    // Common variables
    try{
        h5ppFile.readDataset(sim_state.iteration                      , sim_name + "/sim_state/iteration");
        h5ppFile.readDataset(sim_state.step                           , sim_name + "/sim_state/step");
        h5ppFile.readDataset(sim_state.position                       , sim_name + "/sim_state/position");
        h5ppFile.readDataset(sim_state.chi_temp                       , sim_name + "/sim_state/chi_temp");
        h5ppFile.readDataset(sim_state.min_sweeps                     , sim_name + "/sim_state/min_sweeps");
        h5ppFile.readDataset(sim_state.energy_min                     , sim_name + "/sim_state/energy_min");
        h5ppFile.readDataset(sim_state.energy_max                     , sim_name + "/sim_state/energy_max");
        h5ppFile.readDataset(sim_state.energy_target                  , sim_name + "/sim_state/energy_target");
        h5ppFile.readDataset(sim_state.energy_ubound                  , sim_name + "/sim_state/energy_ubound");
        h5ppFile.readDataset(sim_state.energy_lbound                  , sim_name + "/sim_state/energy_lbound");
        h5ppFile.readDataset(sim_state.energy_now                     , sim_name + "/sim_state/energy_now");
        h5ppFile.readDataset(sim_state.phys_time                      , sim_name + "/sim_state/phys_time");
        h5ppFile.readDataset(sim_state.delta_t                        , sim_name + "/sim_state/delta_t");
        h5ppFile.readDataset(sim_state.time_step_has_converged        , sim_name + "/sim_state/time_step_has_converged");
        h5ppFile.readDataset(sim_state.simulation_has_converged       , sim_name + "/sim_state/simulation_has_converged");
        h5ppFile.readDataset(sim_state.simulation_has_to_stop         , sim_name + "/sim_state/simulation_has_to_stop");
        h5ppFile.readDataset(sim_state.bond_dimension_has_reached_max , sim_name + "/sim_state/bond_dimension_has_reached_max");
        h5ppFile.readDataset(sim_state.entanglement_has_converged     , sim_name + "/sim_state/entanglement_has_converged");
        h5ppFile.readDataset(sim_state.entanglement_has_saturated     , sim_name + "/sim_state/entanglement_has_saturated");
        h5ppFile.readDataset(sim_state.variance_mpo_has_converged     , sim_name + "/sim_state/variance_mpo_has_converged");
        h5ppFile.readDataset(sim_state.variance_mpo_has_saturated     , sim_name + "/sim_state/variance_mpo_has_saturated");
        h5ppFile.readDataset(sim_state.variance_ham_has_converged     , sim_name + "/sim_state/variance_ham_has_converged");
        h5ppFile.readDataset(sim_state.variance_ham_has_saturated     , sim_name + "/sim_state/variance_ham_has_saturated");
        h5ppFile.readDataset(sim_state.variance_mom_has_converged     , sim_name + "/sim_state/variance_mom_has_converged");
        h5ppFile.readDataset(sim_state.variance_mom_has_saturated     , sim_name + "/sim_state/variance_mom_has_saturated");
        h5ppFile.readDataset(sim_state.variance_mpo_saturated_for     , sim_name + "/sim_state/variance_mpo_saturated_for");
        h5ppFile.readDataset(sim_state.variance_ham_saturated_for     , sim_name + "/sim_state/variance_ham_saturated_for");
        h5ppFile.readDataset(sim_state.variance_mom_saturated_for     , sim_name + "/sim_state/variance_mom_saturated_for");
    }catch(std::exception &ex){
        throw std::runtime_error("Failed to load sim_state from hdf5: " + std::string(ex.what()));
    }
}
