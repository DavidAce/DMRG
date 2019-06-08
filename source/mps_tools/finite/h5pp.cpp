//
// Created by david on 2019-03-09.
//


#include <mps_state/nmspc_mps_tools.h>
#include <mps_state/class_finite_chain_state.h>
#include <mps_state/class_superblock.h>
#include <model/class_hamiltonian_factory.h>
#include <algorithms/class_simulation_state.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <spdlog/spdlog.h>
#include <stdexcept>
#include <h5pp/h5pp.h>

using Scalar    = std::complex<double>;



void MPS_Tools::Finite::H5pp::write_all_state(const class_finite_chain_state &state, h5pp::File & h5ppFile,std::string sim_name) {
    switch(settings::hdf5::storage_level){
        case StorageLevel::NONE:
            break;
        case StorageLevel::LIGHT:
            write_bond_matrix(state,h5ppFile,sim_name);
            write_hamiltonian_params(state,h5ppFile,sim_name);
            break;
        case StorageLevel::NORMAL:
            write_bond_matrices(state,h5ppFile,sim_name);
            write_hamiltonian_params(state,h5ppFile,sim_name);
            break;
        case StorageLevel::FULL:
            write_full_mps(state,h5ppFile,sim_name);
            write_full_mpo(state,h5ppFile,sim_name);
            write_hamiltonian_params(state,h5ppFile,sim_name);
            break;
    }

    h5ppFile.writeDataset(state.get_position()       ,sim_name + "/state/position");
    h5ppFile.writeDataset(state.get_length()         ,sim_name + "/state/sites");
    h5ppFile.writeDataset(settings::model::model_type,sim_name + "/model/model_type");
}


void MPS_Tools::Finite::H5pp::write_bond_matrices(const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name)
/*! Writes down all the "Lambda" bond matrices (singular value matrices), so we can obtain the entanglement spectrum easily.
 */
{
    auto middle = (size_t) (state.get_length() / 2);
    for (size_t i = 0; i < state.get_length(); i++){
        if (i == middle){
            h5ppFile.writeDataset(state.get_L(i),sim_name + "/state/mps/L_C");
        }else{
            h5ppFile.writeDataset(state.get_L(i),sim_name + "/state/mps/L_" + std::to_string(i));
        }
    }
}


void MPS_Tools::Finite::H5pp::write_bond_matrix(const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name)
/*! Writes down all the "Lambda" bond matrices (singular value matrices), so we can obtain the entanglement spectrum easily.
 */
{
    auto middle = (size_t) (state.get_length() / 2);
    h5ppFile.writeDataset(state.get_L(middle),sim_name + "/state/mps/L_C");
}

void MPS_Tools::Finite::H5pp::write_full_mps(const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name)
/*!
 * Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation.
 *
 */
{
    write_bond_matrices(state,h5ppFile,sim_name);
    for (size_t i = 0; i < state.get_length(); i++){
        h5ppFile.writeDataset(state.get_G(i) ,sim_name + "/state/mps/G_" + std::to_string(i));
    }
}




void MPS_Tools::Finite::H5pp::write_full_mpo(const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name) {
    // Write all the MPO's
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    for (auto site = 0ul; site < state.get_length(); site++){
        h5ppFile.writeDataset(state.get_MPO(site).MPO(), sim_name + "/state/mpo/H_" + std::to_string(site));
        //Write MPO properties as attributes
        auto values = state.get_MPO(site).get_parameter_values();
        auto names  = state.get_MPO(site).get_parameter_names();
        for (size_t i = 0; i < std::min(values.size(), names.size()); i++){
            h5ppFile.writeAttributeToLink(values[i], names[i],sim_name + "/state/mpo/H_" + std::to_string(site));
        }
    }
}

void MPS_Tools::Finite::H5pp::write_hamiltonian_params(const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name){
    // Write down the Hamiltonian metadata as a table
    // Remember to write tensors in row-major state order because that's what hdf5 uses.
    Eigen::MatrixXd hamiltonian_props;
    for (auto site = 0ul ; site < state.get_length(); site++){
        auto props = state.get_MPO(site).get_parameter_values();
        Eigen::ArrayXd  temp_row  = Eigen::Map<Eigen::ArrayXd> (props.data(),props.size());
        hamiltonian_props.conservativeResize(hamiltonian_props.rows()+1, temp_row.size());
        hamiltonian_props.bottomRows(1) = temp_row.transpose();
    }
    h5ppFile.writeDataset(hamiltonian_props,sim_name + "/model/Hamiltonian");

    int col = 0;
    for (auto &name : state.MPO_L.front()->get_parameter_names()){
        std::string attr_value = name;
        std::string attr_name  = "FIELD_" + std::to_string(col) + "_NAME";
        h5ppFile.writeAttributeToLink(attr_value, attr_name,sim_name + "/model/Hamiltonian");
        col++;
    }
}

void MPS_Tools::Finite::H5pp::write_all_measurements(const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name){
    state.do_all_measurements();
    h5ppFile.writeDataset(state.measurements.length.value()                      , sim_name + "/measurements/length");
    h5ppFile.writeDataset(state.measurements.norm.value()                        , sim_name + "/measurements/norm");
    h5ppFile.writeDataset(state.measurements.bond_dimensions.value()             , sim_name + "/measurements/bond_dimensions");
    h5ppFile.writeDataset(state.measurements.energy_per_site_mpo.value()         , sim_name + "/measurements/energy_per_site_mpo");
    h5ppFile.writeDataset(state.measurements.energy_variance_per_site_mpo.value(), sim_name + "/measurements/energy_variance_per_site_mpo");
    h5ppFile.writeDataset(state.measurements.entanglement_entropies.value()      , sim_name + "/measurements/entanglement_entropies");
    h5ppFile.writeDataset(state.measurements.spin_components.value()             , sim_name + "/measurements/spin_components");
    h5ppFile.writeDataset(state.measurements.spin_component_sx.value()           , sim_name + "/measurements/spin_component_sx");
    h5ppFile.writeDataset(state.measurements.spin_component_sy.value()           , sim_name + "/measurements/spin_component_sy");
    h5ppFile.writeDataset(state.measurements.spin_component_sz.value()           , sim_name + "/measurements/spin_component_sz");

}


void MPS_Tools::Finite::H5pp::write_closest_parity_projection(const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name, std::string paulistring){
    auto state_projected = MPS_Tools::Finite::Ops::get_closest_parity_state(state,paulistring);
    state_projected.unset_measurements();
    state_projected.do_all_measurements();
    MPS_Tools::Finite::H5pp::write_all_state(state_projected,h5ppFile, sim_name + "/projections/" + paulistring);
    MPS_Tools::Finite::H5pp::write_all_measurements(state_projected,h5ppFile, sim_name + "/projections/" + paulistring);
}


void MPS_Tools::Finite::H5pp::load_from_hdf5(const h5pp::File & h5ppFile, class_finite_chain_state & state, class_superblock & superblock, class_simulation_state &sim_state, std::string sim_name){
    // Load into state
    try{
        load_sim_state_from_hdf5(h5ppFile,sim_state,sim_name);
        load_state_from_hdf5(h5ppFile,state,sim_name);
        state.set_sweeps(sim_state.iteration);
        MPS_Tools::Finite::Ops::rebuild_superblock(state,superblock);
        MPS_Tools::Finite::Debug::check_integrity(state,superblock,sim_state);
    }catch(std::exception &ex){
        throw std::runtime_error("Failed to load from hdf5: " + std::string(ex.what()));
    }
}

void MPS_Tools::Finite::H5pp::load_state_from_hdf5(const h5pp::File & h5ppFile, class_finite_chain_state & state, std::string sim_name){
    int    position = 0;
    size_t sites   = 0;
    Eigen::Tensor<Scalar,3> G;
    Eigen::Tensor<Scalar,1> L;
    Eigen::Tensor<Scalar,4> H;
    Eigen::MatrixXd Hamiltonian_params;
    std::string model_type;
    try{
        h5ppFile.readDataset(position             , sim_name + "/state/position");
        h5ppFile.readDataset(sites                , sim_name + "/state/sites");
        h5ppFile.readDataset(Hamiltonian_params   , sim_name + "/model/Hamiltonian");
        h5ppFile.readDataset(model_type           , sim_name + "/model/model_type");
    }catch (std::exception &ex){
        throw std::runtime_error("Couldn't read necessary model parameters: " + std::string(ex.what()));
    }
    state.clear();
    state.set_max_sites(sites);

    try {
        for(size_t i = 0; i < sites; i++){
            h5ppFile.readDataset(G, sim_name + "/state/mps/G_" + std::to_string(i));
            h5ppFile.readDataset(L, sim_name + "/state/mps/L_" + std::to_string(i));
            h5ppFile.readDataset(H, sim_name + "/state/mpo/H_" + std::to_string(i));
            if(i <= (size_t)position ) {
                if(not state.MPS_L.empty() and state.MPS_L.back().get_chiR() != G.dimension(1)){
                    throw std::runtime_error("Mismatch in adjacent MPS dimensions");
                }
                state.MPS_L.emplace_back(G,L,i);
                state.MPO_L.emplace_back(class_hamiltonian_factory::create_mpo(model_type,Hamiltonian_params.row(i)));
            }
            else{
                if(not state.MPS_R.empty() and state.MPS_R.back().get_chiR() != G.dimension(1)){
                    throw std::runtime_error("Mismatch in adjacent MPS dimensions");
                }
                state.MPS_R.emplace_back(G,L,i);
                state.MPO_R.emplace_back(class_hamiltonian_factory::create_mpo(model_type,Hamiltonian_params.row(i)));
            }
        }
        h5ppFile.readDataset(state.MPS_C    , sim_name + "/state/mps/L_C");
        if (state.MPS_L.size() + state.MPS_R.size() != (size_t)sites){
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

void MPS_Tools::Finite::H5pp::load_sim_state_from_hdf5 (const h5pp::File & h5ppFile, class_simulation_state &sim_state, std::string sim_name){
    sim_state.clear();
    // Common variables
    try{
        h5ppFile.readDataset(sim_state.iteration                      , sim_name + "/sim_state/iteration");
        h5ppFile.readDataset(sim_state.step                           , sim_name + "/sim_state/step");
        h5ppFile.readDataset(sim_state.position                       , sim_name + "/sim_state/position");
        h5ppFile.readDataset(sim_state.chi_temp                       , sim_name + "/sim_state/chi_temp");
        h5ppFile.readDataset(sim_state.chi_max                        , sim_name + "/sim_state/chi_temp");
        h5ppFile.readDataset(sim_state.min_sweeps                     , sim_name + "/sim_state/min_sweeps");
        h5ppFile.readDataset(sim_state.energy_min                     , sim_name + "/sim_state/energy_min");
        h5ppFile.readDataset(sim_state.energy_max                     , sim_name + "/sim_state/energy_max");
        h5ppFile.readDataset(sim_state.energy_target                  , sim_name + "/sim_state/energy_target");
        h5ppFile.readDataset(sim_state.energy_ubound                  , sim_name + "/sim_state/energy_ubound");
        h5ppFile.readDataset(sim_state.energy_lbound                  , sim_name + "/sim_state/energy_lbound");
        h5ppFile.readDataset(sim_state.energy_now                     , sim_name + "/sim_state/energy_now");
        h5ppFile.readDataset(sim_state.energy_dens                    , sim_name + "/sim_state/energy_dens");
        h5ppFile.readDataset(sim_state.energy_dens_target             , sim_name + "/sim_state/energy_dens_target");
        h5ppFile.readDataset(sim_state.energy_dens_window             , sim_name + "/sim_state/energy_dens_window");
        h5ppFile.readDataset(sim_state.phys_time                      , sim_name + "/sim_state/phys_time");
        h5ppFile.readDataset(sim_state.wall_time                      , sim_name + "/sim_state/wall_time");
        h5ppFile.readDataset(sim_state.simu_time                      , sim_name + "/sim_state/simu_time");
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
