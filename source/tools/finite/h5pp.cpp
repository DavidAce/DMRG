//
// Created by david on 2019-03-09.
//


#include <tools/nmspc_tools.h>
#include <state/class_finite_state.h>
#include <state/class_infinite_state.h>
#include <model/class_model_factory.h>
#include <simulation/class_simulation_status.h>
#include <simulation/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <stdexcept>
#include <h5pp/h5pp.h>

using Scalar    = std::complex<double>;


bool tools::finite::io::internals::make_extendable_dataset(const std::string & prefix_path) {
    std::string log_string = "/log";  // A logged dataset is not supposed to be extended, just written once.
    return prefix_path.find(log_string) == std::string::npos;
}




void tools::finite::io::write_all_state(const class_finite_state &state, h5pp::File & h5ppFile,const std::string & prefix_path) {
    tools::common::profile::t_hdf.tic();
    switch(settings::output::storage_level){
        case StorageLevel::NONE:
            break;
        case StorageLevel::LIGHT:
            write_bond_matrix(state,h5ppFile,prefix_path);
            break;
        case StorageLevel::NORMAL:
            write_bond_matrices(state,h5ppFile,prefix_path);
            break;
        case StorageLevel::FULL:
            write_full_mps(state,h5ppFile,prefix_path);
            write_full_mpo(state,h5ppFile,prefix_path);
            break;
    }
    h5ppFile.writeDataset(state.get_position()       ,prefix_path + "/state/position");
    h5ppFile.writeDataset(state.get_length()         ,prefix_path + "/state/sites");
    h5ppFile.writeDataset(settings::model::model_type,prefix_path + "/model/model_type");
    tools::common::profile::t_hdf.toc();

}


void tools::finite::io::write_bond_matrices(const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path)
/*! Writes down all the "Lambda" bond matrices (singular value matrices), so we can obtain the entanglement spectrum easily.
 */
{
    bool extendable = internals::make_extendable_dataset(prefix_path);
    auto middle = (size_t) (state.get_length() / 2);
    for (size_t i = 0; i < state.get_length(); i++){
        h5ppFile.writeDataset(state.get_MPS(i).get_L(),prefix_path + "/state/mps/L_" + std::to_string(i),extendable);
        if (i == middle){
            h5ppFile.writeDataset(state.midchain_bond(), prefix_path + "/state/mps/L_C", extendable);
        }
    }
    h5ppFile.writeDataset(state.get_truncation_error(), prefix_path + "/state/truncation_error",extendable);
}


void tools::finite::io::write_bond_matrix(const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path)
/*! Writes down all the "Lambda" bond matrices (singular value matrices), so we can obtain the entanglement spectrum easily.
 */
{
    bool extendable = internals::make_extendable_dataset(prefix_path);
    auto middle = (size_t) (state.get_length() / 2);
    h5ppFile.writeDataset(state.get_MPS(middle).get_LC(),prefix_path + "/state/mps/L_C",extendable);
    h5ppFile.writeDataset(state.get_truncation_error(middle), prefix_path + "/state/truncation_error",extendable);
}

void tools::finite::io::write_full_mps(const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path)
/*!
 * Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation.
 *
 */
{
    bool extendable = internals::make_extendable_dataset(prefix_path);
    write_bond_matrices(state,h5ppFile,prefix_path);
    for (size_t i = 0; i < state.get_length(); i++){
        h5ppFile.writeDataset(state.get_MPS(i).get_M() ,prefix_path + "/state/mps/M_" + std::to_string(i),extendable);
    }
}




void tools::finite::io::write_full_mpo(const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path) {
    // Write all the MPO's
    // Remember to write tensors in row-major state order because that's what output uses.
    bool extendable = internals::make_extendable_dataset(prefix_path);
    for (auto site = 0ul; site < state.get_length(); site++){
        h5ppFile.writeDataset(state.get_MPO(site).MPO(), prefix_path + "/state/mpo/H_" + std::to_string(site),extendable);
        //Write MPO properties as attributes
        auto values = state.get_MPO(site).get_parameter_values();
        auto names  = state.get_MPO(site).get_parameter_names();
        for (size_t i = 0; i < std::min(values.size(), names.size()); i++){
            h5ppFile.writeAttributeToLink(values[i], names[i],prefix_path + "/state/mpo/H_" + std::to_string(site));
        }
    }
}

void tools::finite::io::write_model(const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path){
    // Write down the Hamiltonian metadata as a table
    // Remember to write tensors in row-major state order because that's what output uses.
    if(settings::output::storage_level == StorageLevel::NONE) return;
    Eigen::MatrixXd hamiltonian_props;
    for (auto site = 0ul ; site < state.get_length(); site++){
        auto props = state.get_MPO(site).get_parameter_values();
        Eigen::ArrayXd  temp_row  = Eigen::Map<Eigen::ArrayXd> (props.data(),props.size());
        hamiltonian_props.conservativeResize(hamiltonian_props.rows()+1, temp_row.size());
        hamiltonian_props.bottomRows(1) = temp_row.transpose();
    }
    bool extendable = internals::make_extendable_dataset(prefix_path);
    h5ppFile.writeDataset(hamiltonian_props,prefix_path + "/model/Hamiltonian",extendable);
    int col = 0;
    for (auto &name : state.MPO_L.front()->get_parameter_names()){
        std::string attr_value = name;
        std::string attr_name  = "FIELD_" + std::to_string(col) + "_NAME";
        h5ppFile.writeAttributeToLink(attr_value, attr_name,prefix_path + "/model/Hamiltonian");
        col++;
    }
}

void tools::finite::io::write_all_measurements(const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path){
    state.do_all_measurements();
    tools::common::profile::t_hdf.tic();
    h5ppFile.writeDataset(state.measurements.length.value()                        , prefix_path + "/measurements/length");
    h5ppFile.writeDataset(state.measurements.norm.value()                          , prefix_path + "/measurements/norm");
    h5ppFile.writeDataset(state.measurements.bond_dimensions.value()               , prefix_path + "/measurements/bond_dimensions");
    h5ppFile.writeDataset(state.measurements.bond_dimension_midchain.value()       , prefix_path + "/measurements/bond_dimension_midchain");
    h5ppFile.writeDataset(state.measurements.energy.value()                        , prefix_path + "/measurements/energy");
    h5ppFile.writeDataset(state.measurements.energy_per_site.value()               , prefix_path + "/measurements/energy_per_site");
    h5ppFile.writeDataset(state.measurements.energy_variance_mpo.value()           , prefix_path + "/measurements/energy_variance_mpo");
    h5ppFile.writeDataset(state.measurements.energy_variance_per_site.value()      , prefix_path + "/measurements/energy_variance_per_site");
    h5ppFile.writeDataset(state.measurements.entanglement_entropies.value()        , prefix_path + "/measurements/entanglement_entropies");
    h5ppFile.writeDataset(state.measurements.entanglement_entropy_midchain.value() , prefix_path + "/measurements/entanglement_entropy_midchain");
    h5ppFile.writeDataset(state.measurements.spin_components.value()               , prefix_path + "/measurements/spin_components");
    h5ppFile.writeDataset(state.measurements.spin_component_sx.value()             , prefix_path + "/measurements/spin_component_sx");
    h5ppFile.writeDataset(state.measurements.spin_component_sy.value()             , prefix_path + "/measurements/spin_component_sy");
    h5ppFile.writeDataset(state.measurements.spin_component_sz.value()             , prefix_path + "/measurements/spin_component_sz");
    tools::common::profile::t_hdf.toc();

}


void tools::finite::io::write_projection_to_closest_parity_sector(const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path, std::string parity_sector){
    if (parity_sector == "none") return;
    auto state_projected = tools::finite::ops::get_projection_to_closest_parity_sector(state, parity_sector);
    state_projected.unset_measurements();
    state_projected.do_all_measurements();
    tools::finite::io::write_all_state(state_projected,h5ppFile, prefix_path + "/projections/" + parity_sector);
    tools::finite::io::write_all_measurements(state_projected,h5ppFile, prefix_path + "/projections/" + parity_sector);
}


void tools::finite::io::load_from_hdf5(const h5pp::File & h5ppFile, class_finite_state & state, class_simulation_status &sim_status, const std::string & prefix_path){
    // Load into state
    try{
        sim_status = tools::common::io::load_sim_status_from_hdf5(h5ppFile,prefix_path);
        state     = load_state_from_hdf5(h5ppFile,prefix_path);
        state.set_sweeps(sim_status.iteration);
        tools::finite::debug::check_integrity(state);
    }catch(std::exception &ex){
        throw std::runtime_error("Failed to load from output: " + std::string(ex.what()));
    }
}

class_finite_state tools::finite::io::load_state_from_hdf5(const h5pp::File & h5ppFile, const std::string & prefix_path){
    class_finite_state state;
    size_t position = 0;
    size_t sites   = 0;
    Eigen::Tensor<Scalar,3> G;
    Eigen::Tensor<Scalar,1> L;
    Eigen::Tensor<Scalar,4> H;
    Eigen::MatrixXd Hamiltonian_params;
    std::string model_type;
    try{
        h5ppFile.readDataset(position             , prefix_path + "/state/position");
        h5ppFile.readDataset(sites                , prefix_path + "/state/sites");
        h5ppFile.readDataset(Hamiltonian_params   , prefix_path + "/model/Hamiltonian");
        h5ppFile.readDataset(model_type           , prefix_path + "/model/model_type");
    }catch (std::exception &ex){
        throw std::runtime_error("Couldn't read necessary model parameters: " + std::string(ex.what()));
    }

    try {
        for(size_t i = 0; i < sites; i++){
            h5ppFile.readDataset(G, prefix_path + "/state/mps/G_" + std::to_string(i));
            h5ppFile.readDataset(L, prefix_path + "/state/mps/L_" + std::to_string(i));
            h5ppFile.readDataset(H, prefix_path + "/state/mpo/H_" + std::to_string(i));
            if(i <= (size_t)position ) {
                if(not state.MPS_L.empty() and state.MPS_L.back().get_chiR() != G.dimension(1)){
                    throw std::runtime_error("Mismatch in adjacent MPS dimensions");
                }
                state.MPS_L.emplace_back(G,L,i);
                state.MPO_L.emplace_back(class_model_factory::create_mpo(i,model_type,Hamiltonian_params.row(i)));
            }
            else{
                if(not state.MPS_R.empty() and state.MPS_R.back().get_chiR() != G.dimension(1)){
                    throw std::runtime_error("Mismatch in adjacent MPS dimensions");
                }
                state.MPS_R.emplace_back(G,L,i);
                state.MPO_R.emplace_back(class_model_factory::create_mpo(i,model_type,Hamiltonian_params.row(i)));
            }
        }
        h5ppFile.readDataset(state.midchain_bond()    , prefix_path + "/state/mps/L_C");
        if (state.MPS_L.size() + state.MPS_R.size() != (size_t)sites){
            throw std::runtime_error("Number of sites loaded does not match the number of sites advertised by the output file");
        }
        if (position != state.get_position()){
            throw std::runtime_error("Position loaded does not match the position read from the output file");
        }

    }catch (std::exception &ex){
        throw std::runtime_error("Could not read MPS/MPO tensors from file: " + std::string(ex.what()));
    }
    tools::finite::mps::rebuild_environments(state);
    return state;
}
