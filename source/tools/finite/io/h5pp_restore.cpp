//
// Created by david on 2019-11-07.
//
#include <complex>
#include <tools/finite/mps.h>
#include <tools/finite/io.h>
#include <tools/finite/debug.h>
#include <tools/common/io.h>
#include <h5pp/h5pp.h>
#include <model/class_model_factory.h>
#include <simulation/class_simulation_status.h>
#include <state/class_state_finite.h>
#include <typeindex>

using Scalar = std::complex<double>;

void tools::finite::io::h5restore::load_from_hdf5(const h5pp::File &h5ppFile, class_state_finite &state, class_simulation_status &sim_status, const std::string &prefix_path) {
    // Load into state
    try {
        sim_status = tools::common::io::h5restore::load_sim_status_from_hdf5(h5ppFile, prefix_path);
        state      = tools::finite::io::h5restore::load_state_from_hdf5(h5ppFile, prefix_path);
        state.set_sweeps(sim_status.iteration);
        tools::finite::debug::check_integrity(state);
    } catch(std::exception &ex) { throw std::runtime_error("Failed to load from hdf5 file: " + std::string(ex.what())); }
}


class_model_base::Parameters get_parameters(const h5pp::File &h5ppFile, const std::string & dsetName){
    class_model_base::Parameters parameters;
    auto attrTypeInfo = h5ppFile.getAttributeTypeInfoAll(dsetName);
    for (auto & attr : attrTypeInfo){
        if(attr.type() == typeid(int) and attr.size() == 1)
            parameters.push_back({attr.name(), h5ppFile.readAttribute<int>(attr.name(),dsetName)});
        if(attr.type() == typeid(double) and attr.size() == 1)
            parameters.push_back({attr.name(), h5ppFile.readAttribute<double>(attr.name(),dsetName)});
        if(attr.type() == typeid(size_t) and attr.size() == 1)
            parameters.push_back({attr.name(), h5ppFile.readAttribute<size_t>(attr.name(),dsetName)});
        if(attr.type() == typeid(std::string) and attr.size() == 1)
            parameters.push_back({attr.name(), h5ppFile.readAttribute<std::string>(attr.name(),dsetName)});
        if(attr.type() == typeid(bool) and attr.size() == 1)
            parameters.push_back({attr.name(), h5ppFile.readAttribute<bool>(attr.name(),dsetName)});
    }
    return parameters;
}


class_state_finite tools::finite::io::h5restore::load_state_from_hdf5(const h5pp::File &h5ppFile, const std::string &prefix_path) {
    class_state_finite       state;
    size_t                   position = 0;
    size_t                   sites    = 0;
    Eigen::Tensor<Scalar, 3> G;
    Eigen::Tensor<Scalar, 1> L;
    Eigen::Tensor<Scalar, 4> H;
    std::vector<class_model_base::Parameters> hamiltonian_parameters;
    std::string              model_type;
    try {
        h5ppFile.readDataset(position, prefix_path + "/state/position");
        h5ppFile.readDataset(sites, prefix_path + "/state/sites");
        h5ppFile.readDataset(model_type, prefix_path + "/model/model_type");
        for(size_t site = 0; site < sites; site++){
            hamiltonian_parameters.push_back(get_parameters(h5ppFile,"/model/hamiltonian_" + std::to_string(site)));
        }

    } catch(std::exception &ex) { throw std::runtime_error("Couldn't read necessary model parameters: " + std::string(ex.what())); }

    try {
        for(size_t i = 0; i < sites; i++) {
            h5ppFile.readDataset(G, prefix_path + "/state/mps/G_" + std::to_string(i));
            h5ppFile.readDataset(L, prefix_path + "/state/mps/L_" + std::to_string(i));
            h5ppFile.readDataset(H, prefix_path + "/state/mpo/H_" + std::to_string(i));
            if(i <= (size_t) position) {
                if(not state.MPS_L.empty() and state.MPS_L.back().get_chiR() != G.dimension(1)) { throw std::runtime_error("Mismatch in adjacent MPS dimensions"); }
                state.MPS_L.emplace_back(G, L, i);
            } else {
                if(not state.MPS_R.empty() and state.MPS_R.back().get_chiR() != G.dimension(1)) { throw std::runtime_error("Mismatch in adjacent MPS dimensions"); }
                state.MPS_R.emplace_back(G, L, i);
            }
        }
        state.set_positions();
        size_t center_pos = (state.get_length() - 1) / 2;
        h5ppFile.readDataset(state.get_MPS(center_pos).get_LC(), prefix_path + "/state/mps/L_C");
        if(state.MPS_L.size() + state.MPS_R.size() != (size_t) sites) {
            throw std::runtime_error("Number of sites loaded does not match the number of sites advertised by the output file");
        }
        if(position != state.get_position()) { throw std::runtime_error("Position loaded does not match the position read from the output file"); }

    } catch(std::exception &ex) { throw std::runtime_error("Could not read MPS/MPO tensors from file: " + std::string(ex.what())); }
    tools::finite::mps::rebuild_environments(state);
    return state;
}