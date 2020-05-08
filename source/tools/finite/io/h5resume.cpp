//
// Created by david on 2019-11-07.
//
#include <complex>
#include <h5pp/h5pp.h>
#include <model/class_model_factory.h>
#include <simulation/class_simulation_status.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_finite.h>
#include <tools/common/io.h>
#include <tools/finite/debug.h>
#include <tools/finite/io.h>
#include <tools/finite/mpo.h>
#include <tools/finite/mps.h>
#include <typeindex>

using Scalar = std::complex<double>;

// Load model, state and simulation status from HDF5
void tools::finite::io::h5resume::load_all(const h5pp::File &h5ppFile, const std::string &prefix, class_simulation_status &sim_status,
                                           class_state_finite &state) {
    try {
        if(h5ppFile.readAttribute<std::string>("storage_level", prefix) != enum2str(StorageLevel::FULL))
            throw std::runtime_error("Given prefix to simulation data with StorageLevel < FULL. The simulation can only be resumed from FULL storage");

        tools::finite::io::h5resume::load_mpo(h5ppFile, prefix, state);
        tools::finite::io::h5resume::load_mps(h5ppFile, prefix, state);
        tools::common::io::h5resume::load_sim_status_from_hdf5(h5ppFile, prefix, sim_status);
        tools::common::io::h5resume::load_profiling_from_hdf5(h5ppFile, prefix);

        tools::finite::mps::rebuild_environments(state);
        tools::finite::mps::normalize(state);

        tools::finite::debug::check_integrity(state);
    } catch(std::exception &ex) {
        throw std::runtime_error("Failed to load simulation from hdf5 file: " + std::string(ex.what()));
    }
}

void tools::finite::io::h5resume::load_mpo(const h5pp::File &h5ppFile, const std::string &prefix, class_state_finite &state) {
    if(h5ppFile.readAttribute<std::string>("storage_level", prefix) != enum2str(StorageLevel::FULL))
        throw std::runtime_error("Given prefix to MPO data with StorageLevel < FULL. The MPO's can only be resumed from FULL storage");
    settings::model::model_type = h5ppFile.readAttribute<std::string>("model_type", prefix);
    auto model_path             = h5ppFile.readAttribute<std::string>("model_path", prefix);
    auto num_sites              = h5ppFile.readAttribute<size_t>("sites", prefix+"/mpo");
    auto position               = h5ppFile.readAttribute<size_t>("position", prefix+"/mpo");
    tools::finite::mpo::initialize(state, settings::model::model_type, num_sites, position);
    for(size_t pos = 0; pos < num_sites; pos++) state.get_MPO(pos).read_parameters(h5ppFile, model_path);
}

void tools::finite::io::h5resume::load_mps(const h5pp::File &h5ppFile, const std::string &prefix, class_state_finite &state) {
    if(h5ppFile.readAttribute<std::string>("storage_level", prefix) != enum2str(StorageLevel::FULL))
        throw std::runtime_error("Given prefix to MPS data with StorageLevel < FULL. The MPS's can only be resumed from FULL storage");
    settings::model::model_type = h5ppFile.readAttribute<std::string>("model_type", prefix);
    auto model_path             = h5ppFile.readAttribute<std::string>("model_path", prefix);
    auto num_sites              = h5ppFile.readAttribute<size_t>("sites", prefix+"/mps");
    auto position               = h5ppFile.readAttribute<size_t>("position", prefix+"/mps");
    tools::finite::mps::initialize(state, settings::model::model_type, num_sites, position);
    for(size_t pos = 0; pos < num_sites; pos++) {
        if( state.get_MPS(pos).isCenter()) {
            if(not h5ppFile.linkExists(prefix + "/mps/L_C")) throw std::runtime_error("Dataset does not exist: " + prefix + "/mps/L_C");
            auto LC = h5ppFile.readDataset<Eigen::Tensor<Scalar, 1>>(prefix + "/mps/L_C");
            state.get_MPS(pos).set_LC(LC);
        }
        if(not h5ppFile.linkExists(prefix + "/mps/L_" + std::to_string(pos))) throw std::runtime_error("Dataset does not exist: " + prefix + "/mps/L_" + std::to_string(pos));
        if(not h5ppFile.linkExists(prefix + "/mps/M_" + std::to_string(pos))) throw std::runtime_error("Dataset does not exist: " + prefix + "/mps/M_" + std::to_string(pos));
        auto  L   = h5ppFile.readDataset<Eigen::Tensor<Scalar, 1>>(prefix + "/mps/L_" + std::to_string(pos));
        auto  M   = h5ppFile.readDataset<Eigen::Tensor<Scalar, 3>>(prefix + "/mps/M_" + std::to_string(pos));
        state.get_MPS(pos).set_mps(M, L);
        //Sanity checks
        if(pos == position and not  state.get_MPS(pos).isCenter())
            throw std::logic_error("Given position is not a a center. State may be wrongly initialized or something is wrong with the resumed file");
        if(pos != position and  state.get_MPS(pos).isCenter())
            throw std::logic_error("A site not at current position claims to be a state center");
//        if(passed_LC > 1) throw std::logic_error("Multiple centers encountered");

    }
    if(state.MPS_L.size() + state.MPS_R.size() != num_sites) throw std::logic_error("Initialized MPS with the wrong number of sites");
    if(not state.get_MPS(position).isCenter()) throw std::logic_error("Initialized MPS with the center at the wrong position");
    state.set_iter(h5ppFile.readAttribute<size_t>("iteration", prefix + "/mps"));
    state.set_step(h5ppFile.readAttribute<size_t>("step", prefix + "/mps"));
    state.set_chi_lim(h5ppFile.readAttribute<long>("chi_lim", prefix + "/mps"));
    state.set_chi_max(h5ppFile.readAttribute<long>("chi_max", prefix + "/mps"));
}
