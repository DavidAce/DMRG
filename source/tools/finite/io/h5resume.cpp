//
// Created by david on 2019-11-07.
//
#include <complex>
#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <model/class_model_factory.h>
#include <simulation/class_simulation_status.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_finite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/finite/debug.h>
#include <tools/finite/io.h>
#include <tools/finite/mpo.h>
#include <tools/finite/mps.h>
#include <typeindex>

using Scalar = std::complex<double>;

// Load model, state and simulation status from HDF5
void tools::finite::io::h5resume::load_all(const h5pp::File &h5ppFile, const std::string &state_prefix, class_state_finite &state,
                                           class_simulation_status &sim_status) {
    try {
        state.clear_measurements();
        state.clear_cache();

        if(h5ppFile.readAttribute<std::string>("storage_level", state_prefix) != enum2str(StorageLevel::FULL))
            throw std::runtime_error("Given prefix to simulation data with StorageLevel < FULL. The simulation can only be resumed from FULL storage");
        tools::common::io::h5resume::load_sim_status_from_hdf5(h5ppFile, state_prefix, sim_status);
        tools::finite::io::h5resume::load_mpo(h5ppFile, state_prefix, state, sim_status);
        tools::finite::io::h5resume::load_mps(h5ppFile, state_prefix, state, sim_status);
        tools::common::io::h5resume::load_profiling_from_hdf5(h5ppFile, state_prefix);
        tools::finite::mps::rebuild_environments(state);
        tools::finite::mps::normalize(state);
        tools::finite::io::h5resume::validate(h5ppFile, state_prefix, state, sim_status);
    } catch(std::exception &ex) {
        throw std::runtime_error("Failed to load simulation from hdf5 file: " + std::string(ex.what()));
    }
}

void tools::finite::io::h5resume::load_mpo(const h5pp::File &h5ppFile, const std::string &state_prefix, class_state_finite &state,
                                           const class_simulation_status &sim_status) {
    if(h5ppFile.readAttribute<std::string>("storage_level", state_prefix) != enum2str(StorageLevel::FULL))
        throw std::runtime_error("Given prefix to MPO data with StorageLevel < FULL. The MPO's can only be resumed from FULL storage");
    auto position = h5ppFile.readAttribute<size_t>("position", state_prefix);
    if(position != sim_status.position)
        throw std::runtime_error(fmt::format("Mismatch when loading MPO: State position [{}] != sim_status.position [{}]", position, sim_status.position));

    // Find the path to the MPO
    auto mpo_path = h5ppFile.readAttribute<std::string>(state_prefix, "common/mpo_path");
    auto ham_path = h5ppFile.readAttribute<std::string>(state_prefix, "common/ham_path");
    if(h5ppFile.linkExists(ham_path)) {
        tools::log->trace("Initializing MPOs from Hamiltonian table on file: [{}]", ham_path);
        auto model_type = h5ppFile.readAttribute<std::string>("model_type", ham_path);
        auto model_size = h5ppFile.readAttribute<size_t>("model_size", ham_path);
        if(model_type != settings::model::model_type)
            throw std::runtime_error(
                fmt::format("Mismatch when loading MPO: model_type [{}] != settings::model::model_type [{}]", model_type, settings::model::model_type));
        if(model_size != settings::model::model_size)
            throw std::runtime_error(
                fmt::format("Mismatch when loading MPO: model_size [{}] != settings::model::model_size [{}]", model_type, settings::model::model_size));

        tools::finite::mpo::initialize(state, model_type, model_size, position);
        for(size_t pos = 0; pos < model_size; pos++) state.get_MPO(pos).read_hamiltonian(h5ppFile, ham_path);

    } else if(h5ppFile.linkExists(mpo_path)) {
        tools::log->trace("Initializing MPOs from MPO's on file: [{}]", mpo_path);
        auto model_type = h5ppFile.readAttribute<std::string>("model_type", mpo_path);
        auto model_size = h5ppFile.readAttribute<size_t>("model_size", mpo_path);
        if(model_type != settings::model::model_type)
            throw std::runtime_error(fmt::format("Mismatch when loading MPO: Hamiltonian model_type [{}] != settings::model::model_type [{}]", model_type,
                                                 settings::model::model_type));
        if(model_size != settings::model::model_size)
            throw std::runtime_error(
                fmt::format("Mismatch when loading MPO: Hamiltonian sites [{}] != settings::model::model_size [{}]", model_type, settings::model::model_size));

        tools::finite::mpo::initialize(state, model_type, model_size, position);
        for(size_t pos = 0; pos < model_size; pos++) state.get_MPO(pos).read_mpo(h5ppFile, mpo_path);
    }
}

void tools::finite::io::h5resume::load_mps(const h5pp::File &h5ppFile, const std::string &state_prefix, class_state_finite &state,
                                           const class_simulation_status &sim_status) {
    if(h5ppFile.readAttribute<std::string>("storage_level", state_prefix) != enum2str(StorageLevel::FULL))
        throw std::runtime_error("Given prefix to MPS data with StorageLevel < FULL. The MPS's can only be resumed from FULL storage");
    std::string mps_path   = state_prefix + "/mps";
    auto        model_type = h5ppFile.readAttribute<std::string>("model_type", state_prefix);
    auto        model_size = h5ppFile.readAttribute<size_t>("model_size", mps_path);
    auto        position   = h5ppFile.readAttribute<size_t>("position", mps_path);
    if(position != sim_status.position)
        throw std::runtime_error(fmt::format("Mismatch when loading MPS: position [{}] != sim_status.position [{}]", position, sim_status.position));
    if(model_type != settings::model::model_type)
        throw std::runtime_error(
            fmt::format("Mismatch when loading MPS: model_type [{}] != settings::model::model_type [{}]", model_type, settings::model::model_type));
    if(model_size != settings::model::model_size)
        throw std::runtime_error(
            fmt::format("Mismatch when loading MPS: model_size [{}] != settings::model::model_size [{}]", model_size, settings::model::model_size));

    tools::finite::mps::initialize(state, model_type, model_size, position);
    for(size_t pos = 0; pos < model_size; pos++) {
        std::string pos_str = std::to_string(pos);
        if(state.get_MPS(pos).isCenter()) {
            if(not h5ppFile.linkExists(mps_path + "/L_C")) throw std::runtime_error("Dataset does not exist: " + mps_path + "/L_C");
            auto LC          = h5ppFile.readDataset<Eigen::Tensor<Scalar, 1>>(mps_path + "/L_C");
            auto pos_on_file = h5ppFile.readAttribute<size_t>("position", mps_path + "/L_C");
            if(pos != pos_on_file) throw std::runtime_error(fmt::format("Center bond position mismatch: pos [{}] != pos on file [{}]", pos, pos_on_file));
            state.get_MPS(pos).set_LC(LC);
        }
        if(not h5ppFile.linkExists(mps_path + "/L_" + std::to_string(pos))) throw std::runtime_error("Dataset does not exist: " + mps_path + "/L_" + pos_str);
        if(not h5ppFile.linkExists(mps_path + "/M_" + std::to_string(pos))) throw std::runtime_error("Dataset does not exist: " + mps_path + "/M_" + pos_str);

        auto L = h5ppFile.readDataset<Eigen::Tensor<Scalar, 1>>(mps_path + "/L_" + pos_str);
        auto M = h5ppFile.readDataset<Eigen::Tensor<Scalar, 3>>(mps_path + "/M_" + pos_str);
        state.get_MPS(pos).set_mps(M, L);
        // Sanity checks
        if(pos == position and not state.get_MPS(pos).isCenter())
            throw std::logic_error("Given position is not a a center. State may be wrongly initialized or something is wrong with the resumed file");
        if(pos != position and state.get_MPS(pos).isCenter()) throw std::logic_error("A site not at current position claims to be a state center");
        //        if(passed_LC > 1) throw std::logic_error("Multiple centers encountered");
    }
    if(state.MPS_L.size() + state.MPS_R.size() != model_size) throw std::logic_error("Initialized MPS with the wrong number of sites");
    if(not state.get_MPS(position).isCenter()) throw std::logic_error("Initialized MPS with the center at the wrong position");
    state.set_iter(sim_status.iter);
    state.set_step(sim_status.step);
    state.set_chi_lim(sim_status.chi_lim);
    state.set_chi_max(sim_status.chi_max);
}

void compare(double val1, double val2, double tol, const std::string &tag) {
    if(std::abs(val1 - val2) > tol) tools::log->warn("Value mismatch after resume: {} = {:.16f} | expected {:.16f}", tag, val1, val2);
}

void tools::finite::io::h5resume::validate(const h5pp::File &h5ppFile, const std::string &prefix, const class_state_finite &state,
                                           const class_simulation_status &sim_status) {
    tools::finite::debug::check_integrity(state);
    auto expected_measurements = h5ppFile.readTableEntries<h5pp_table_measurements_finite::table>(prefix + "/measurements");
    state.do_all_measurements();
    compare(state.measurements.energy_variance_per_site.value(), expected_measurements.energy_variance_per_site, 1e-8, "Energy variance per site");
    compare(state.measurements.energy_per_site.value(), expected_measurements.energy_per_site, expected_measurements.energy_variance_per_site,
            "Energy per site");
}



std::vector<ResumeReason> tools::finite::io::h5resume::findResumeReasons (const h5pp::File &h5ppFile, const std::string &sim_name, const std::string &state_prefix,const class_state_finite & state, const class_simulation_status & sim_status){
    // Here we try to find reasons for resuming the previous simulation.
    // The simplest one: it may not have finished

    std::vector<ResumeReason> reasons;
    if(not sim_status.simulation_has_finished) reasons.push_back(ResumeReason::NOT_FINISHED);
    if(not sim_status.simulation_has_converged) reasons.push_back(ResumeReason::NOT_CONVERGED);

    // There may also be
    return reasons;

}
