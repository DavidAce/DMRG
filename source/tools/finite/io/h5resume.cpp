//
// Created by david on 2019-11-07.
//
#include <algorithms/class_algorithm_status.h>
#include <complex>
#include <config/nmspc_settings.h>
#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/finite/debug.h>
#include <tools/finite/env.h>
#include <tools/finite/io.h>
#include <tools/finite/mpo.h>
#include <tools/finite/mps.h>
#include <typeindex>

using Scalar = std::complex<double>;

// Load model, state and simulation status from HDF5
void tools::finite::io::h5resume::load_tensors(const h5pp::File &h5ppFile, const std::string &state_prefix, class_tensors_finite & tensors, class_algorithm_status &status) {
    try {
        tensors.clear_measurements();
        tensors.clear_cache();
        if(h5ppFile.readAttribute<std::string>("storage_level", state_prefix) != enum2str(StorageLevel::FULL))
            throw std::runtime_error("Given prefix to simulation data with StorageLevel < FULL. The simulation can only be resumed from FULL storage");
        tools::common::io::h5resume::load_sim_status_from_hdf5(h5ppFile, state_prefix, status);
        tools::finite::io::h5resume::load_model(h5ppFile, state_prefix, *tensors.model, status);
        tools::finite::io::h5resume::load_state(h5ppFile, state_prefix, *tensors.state, status);
        tools::common::io::h5resume::load_profiling_from_hdf5(h5ppFile, state_prefix);
        tensors.rebuild_edges();
        tools::finite::io::h5resume::validate(h5ppFile, state_prefix, tensors, status);
    } catch(std::exception &ex) {
        throw std::runtime_error("Failed to load simulation from hdf5 file: " + std::string(ex.what()));
    }
}

void tools::finite::io::h5resume::load_model(const h5pp::File &h5ppFile, const std::string &state_prefix, class_model_finite &model,
                                             const class_algorithm_status &status) {
    if(h5ppFile.readAttribute<std::string>("storage_level", state_prefix) != enum2str(StorageLevel::FULL))
        throw std::runtime_error("Given prefix to MPO data with StorageLevel < FULL. The MPO's can only be resumed from FULL storage");
    //    auto position = h5pp_file.readAttribute<size_t>("position", state_prefix);
    //    if(position != status.position)
    //        throw std::runtime_error(fmt::format("Mismatch when loading MPO: State position [{}] != status.position [{}]", position, status.position));

    // Find the path to the MPO
    auto mpo_path = h5ppFile.readAttribute<std::string>(state_prefix, "common/mpo_path");
    auto ham_path = h5ppFile.readAttribute<std::string>(state_prefix, "common/ham_path");
    if(h5ppFile.linkExists(ham_path)) {
        tools::log->trace("Initializing MPOs from Hamiltonian table on file: [{}]", ham_path);
        auto model_type = h5ppFile.readAttribute<std::string>("model_type", ham_path);
        auto model_size = h5ppFile.readAttribute<size_t>("model_size", ham_path);
        if(str2enum<ModelType>(model_type) != settings::model::model_type)
            throw std::runtime_error(fmt::format("Mismatch when loading MPO: model_type [{}] != settings::model::model_type [{}]", model_type,
                                                 enum2str(settings::model::model_type)));
        if(model_size != settings::model::model_size)
            throw std::runtime_error(
                fmt::format("Mismatch when loading MPO: model_size [{}] != settings::model::model_size [{}]", model_type, settings::model::model_size));
        model.initialize(str2enum<ModelType>(model_type),model_size);
        for(size_t pos = 0; pos < model_size; pos++) model.get_mpo(pos).read_hamiltonian(h5ppFile, ham_path);

    } else if(h5ppFile.linkExists(mpo_path)) {
        tools::log->trace("Initializing MPOs from MPO's on file: [{}]", mpo_path);
        auto model_type = h5ppFile.readAttribute<std::string>("model_type", mpo_path);
        auto model_size = h5ppFile.readAttribute<size_t>("model_size", mpo_path);
        if(str2enum<ModelType>(model_type) != settings::model::model_type)
            throw std::runtime_error(fmt::format("Mismatch when loading MPO: Hamiltonian model_type [{}] != settings::model::model_type [{}]", model_type,
                                                 enum2str(settings::model::model_type)));
        if(model_size != settings::model::model_size)
            throw std::runtime_error(
                fmt::format("Mismatch when loading MPO: Hamiltonian sites [{}] != settings::model::model_size [{}]", model_type, settings::model::model_size));
        model.initialize(str2enum<ModelType>(model_type),model_size);
        for(size_t pos = 0; pos < model_size; pos++) model.get_mpo(pos).read_mpo(h5ppFile, mpo_path);
    }
}

void tools::finite::io::h5resume::load_state(const h5pp::File &h5ppFile, const std::string &state_prefix, class_state_finite &state,
                                             const class_algorithm_status &status) {
    if(h5ppFile.readAttribute<std::string>("storage_level", state_prefix) != enum2str(StorageLevel::FULL))
        throw std::runtime_error("Given prefix to MPS data with StorageLevel < FULL. The MPS's can only be resumed from FULL storage");
    std::string mps_path   = state_prefix + "/mps";
    auto        model_type = h5ppFile.readAttribute<std::string>("model_type", state_prefix);
    auto        model_size = h5ppFile.readAttribute<size_t>("model_size", mps_path);
    auto        position   = h5ppFile.readAttribute<size_t>("position", mps_path);
    if(position != status.position)
        throw std::runtime_error(fmt::format("Mismatch when loading MPS: position [{}] != status.position [{}]", position, status.position));
    if(str2enum<ModelType>(model_type) != settings::model::model_type)
        throw std::runtime_error(
            fmt::format("Mismatch when loading MPS: model_type [{}] != settings::model::model_type [{}]", model_type, enum2str(settings::model::model_type)));
    if(model_size != settings::model::model_size)
        throw std::runtime_error(
            fmt::format("Mismatch when loading MPS: model_size [{}] != settings::model::model_size [{}]", model_size, settings::model::model_size));
    state.initialize(str2enum<ModelType>(model_type),model_size,position);
    for(size_t pos = 0; pos < model_size; pos++) {
        std::string pos_str     = std::to_string(pos);
        std::string dset_L_name = fmt::format("{}/{}_{}", mps_path, "L", pos);
        std::string dset_M_name = fmt::format("{}/{}_{}", mps_path, "M", pos);
        if(state.get_mps_site(pos).isCenter()) {
            std::string dset_LC_name = fmt::format("{}/{}", mps_path, "L_C");
            if(not h5ppFile.linkExists(mps_path + "/L_C")) throw std::runtime_error(fmt::format("Dataset does not exist: {}", dset_LC_name));
            auto LC          = h5ppFile.readDataset<Eigen::Tensor<Scalar, 1>>(dset_LC_name);
            auto pos_on_file = h5ppFile.readAttribute<size_t>("position", dset_LC_name);
            if(pos != pos_on_file) throw std::runtime_error(fmt::format("Center bond position mismatch: pos [{}] != pos on file [{}]", pos, pos_on_file));
            state.get_mps_site(pos).set_LC(LC);
        }
        if(not h5ppFile.linkExists(dset_L_name)) throw std::runtime_error(fmt::format("Dataset does not exist: {} ", dset_L_name));
        if(not h5ppFile.linkExists(dset_M_name)) throw std::runtime_error(fmt::format("Dataset does not exist: {} ", dset_M_name));
        auto L = h5ppFile.readDataset<Eigen::Tensor<Scalar, 1>>(dset_L_name);
        auto M = h5ppFile.readDataset<Eigen::Tensor<Scalar, 3>>(dset_M_name);
        state.get_mps_site(pos).set_mps(M, L);
        // Sanity checks
        if(pos == position and not state.get_mps_site(pos).isCenter())
            throw std::logic_error("Given position is not a a center. State may be wrongly initialized or something is wrong with the resumed file");
        if(pos != position and state.get_mps_site(pos).isCenter()) throw std::logic_error("A site not at current position claims to be a state center");
        //        if(passed_LC > 1) throw std::logic_error("Multiple centers encountered");
    }

    state.set_iter(status.iter);
    state.set_step(status.step);
}

void compare(double val1, double val2, double tol, const std::string &tag) {
    if(std::abs(val1 - val2) > tol) tools::log->warn("Value mismatch after resume: {} = {:.16f} | expected {:.16f}", tag, val1, val2);
}

void tools::finite::io::h5resume::validate(const h5pp::File &h5ppFile, const std::string &prefix, class_tensors_finite &tensors,
                                           const class_algorithm_status &status) {
    tools::finite::debug::check_integrity(tensors);
    auto expected_measurements = h5ppFile.readTableEntries<h5pp_table_measurements_finite::table>(prefix + "/measurements");
    tensors.do_all_measurements();
    compare(tensors.measurements.energy_variance_per_site.value(), expected_measurements.energy_variance_per_site, 1e-8, "Energy variance per site");
    compare(tensors.measurements.energy_per_site.value(), expected_measurements.energy_per_site, expected_measurements.energy_variance_per_site,
            "Energy per site");
}

//std::list<SimulationTask> tools::finite::io::h5resume::getTaskList(const h5pp::File &h5ppFile, const std::string &sim_name, const std::string &state_prefix,
//                                                                   const class_tensors_finite &tensors, const class_algorithm_status &status) {
    // Here we try to find reasons for resuming the previous simulation.
    // The simplest one: it may not have finished

//    std::list<SimulationTask> tasks;
    //    if(not status.algorithm_has_succeeded) reasons.push_back(ResumeTasks::FINISH);
    //    if(not status.algorithm_has_converged) reasons.push_back(ResumeTasks::NOT_CONVERGED);
    //    if(status.excited_state_number < settings::xdmrg::max_states) reasons.push_back(ResumeTasks::APPEND_STATES);
    //    if(status.chi_lim_has_reached_chi_max and settings::xdmrg::cfg_chi_lim_max)
    //        ResumeTasks::
    // There may also be
//    return tasks;
//}
