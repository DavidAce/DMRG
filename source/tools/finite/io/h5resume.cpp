//
// Created by david on 2019-11-07.
//
#include <algorithms/class_algorithm_status.h>
#include <complex>
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_exceptions.h>
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
#include <tools/finite/mps.h>
#include <typeindex>

using Scalar = std::complex<double>;

// Load model, state and simulation status from HDF5
void tools::finite::io::h5resume::load_simulation(const h5pp::File &h5ppFile, const std::string &state_prefix, class_tensors_finite &tensors,
                                                  class_algorithm_status &status) {
    try {
        tensors.clear_measurements();
        tensors.clear_cache();
        if(h5ppFile.readAttribute<std::string>(state_prefix, "common/storage_level") != enum2str(StorageLevel::FULL))
            throw std::runtime_error("Given prefix to simulation data with StorageLevel < FULL. The simulation can only be resumed from FULL storage");
        tools::common::io::h5table::load_sim_status(h5ppFile, state_prefix, status);
        tools::common::io::h5table::load_profiling(h5ppFile, state_prefix);
        tools::finite::io::h5resume::load_model(h5ppFile, state_prefix, *tensors.model);
        tools::finite::io::h5resume::load_state(h5ppFile, state_prefix, *tensors.state, status);
        tensors.activate_sites(settings::precision::max_size_full_diag, 2);
        tensors.reduce_mpo_energy();
        tools::finite::io::h5resume::validate(h5ppFile, state_prefix, tensors);
    } catch(const std::exception &ex) { throw except::load_error("Failed to load simulation from hdf5 file: " + std::string(ex.what())); }
}

void tools::finite::io::h5resume::load_model(const h5pp::File &h5ppFile, const std::string &state_prefix, class_model_finite &model) {
    if(h5ppFile.readAttribute<std::string>(state_prefix, "common/storage_level") != enum2str(StorageLevel::FULL))
        throw std::runtime_error("Given prefix to model data with StorageLevel < FULL. The model can only be resumed from FULL storage");

    // Find the path to the MPO
    auto model_prefix = h5ppFile.readAttribute<std::string>(state_prefix, "common/model_prefix");
    if(h5ppFile.linkExists(model_prefix)) {
        auto table_path = fmt::format("{}/hamiltonian", model_prefix);
        if(not h5ppFile.linkExists(table_path)) throw std::runtime_error(fmt::format("Hamiltonian table does not exist: [{}]",table_path));
        tools::log->info("Loading model data from hamiltonian table: [{}]", table_path);
        auto model_type = h5ppFile.readAttribute<std::string>("model_type", table_path);
        auto model_size = h5ppFile.readAttribute<size_t>("model_size", table_path);
        if(str2enum<ModelType>(model_type) != settings::model::model_type)
            throw std::runtime_error(fmt::format("Mismatch when loading model: model_type [{}] != settings::model::model_type [{}]", model_type,
                                                 enum2str(settings::model::model_type)));
        if(model_size != settings::model::model_size)
            throw std::runtime_error(
                fmt::format("Mismatch when loading model: model_size [{}] != settings::model::model_size [{}]", model_size, settings::model::model_size));
        for(size_t pos = 0; pos < model_size; pos++) {
            model.get_mpo(pos).load_hamiltonian(h5ppFile, model_prefix);
        }
    } else {
        throw std::runtime_error(fmt::format("Could not find model data in [{}]", model_prefix));
    }
}

void tools::finite::io::h5resume::load_state(const h5pp::File &h5ppFile, const std::string &state_prefix, class_state_finite &state,
                                             const class_algorithm_status &status) {
    if(h5ppFile.readAttribute<std::string>(state_prefix, "common/storage_level") != enum2str(StorageLevel::FULL))
        throw std::runtime_error("Given prefix to MPS data with StorageLevel < FULL. The MPS's can only be resumed from FULL storage");
    auto mps_prefix = h5ppFile.readAttribute<std::string>(state_prefix, "common/mps_prefix");
    auto model_type = h5ppFile.readAttribute<std::string>(state_prefix, "common/model_type");
    auto model_size = h5ppFile.readAttribute<size_t>(state_prefix, "common/model_size");
    auto position   = h5ppFile.readAttribute<long>(state_prefix, "common/position");
    if(position != status.position)
        throw std::runtime_error(fmt::format("Mismatch when loading MPS: position [{}] != status.position [{}]", position, status.position));
    if(str2enum<ModelType>(model_type) != settings::model::model_type)
        throw std::runtime_error(
            fmt::format("Mismatch when loading MPS: model_type [{}] != settings::model::model_type [{}]", model_type, enum2str(settings::model::model_type)));
    if(model_size != settings::model::model_size)
        throw std::runtime_error(
            fmt::format("Mismatch when loading MPS: model_size [{}] != settings::model::model_size [{}]", model_size, settings::model::model_size));
    state.initialize(str2enum<ModelType>(model_type), model_size, position);
    tools::log->debug("Loading state data from MPS in [{}]", mps_prefix);
    for(size_t pos = 0; pos < model_size; pos++) {
        std::string pos_str     = std::to_string(pos);
        std::string dset_L_name = fmt::format("{}/{}_{}", mps_prefix, "L", pos);
        std::string dset_M_name = fmt::format("{}/{}_{}", mps_prefix, "M", pos);
        if(state.get_mps_site(pos).isCenter()) {
            std::string dset_LC_name = fmt::format("{}/{}", mps_prefix, "L_C");
            if(not h5ppFile.linkExists(dset_LC_name)) throw std::runtime_error(fmt::format("Dataset does not exist: {}", dset_LC_name));
            auto LC          = h5ppFile.readDataset<Eigen::Tensor<Scalar, 1>>(dset_LC_name);
            auto pos_on_file = h5ppFile.readAttribute<size_t>("position", dset_LC_name);
            if(pos != pos_on_file) throw std::runtime_error(fmt::format("Center bond position mismatch: pos [{}] != pos on file [{}]", pos, pos_on_file));
            state.get_mps_site(pos).set_LC(LC);
        }
        if(not h5ppFile.linkExists(dset_L_name)) throw std::runtime_error(fmt::format("Dataset does not exist: {} ", dset_L_name));
        if(not h5ppFile.linkExists(dset_M_name)) throw std::runtime_error(fmt::format("Dataset does not exist: {} ", dset_M_name));
        auto L = h5ppFile.readDataset<Eigen::Tensor<Scalar, 1>>(dset_L_name);
        auto M = h5ppFile.readDataset<Eigen::Tensor<Scalar, 3>>(dset_M_name);
        auto error = h5ppFile.readAttribute<double>("truncation_error", dset_L_name);
        auto label = h5ppFile.readAttribute<std::string>("label", dset_M_name);
        state.get_mps_site(pos).set_mps(M, L, error, label);
        // Sanity checks
        if(static_cast<long>(pos) == position and not state.get_mps_site(pos).isCenter())
            throw std::logic_error("Given position is not a a center. State may be wrongly initialized or something is wrong with the resumed file");
        if(static_cast<long>(pos) != position and state.get_mps_site(pos).isCenter()) throw std::logic_error("A site not at current position claims to be a state center");
        //        if(passed_LC > 1) throw std::logic_error("Multiple centers encountered");
    }
}

void compare(double val1, double val2, double tol, const std::string &tag) {
    if(std::abs(val1 - val2) > tol)
        throw std::runtime_error(fmt::format("Value mismatch after resume: {} = {:.16f} | expected {:.16f} | diff = {:.16f} | tol = {:.16f}", tag, val1, val2,
                                             std::abs(val1 - val2), tol));
    else
        tools::log->debug("{} matches measurement on file: {:.16f} == {:.16f} | diff = {:.16f} | tol = {:.16f}", tag, val1, val2, std::abs(val1 - val2), tol);
}

void tools::finite::io::h5resume::validate(const h5pp::File &h5ppFile, const std::string &state_prefix, class_tensors_finite &tensors) {
    tools::finite::debug::check_integrity(tensors);
    tools::log->debug("Validating resumed state: Measurements should match those on file");
    auto expected_measurements = h5ppFile.readTableRecords<h5pp_table_measurements_finite::table>(state_prefix + "/measurements");
    tensors.do_all_measurements();
    compare(tensors.measurements.energy_variance.value(), expected_measurements.energy_variance, 1e-8, "Energy variance");
    compare(tensors.measurements.energy.value(), expected_measurements.energy, 1e-8, "Energy");
}
