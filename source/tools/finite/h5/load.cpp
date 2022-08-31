#include "config/settings.h"
#include "debug/exceptions.h"
#include "io/table_types.h"
#include "math/num.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/site/mpo/MpoSite.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include "tools/finite/env.h"
#include "tools/finite/h5.h"
#include "tools/finite/print.h"
#include <algorithms/AlgorithmStatus.h>
#include <complex>
#include <h5pp/h5pp.h>
#include <typeindex>

namespace tools::finite::h5 {
    // Load model, state and simulation status from HDF5
    void load::simulation(const h5pp::File &h5file, std::string_view state_prefix, TensorsFinite &tensors, AlgorithmStatus &status, AlgorithmType algo_type) {
        try {
            auto storage_level = h5file.readAttribute<std::string>("common/storage_level", state_prefix);
            if(algo_type == AlgorithmType::fLBIT) {
                if(storage_level == enum2sv(StorageLevel::NONE))
                    throw except::runtime_error(
                        "Given prefix to simulation data with StorageLevel == NONE. The simulation can only be resumed from LIGHT storage or above");

                tensors = TensorsFinite(algo_type, settings::model::model_type, settings::model::model_size, 0);
                tools::common::h5::load::status(h5file, state_prefix, status);
                tools::finite::h5::load::model(h5file, state_prefix, *tensors.model);
                tools::common::h5::load::timer(h5file, state_prefix, status);
            } else {
                if(storage_level != enum2sv(StorageLevel::FULL))
                    throw except::runtime_error(
                        "Given prefix to simulation data with StorageLevel < FULL. The simulation can only be resumed from FULL storage");

                // Reset tensors
                tensors = TensorsFinite(algo_type, settings::model::model_type, settings::model::model_size, 0);

                tools::common::h5::load::status(h5file, state_prefix, status);
                tools::finite::h5::load::model(h5file, state_prefix, *tensors.model);
                tools::finite::h5::load::state(h5file, state_prefix, *tensors.state, status);
                tools::common::h5::load::timer(h5file, state_prefix, status);
                tools::finite::h5::load::validate(h5file, state_prefix, tensors, algo_type);
            }

        } catch(const std::exception &ex) { throw except::load_error("failed to load from state prefix [{}]: {}", state_prefix, ex.what()); }
    }

    void load::model(const h5pp::File &h5file, std::string_view state_prefix, ModelFinite &model) {
        auto model_prefix = h5file.readAttribute<std::string>("common/model_prefix", state_prefix);
        try {
            if(h5file.readAttribute<std::string>("common/storage_level", state_prefix) == enum2sv(StorageLevel::NONE))
                throw std::runtime_error("Given prefix to model data with StorageLevel == NONE. The model can only be resumed from LIGHT storage or above");
            // Find the path to the model
            if(h5file.linkExists(model_prefix)) {
                auto table_path = fmt::format("{}/hamiltonian", model_prefix);
                if(not h5file.linkExists(table_path)) throw except::runtime_error("Hamiltonian table does not exist: [{}]", table_path);
                tools::log->info("Loading model data from hamiltonian table: [{}]", table_path);
                auto model_type = h5file.readAttribute<std::string>(table_path, "model_type");
                auto model_size = h5file.readAttribute<size_t>(table_path, "model_size");
                if(sv2enum<ModelType>(model_type) != settings::model::model_type)
                    throw std::runtime_error(
                        fmt::format("model_type [{}] != settings::model::model_type [{}]", model_type, enum2sv(settings::model::model_type)));
                if(model_size != settings::model::model_size)
                    throw except::runtime_error("model_size [{}] != settings::model::model_size [{}]", model_size, settings::model::model_size);
                for(const auto &mpo : model.MPO) mpo->load_hamiltonian(h5file, model_prefix);
                for(const auto &mpo : model.MPO) tools::log->trace("Loaded mpo: {}({})", enum2sv(mpo->model_type), mpo->get_position());
            } else {
                throw std::runtime_error("Could not find model data");
            }
            tools::log->debug("Finished loading model");
            tools::finite::print::model(model);
        } catch(const std::exception &ex) { throw except::runtime_error("model error [{}]: {}", model_prefix, ex.what()); }
    }

    void load::state(const h5pp::File &h5file, std::string_view state_prefix, StateFinite &state, const AlgorithmStatus &status) {
        using cplx = std::complex<double>;
        try {
            if(h5file.readAttribute<std::string>("common/storage_level", state_prefix) != enum2sv(StorageLevel::FULL))
                throw std::runtime_error("Given prefix to MPS data with StorageLevel < FULL. The MPS's can only be resumed from FULL storage");
            auto bonds_prefix = h5file.readAttribute<std::string>("common/bonds_prefix", state_prefix);
            auto mps_prefix   = h5file.readAttribute<std::string>("common/mps_prefix", state_prefix);
            auto model_type = h5file.readAttribute<std::string>("common/model_type", state_prefix);
            auto model_size = h5file.readAttribute<size_t>("common/model_size", state_prefix);
            auto position   = h5file.readAttribute<long>("common/position", state_prefix);
            if(position != status.position) throw except::runtime_error("MPS: position [{}] != status.position [{}]", position, status.position);
            if(sv2enum<ModelType>(model_type) != settings::model::model_type)
                throw except::runtime_error("MPS: model_type [{}] != settings::model::model_type [{}]", model_type, enum2sv(settings::model::model_type));
            if(model_size != settings::model::model_size)
                throw except::runtime_error("MPS: model_size [{}] != settings::model::model_size [{}]", model_size, settings::model::model_size);
            state.initialize(status.algo_type, sv2enum<ModelType>(model_type), model_size, static_cast<size_t>(position));
            tools::log->debug("Loading state data from MPS in [{}]", mps_prefix);
            for(const auto &mps : state.mps_sites) {
                auto        pos         = mps->get_position<long>();
                std::string pos_str     = std::to_string(pos);
                std::string dset_L_name = fmt::format("{}/{}_{}", bonds_prefix, "L", pos);
                std::string dset_M_name = fmt::format("{}/{}_{}", mps_prefix, "M", pos);
                if(mps->isCenter()) {
                    std::string dset_LC_name = fmt::format("{}/{}", bonds_prefix, "L_C");
                    if(not h5file.linkExists(dset_LC_name)) throw except::runtime_error("Dataset does not exist: {}", dset_LC_name);
                    auto LC          = h5file.readDataset<Eigen::Tensor<cplx, 1>>(dset_LC_name);
                    auto pos_on_file = h5file.readAttribute<long>(dset_LC_name, "position");
                    auto error       = h5file.readAttribute<double>(dset_LC_name, "truncation_error");
                    if(pos != pos_on_file) throw except::runtime_error("Center bond position mismatch: pos [{}] != pos on file [{}]", pos, pos_on_file);
                    mps->set_LC(LC, error);
                }
                if(not h5file.linkExists(dset_L_name)) throw except::runtime_error("Dataset does not exist: {} ", dset_L_name);
                if(not h5file.linkExists(dset_M_name)) throw except::runtime_error("Dataset does not exist: {} ", dset_M_name);
                auto L     = h5file.readDataset<Eigen::Tensor<cplx, 1>>(dset_L_name);
                auto M     = h5file.readDataset<Eigen::Tensor<cplx, 3>>(dset_M_name);
                auto error = h5file.readAttribute<double>(dset_L_name, "truncation_error");
                auto label = h5file.readAttribute<std::string>(dset_M_name, "label");
                mps->set_mps(M, L, error, label);
                tools::log->trace("Loaded mps: {}({})", mps->get_label(), mps->get_position());
                // Sanity checks
                if(num::cmp_equal(pos, position) and not mps->isCenter())
                    throw except::logic_error("Given position is not a a center. State may be wrongly initialized or something is wrong with the resumed file");
                if(num::cmp_not_equal(pos, position) and mps->isCenter())
                    throw except::logic_error("A site not at current position claims to be a state center");
                //        if(passed_LC > 1) throw std::logic_error("Multiple centers encountered");
            }
        } catch(const std::exception &ex) { throw except::runtime_error("state error: {}", ex.what()); }
    }

    void compare(double val1, double val2, double tol, std::string_view tag) {
        if(val2 == 0) {
            tools::log->warn("Value mismatch after resume: {} = {:.16f} | expected {:.16f} | diff = {:.16f} | tol = {:.16f}", tag, val1, val2,
                             std::abs(val1 - val2), tol);
            tools::log->warn("Value 2 is zero: It may not have been initialized in the previous simulation");
            return;
        }
        if(std::abs(val1 - val2) > tol)
            throw except::runtime_error("Value mismatch after resume: {} = {:.16f} | expected {:.16f} | diff = {:.16f} | tol = {:.16f}", tag, val1, val2,
                                        std::abs(val1 - val2), tol);
        else
            tools::log->debug("{} matches measurement on file: {:.16f} == {:.16f} | diff = {:.16f} | tol = {:.16f}", tag, val1, val2, std::abs(val1 - val2),
                              tol);
    }

    void load::validate(const h5pp::File &h5file, std::string_view state_prefix, TensorsFinite &tensors, AlgorithmType algo_type) {
        auto t_val             = tid::tic_scope("validate");
        auto table_prefix      = h5file.readAttribute<std::vector<std::string>>("common/table_prfxs", state_prefix).front();
        auto measurements_path = fmt::format("{}/measurements", table_prefix);
        tensors.rebuild_mpo();
        tensors.rebuild_mpo_squared();
        tensors.activate_sites({tensors.get_position<size_t>()});
        tensors.rebuild_edges();

        if(algo_type == AlgorithmType::fLBIT) {
            // In this case we have loaded state_real from file.
            // However, the MPO's belong to state_lbit, so measuring the energy on state_real w.r.t the lbit-hamiltonian
            // makes no sense.
            // Therefore, we have to validate using entropy, or other state-specific measurements that are independent of the hamlitonian.
            auto expected_measurements = h5file.readTableRecords<h5pp_table_measurements_finite::table>(measurements_path);
            tools::log->debug("Validating resumed state: [{}]", state_prefix);
            tools::log->debug("State labels: {}", tensors.state->get_labels());
            tensors.state->clear_cache();
            tensors.state->clear_measurements();
            tensors.state->do_all_measurements();
            compare(tensors.state->measurements.entanglement_entropy_midchain.value(), expected_measurements.entanglement_entropy_midchain, 1e-8,
                    "Entanglement entropy");
        } else {
            auto expected_measurements = h5file.readTableRecords<h5pp_table_measurements_finite::table>(measurements_path);
            tools::log->debug("Validating resumed state (without energy reduction): [{}]", state_prefix);
            tools::log->debug("State labels: {}", tensors.state->get_labels());
            tensors.clear_measurements();
            tensors.do_all_measurements();
            compare(tensors.measurements.energy.value(), expected_measurements.energy, 1e-8, "Energy");
            compare(tensors.measurements.energy_variance.value(), expected_measurements.energy_variance, 1e-8, "Energy variance");

            if(settings::precision::use_mpo_energy_shift) {
                tensors.shift_mpo_energy();
                tensors.rebuild_mpo();
                tensors.rebuild_mpo_squared();
                tensors.rebuild_edges();
                tools::log->debug("Validating resumed state (after energy reduction): [{}]", state_prefix);
                tools::log->debug("State labels: {}", tensors.state->get_labels());
                tensors.clear_measurements();
                tensors.do_all_measurements();
                compare(tensors.measurements.energy.value(), expected_measurements.energy, 1e-8, "Energy");
                compare(tensors.measurements.energy_variance.value(), expected_measurements.energy_variance, 1e-8, "Energy variance");
            }
        }
    }
}
