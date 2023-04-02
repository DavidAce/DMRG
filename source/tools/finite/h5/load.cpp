#include "config/settings.h"
#include "debug/exceptions.h"
#include "io/hdf5_types.h"
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
            if(algo_type == AlgorithmType::fLBIT) {
                // To successfully load a simulation there has to be an MPS with StorageLevel::FULL.
                tensors = TensorsFinite(algo_type, settings::model::model_type, settings::model::model_size, 0);
                tools::common::h5::load::pattern(h5file, state_prefix, settings::strategy::initial_pattern);
                tools::common::h5::load::status(h5file, state_prefix, status);
                tools::finite::h5::load::model(h5file, algo_type, *tensors.model);
                tools::common::h5::load::timer(h5file, state_prefix, status);
            } else {
                // Reset tensors
                tensors = TensorsFinite(algo_type, settings::model::model_type, settings::model::model_size, 0);
                MpsInfo info;
                tools::finite::h5::load::state(h5file, state_prefix, *tensors.state, info);
                tools::common::h5::load::status(h5file, state_prefix, status, info);
                tools::finite::h5::load::model(h5file, algo_type, *tensors.model);
                tools::common::h5::load::timer(h5file, state_prefix, status);
                tools::finite::h5::load::validate(h5file, state_prefix, tensors, algo_type);
            }
        } catch(const std::exception &ex) { throw except::load_error("failed to load from state prefix [{}]: {}", state_prefix, ex.what()); }
    }

    void load::state(const h5pp::File &h5file, std::string_view state_prefix, StateFinite &state, MpsInfo &info) {
        using cplx = std::complex<double>;
        try {
            // To successfully load a state there has to be an MPS with StorageLevel::FULL
            auto mps_info = tools::common::h5::resume::find_fully_stored_mps(h5file, state_prefix);
            if(mps_info.empty()) throw except::load_error("No mps with StorageLevel::FULL are stored under [{}]", state_prefix);
            info            = mps_info.back();
            auto model_size = h5file.readAttribute<size_t>(info.pfx, "model_size");
            auto algo_type  = h5file.readAttribute<AlgorithmType>(info.pfx, "algo_type");

            state.mps_sites.clear();
            state.set_algorithm(algo_type);
            for(size_t pos = 0; pos < model_size; ++pos) {
                auto dset_M_name = fmt::format("{}/M_{}", info.pfx, pos);
                auto position    = h5file.readAttribute<long>(dset_M_name, "position");
                auto M           = h5file.readDataset<Eigen::Tensor<cplx, 3>>(dset_M_name);
                auto L           = h5file.readAttribute<Eigen::Tensor<double, 1>>(dset_M_name, "L");
                auto error       = h5file.readAttribute<double>(dset_M_name, "truncation_error");
                auto label       = h5file.readAttribute<std::string>(dset_M_name, "label");
                auto isCenter    = h5file.readAttribute<bool>(dset_M_name, "isCenter");
                state.mps_sites.emplace_back(std::make_unique<MpsSite>(M, L, position, error, label));
                if(isCenter) {
                    auto LC       = h5file.readAttribute<Eigen::Tensor<double, 1>>(dset_M_name, "LC");
                    auto error_LC = h5file.readAttribute<double>(dset_M_name, "truncation_error_LC");
                    state.mps_sites.back()->set_LC(LC, error_LC);
                }
                tools::log->trace("Loaded mps: {}({})", state.mps_sites.back()->get_label(), state.mps_sites.back()->get_position());
                // Sanity checks
                if(num::cmp_equal(pos, position) and not state.mps_sites.back()->isCenter())
                    throw except::logic_error("Given position is not a a center. State may be wrongly initialized or something is wrong with the resumed file");
                if(num::cmp_not_equal(pos, position) and state.mps_sites.back()->isCenter())
                    throw except::logic_error("A site not at current position claims to be a state center");
                //        if(passed_LC > 1) throw std::logic_error("Multiple centers encountered");
            }
            state.set_positions();
        } catch(const std::exception &ex) { throw except::load_error("load state error: {}", ex.what()); }
    }

    void load::model(const h5pp::File &h5file, AlgorithmType algo_type, ModelFinite &model) {
        auto model_path = fmt::format("{}/model", enum2sv(algo_type));
        try {
            // Find the path to the model
            auto table_path = fmt::format("{}/hamiltonian", model_path);
            tools::log->info("Loading model data from hamiltonian table: [{}]", table_path);
            auto model_type = h5file.readAttribute<ModelType>(table_path, "model_type");
            auto model_size = h5file.readAttribute<size_t>(table_path, "model_size");
            if(model_type != settings::model::model_type)
                throw std::runtime_error(
                    fmt::format("model_type [{}] != settings::model::model_type [{}]", enum2sv(model_type), enum2sv(settings::model::model_type)));
            if(model_size != settings::model::model_size)
                throw except::runtime_error("model_size [{}] != settings::model::model_size [{}]", model_size, settings::model::model_size);
            for(const auto &mpo : model.MPO) mpo->load_hamiltonian(h5file, model_path);
            for(const auto &mpo : model.MPO) tools::log->trace("Loaded mpo: {}({})", enum2sv(mpo->model_type), mpo->get_position());
            tools::log->debug("Finished loading model");
            tools::finite::print::model(model);
        } catch(const std::exception &ex) { throw except::load_error("failed to load model [{}]: {}", model_path, ex.what()); }
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
            compare(tensors.state->measurements.entanglement_entropy_midchain.value(), expected_measurements.entanglement_entropy, 1e-8,
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
