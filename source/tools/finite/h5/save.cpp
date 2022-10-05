#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
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
#include "tools/finite/h5.h"
#include "tools/finite/measure.h"
#include "tools/finite/ops.h"
#include <complex>
#include <h5pp/h5pp.h>

namespace tools::finite::h5 {

    void tools::finite::h5::save::measurements(h5pp::File &h5file, const StorageInfo &sinfo, const TensorsFinite &tensors, const AlgorithmStatus &status) {
        save::measurements(h5file, sinfo, *tensors.state, *tensors.model, *tensors.edges, status);
    }

    void save::measurements(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges,
                            const AlgorithmStatus &status) {
        if(sinfo.storage_level <= StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;

        sinfo.assert_well_defined();
        std::string table_path = fmt::format("{}/measurements", sinfo.get_state_prefix());

        // Check if the current entry has already been saved
        auto h5_save_point = tools::common::h5::save::get_last_save_point(h5file, table_path);
        auto save_point    = std::make_pair(sinfo.iter, sinfo.step);
        if(h5_save_point and h5_save_point.value() == save_point) return; // Already saved

        tools::log->trace("Appending to table: {}", table_path);
        auto t_hdf = tid::tic_scope("measurements", tid::level::higher);
        h5pp_table_measurements_finite::register_table_type();
        if(not h5file.linkExists(table_path)) h5file.createTable(h5pp_table_measurements_finite::h5_type, table_path, "measurements");

        h5pp_table_measurements_finite::table measurement_entry{};
        measurement_entry.step     = static_cast<uint64_t>(sinfo.step);
        measurement_entry.iter     = static_cast<uint64_t>(sinfo.iter);
        measurement_entry.position = static_cast<long>(sinfo.position);
        measurement_entry.event    = sinfo.storage_event;
        measurement_entry.length   = static_cast<uint64_t>(tools::finite::measure::length(state));
        measurement_entry.bond_lim = sinfo.bond_lim;
        measurement_entry.bond_max = status.bond_max;
        measurement_entry.bond_mid = static_cast<long>(tools::finite::measure::bond_dimension_midchain(state));
        measurement_entry.norm     = tools::finite::measure::norm(state);
        if(sinfo.algo_type != AlgorithmType::fLBIT) {
            measurement_entry.energy                 = tools::finite::measure::energy(state, model, edges);
            measurement_entry.energy_variance        = tools::finite::measure::energy_variance(state, model, edges);
            measurement_entry.energy_variance_lowest = status.energy_variance_lowest;
        }
        if(status.algo_type == AlgorithmType::fLBIT) { measurement_entry.number_entropy_midchain = tools::finite::measure::number_entropy_midchain(state); }
        measurement_entry.entanglement_entropy_midchain = tools::finite::measure::entanglement_entropy_midchain(state);

        measurement_entry.spin_components  = tools::finite::measure::spin_components(state);
        measurement_entry.truncation_error = state.get_truncation_error_midchain();
        measurement_entry.total_time       = status.wall_time;
        measurement_entry.algorithm_time   = status.algo_time;
        measurement_entry.physical_time    = status.phys_time;

        h5file.appendTableRecords(measurement_entry, table_path);
        h5file.writeAttribute(status.iter, table_path, "iter");
        h5file.writeAttribute(status.step, table_path, "step");
        h5file.writeAttribute(status.bond_lim, table_path, "bond_lim");
        h5file.writeAttribute(status.bond_max, table_path, "bond_max");
    }

    void save::correlations(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.algo_type != AlgorithmType::xDMRG) return;
        if(sinfo.storage_level <= StorageLevel::LIGHT) return;
        if(sinfo.storage_event != StorageEvent::LAST_STATE) return;
        tools::finite::measure::correlation_matrix_xyz(state);
        /* clang-format off */
        if(state.measurements.correlation_matrix_sx) save::data(h5file, sinfo, state.measurements.correlation_matrix_sx.value(), "correlation_matrix_sx", sinfo.get_state_prefix());
        if(state.measurements.correlation_matrix_sy) save::data(h5file, sinfo, state.measurements.correlation_matrix_sy.value(), "correlation_matrix_sy", sinfo.get_state_prefix());
        if(state.measurements.correlation_matrix_sz) save::data(h5file, sinfo, state.measurements.correlation_matrix_sz.value(), "correlation_matrix_sz", sinfo.get_state_prefix());
        /* clang-format on */
    }
    void save::structure_factors(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.algo_type != AlgorithmType::xDMRG) return;
        if(sinfo.storage_level <= StorageLevel::LIGHT) return;
        if(sinfo.storage_event != StorageEvent::LAST_STATE) return;
        tools::finite::measure::structure_factors_xyz(state);
        save::data_as_table(h5file, sinfo, state.measurements.structure_factor_x, "structure_factor_x", "structure factor x", "f");
        save::data_as_table(h5file, sinfo, state.measurements.structure_factor_y, "structure_factor_y", "structure factor y", "f");
        save::data_as_table(h5file, sinfo, state.measurements.structure_factor_z, "structure_factor_z", "structure factor z", "f");
    }

    void save::kvornings_marker(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.algo_type != AlgorithmType::xDMRG) return;
        if(sinfo.storage_level <= StorageLevel::LIGHT) return;
        if(sinfo.storage_event != StorageEvent::LAST_STATE) return;
        tools::finite::measure::kvornings_marker(state);
        save::data_as_table(h5file, sinfo, state.measurements.kvornings_marker, "kvornings_marker", "Kvornings marker", "eigval");
    }

    template<typename T>
    void save::data_as_table(h5pp::File &h5file, const StorageInfo &sinfo, const T *const data, size_t size, std::string_view table_name,
                             std::string_view table_title, std::string_view fieldname) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        auto t_hdf      = tid::tic_scope(table_name);
        auto table_path = fmt::format("{}/{}", sinfo.get_state_prefix(), table_name);
        tools::log->trace("Appending to table: {}", table_path);

        // Check if the current entry has already been appended
        auto h5_save_point = tools::common::h5::save::get_last_save_point(h5file, table_path);
        auto save_point    = std::make_pair(sinfo.iter, sinfo.step);
        if(h5_save_point and h5_save_point.value() == save_point) return;

        // Register the table and create if it doesn't exist
        auto h5_type = h5pp_table_data<T>::register_table_type(h5pp::type::getH5Type<T>(), size, fieldname);
        if(not h5file.linkExists(table_path)) h5file.createTable(h5_type, table_path, table_title);

        // Copy the data into an std::vector<std::byte> stream, which will act as a struct for our table entry
        auto entry = h5pp_table_data<T>::make_entry(sinfo.iter, sinfo.step, sinfo.position, sinfo.storage_event, sinfo.bond_lim, data, size);
        h5file.appendTableRecords(entry, table_path);
        h5file.writeAttribute(sinfo.iter, table_path, "iter");
        h5file.writeAttribute(sinfo.step, table_path, "step");
        h5file.writeAttribute(sinfo.bond_lim, table_path, "bond_lim");
        h5file.writeAttribute(sinfo.bond_max, table_path, "bond_max");
    }

    template<typename T>
    void save::data_as_table_vla(h5pp::File &h5file, const StorageInfo &sinfo, const std::vector<T> &data, const h5pp::hid::h5t &h5elem_t,
                                 std::string_view table_name, std::string_view table_title, std::string_view fieldname) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        auto t_hdf      = tid::tic_scope(table_name);
        auto table_path = fmt::format("{}/{}", sinfo.get_state_prefix(), table_name);
        // Check if the current entry has already been appended
        auto h5_save_point = tools::common::h5::save::get_last_save_point(h5file, table_path);
        auto save_point    = std::make_pair(sinfo.iter, sinfo.step);
        if(h5_save_point and h5_save_point.value() == save_point) return;
        tools::log->trace("Appending to table: {}", table_path);

        // Register the table and create if it doesn't exist
        auto h5_type = h5pp_table_data<T>::register_table_type(h5elem_t, data.size(), fieldname);
        if(not h5file.linkExists(table_path)) h5file.createTable(h5_type, table_path, table_title);

        // Copy the data into an std::vector<std::byte> stream, which will act as a struct for our table entry
        auto entry = h5pp_table_data<T>::make_entry(sinfo.iter, sinfo.step, sinfo.position, sinfo.storage_event, sinfo.bond_lim, data.data(), data.size());
        h5file.appendTableRecords(entry, table_path);
        h5file.writeAttribute(sinfo.iter, table_path, "iter");
        h5file.writeAttribute(sinfo.step, table_path, "step");
        h5file.writeAttribute(sinfo.bond_lim, table_path, "bond_lim");
        h5file.writeAttribute(sinfo.bond_max, table_path, "bond_max");
    }

    void save::bond_dimensions(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        data_as_table(h5file, sinfo, tools::finite::measure::bond_dimensions(state), "bond_dims", "Bond Dimensions", "L_");
    }

    void save::schmidt_values(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        auto                     t_hdf      = tid::tic_scope("schmidt_values", tid::level::higher);
        auto                     table_path = fmt::format("{}/{}", sinfo.get_state_prefix(), "schmidt_values");
        const long               cols       = state.get_length<long>() + 1; // Number of bonds matrices (1d arrays)
        const long               rows       = sinfo.bond_max;               // Max number of schmidt values in each bond (singular values)
        long                     past_LC    = 0;                            // Add one after collecting LC
        Eigen::Tensor<double, 2> bonds(rows, cols);                         // Collect all bond matrices.
        bonds.setZero();
        for(const auto &mps : state.mps_sites) {
            auto offset                 = std::array<long, 2>{0, mps->get_position<long>() + past_LC};
            auto extent                 = std::array<long, 2>{mps->get_L().size(), 1};
            bonds.slice(offset, extent) = mps->get_L().reshape(extent).real();
            if(mps->isCenter()) {
                past_LC                     = 1;
                offset                      = std::array<long, 2>{0, mps->get_position<long>() + past_LC};
                extent                      = std::array<long, 2>{mps->get_LC().size(), 1};
                bonds.slice(offset, extent) = mps->get_LC().reshape(extent).real();
            }
        }
        if(not h5file.linkExists(table_path)) {
            auto                 hrows = static_cast<hsize_t>(rows);
            auto                 hcols = static_cast<hsize_t>(cols);
            std::vector<hsize_t> dims  = {hrows, hcols, 0};
            std::vector<hsize_t> chnk  = {hrows, hcols, 10};
            h5file.createDataset(table_path, h5pp::type::getH5Type<double>(), H5D_CHUNKED, dims, chnk, std::nullopt, 9);
        }
        h5file.appendToDataset(bonds, table_path, 2);
        h5file.writeAttribute(state.get_position<long>(), table_path, "position");
        h5file.writeAttribute(sinfo.iter, table_path, "iter");
        h5file.writeAttribute(sinfo.step, table_path, "step");
        h5file.writeAttribute(sinfo.bond_lim, table_path, "bond_lim");
        h5file.writeAttribute(sinfo.bond_max, table_path, "bond_max");
    }

    void save::truncation_errors(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        auto t_hdf = tid::tic_scope("truncation_errors", tid::level::higher);
        data_as_table(h5file, sinfo, tools::finite::measure::truncation_errors(state), "truncation_errors", "Truncation errors", "L_");
    }

    void save::entropies_neumann(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        auto t_hdf = tid::tic_scope("entropies", tid::level::higher);
        data_as_table(h5file, sinfo, tools::finite::measure::entanglement_entropies(state), "entanglement_entropies", "Entanglement Entropies", "L_");
    }

    void save::number_probabilities(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.algo_type != AlgorithmType::fLBIT) return;
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        if(not state.measurements.number_probabilities) return;
        auto t_hdf      = tid::tic_scope("probabilities", tid::level::higher);
        auto table_path = fmt::format("{}/{}", sinfo.get_state_prefix(), "number_probabilities");
        tools::log->trace("Appending to table: {}", table_path);
        // Check if the current entry has already been appended
        auto h5_save_point = tools::common::h5::save::get_last_save_point(h5file, table_path);
        auto save_point    = std::make_pair(sinfo.iter, sinfo.step);
        if(h5_save_point and h5_save_point.value() == save_point) return;

        if(not h5file.linkExists(table_path)) {
            auto                 rows = static_cast<hsize_t>(state.measurements.number_probabilities->dimension(0));
            auto                 cols = static_cast<hsize_t>(state.measurements.number_probabilities->dimension(1));
            std::vector<hsize_t> dims = {rows, cols, 0};
            std::vector<hsize_t> chnk = {rows, cols, 10};
            h5file.createDataset(table_path, h5pp::type::getH5Type<double>(), H5D_CHUNKED, dims, chnk);
        }
        h5file.appendToDataset(state.measurements.number_probabilities.value(), table_path, 2);
        h5file.writeAttribute(sinfo.iter, table_path, "iter");
        h5file.writeAttribute(sinfo.step, table_path, "step");
        h5file.writeAttribute(sinfo.bond_lim, table_path, "bond_lim");
        h5file.writeAttribute(sinfo.bond_max, table_path, "bond_max");
    }

    void save::entropies_renyi(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level <= StorageLevel::LIGHT) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        auto t_hdf = tid::tic_scope("entropies", tid::level::higher);
        auto inf   = std::numeric_limits<double>::infinity();
        data_as_table(h5file, sinfo, tools::finite::measure::renyi_entropies(state, 2), "renyi_entropies_2", "Renyi Entropy 2", "L_");
        data_as_table(h5file, sinfo, tools::finite::measure::renyi_entropies(state, 3), "renyi_entropies_3", "Renyi Entropy 3", "L_");
        data_as_table(h5file, sinfo, tools::finite::measure::renyi_entropies(state, 4), "renyi_entropies_4", "Renyi Entropy 4", "L_");
        data_as_table(h5file, sinfo, tools::finite::measure::renyi_entropies(state, inf), "renyi_entropies_inf", "Renyi Entropy inf", "L_");
    }

    void save::entropies_number(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level <= StorageLevel::LIGHT) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        if(sinfo.algo_type != AlgorithmType::fLBIT) return;
        auto t_hdf = tid::tic_scope("entropies", tid::level::higher);
        data_as_table(h5file, sinfo, tools::finite::measure::number_entropies(state), "number_entropies", "Number entropies", "L_");
    }

    void save::expectations(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level <= StorageLevel::LIGHT) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        if(sinfo.algo_type != AlgorithmType::xDMRG) return;
        auto t_hdf = tid::tic_scope("expectations", tid::level::higher);
        tools::finite::measure::expectation_values_xyz(state);
        tools::log->trace("Saving expectations to {}", sinfo.get_state_prefix());
        /* clang-format off */
        save::data_as_table(h5file, sinfo, state.measurements.expectation_values_sx, "expectation_values_sx", "<sigma x>", "L_");
        save::data_as_table(h5file, sinfo, state.measurements.expectation_values_sy, "expectation_values_sy", "<sigma y>", "L_");
        save::data_as_table(h5file, sinfo, state.measurements.expectation_values_sz, "expectation_values_sz", "<sigma z>", "L_");
        /* clang-format on */
    }

    void save::bonds(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        auto t_hdf = tid::tic_scope("bonds", tid::level::higher);
        // Transform from cplx to real to save space
        using real      = MpsSite::real;
        auto h5real     = h5pp::type::getH5Type<real>();
        auto bonds_real = std::vector<h5pp::varr_t<real>>();

        std::string column_name = "L_";
        std::string table_title = "Bonds (singular values)";

        if(settings::storage::storage_level_tables == StorageLevel::LIGHT) {
            bonds_real.emplace_back(Eigen::Tensor<real, 1>(state.get_midchain_bond().real()));
            column_name = "L";
            table_title = "Midchain Bond (singular values)";
        } else {
            bonds_real.reserve(state.get_length<size_t>() + 1);
            for(const auto &mps : state.mps_sites) {
                bonds_real.emplace_back(Eigen::Tensor<real, 1>(mps->get_L().real()));
                if(mps->isCenter()) { bonds_real.emplace_back(Eigen::Tensor<real, 1>(mps->get_LC().real())); }
            }
        }
        data_as_table_vla(h5file, sinfo, bonds_real, h5real, "bonds", table_title, column_name);
    }

    void save::state(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level < StorageLevel::FULL) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        auto t_hdf      = tid::tic_scope("state", tid::level::higher);
        auto mps_prefix = sinfo.get_mps_prefix();

        // Check if the current entry has already been saved
        auto h5_save_point_mps = tools::common::h5::save::get_last_save_point(h5file, mps_prefix);
        auto save_point        = std::make_pair(sinfo.iter, sinfo.step);
        if(h5_save_point_mps and h5_save_point_mps.value() == save_point) return;

        tools::log->trace("Storing [{: ^6}]: mps tensors", enum2sv(sinfo.storage_level));
        for(const auto &mps : state.mps_sites) {
            auto dsetName = fmt::format("{}/M_{}", mps_prefix, mps->get_position<long>());
            h5file.writeDataset(mps->get_M_bare(), dsetName, H5D_CHUNKED); // Important to write bare matrices!!
            // TODO: what happens if L is shorter than an existing L? Does it overwrite?
            h5file.deleteAttribute(dsetName, "L");
            h5file.deleteAttribute(dsetName, "LC");
            h5file.deleteAttribute(dsetName, "truncation_error_LC");
            h5file.writeAttribute(Eigen::Tensor<double, 1>(mps->get_L().real()), dsetName, "L");
            h5file.writeAttribute(mps->get_position<long>(), dsetName, "position");
            h5file.writeAttribute(mps->get_M_bare().dimensions(), dsetName, "dimensions");
            h5file.writeAttribute(mps->get_label(), dsetName, "label");
            h5file.writeAttribute(mps->get_truncation_error(), dsetName, "truncation_error");
            h5file.writeAttribute(mps->get_unique_id(), dsetName, "unique_id");
            h5file.writeAttribute(mps->isCenter(), dsetName, "isCenter");
            if(mps->isCenter()) {
                h5file.writeAttribute(mps->get_LC(), dsetName, "LC");
                h5file.writeAttribute(mps->get_truncation_error_LC(), dsetName, "truncation_error_LC");
            }
        }
        h5file.writeAttribute(state.get_length<size_t>(), mps_prefix, "model_size");
        h5file.writeAttribute(state.get_position<long>(), mps_prefix, "position");
        h5file.writeAttribute(sinfo.iter, mps_prefix, "iter");
        h5file.writeAttribute(sinfo.step, mps_prefix, "step");
        h5file.writeAttribute(sinfo.bond_lim, mps_prefix, "bond_lim");
        h5file.writeAttribute(sinfo.bond_max, mps_prefix, "bond_max");
        // Save the storage level as an attribute, to decide if we can resume later
        h5file.writeAttribute(sinfo.algo_type, mps_prefix, "algo_type", std::nullopt, h5_enum_algo_type::get_h5t());
        h5file.writeAttribute(sinfo.storage_level, mps_prefix, "storage_level", std::nullopt, h5_enum_storage_level::get_h5t());
        h5file.writeAttribute(sinfo.storage_event, mps_prefix, "storage_event", std::nullopt, h5_enum_storage_event::get_h5t());
    }

    /*! Write down the Hamiltonian model type and site info as attributes */
    void save::model(h5pp::File &h5file, const StorageInfo &sinfo, const ModelFinite &model) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event != StorageEvent::MODEL) return;
        std::string table_path = fmt::format("{}/model/hamiltonian", sinfo.algo_name);
        if(h5file.linkExists(table_path)) return tools::log->debug("The hamiltonian has already been written to [{}]", table_path);

        tools::log->trace("Storing table: [{}]", table_path);
        auto t_hdf = tid::tic_scope("model", tid::level::higher);
        for(auto site = 0ul; site < model.get_length(); site++) model.get_mpo(site).save_hamiltonian(h5file, table_path);
        h5file.writeAttribute(settings::model::model_size, table_path, "model_size");
        h5file.writeAttribute(settings::model::model_type, table_path, "model_type");
        h5file.writeAttribute(enum2sv(settings::model::model_type), table_path, "model_name");
    }

    /*! Write all the MPO's with site info in attributes */
    void tools::finite::h5::save::mpo(h5pp::File &h5file, const StorageInfo &sinfo, const ModelFinite &model) {
        if(sinfo.storage_level < StorageLevel::FULL) return;
        if(sinfo.storage_event != StorageEvent::MODEL) return;
        std::string mpo_prefix = fmt::format("{}/model/mpo", sinfo.get_state_prefix());
        // We do not expect the MPO's to change. Therefore, if they exist, there is nothing else to do here
        if(h5file.linkExists(mpo_prefix)) return tools::log->trace("The MPO's have already been written to [{}]", mpo_prefix);
        auto t_hdf = tid::tic_scope("mpo", tid::level::higher);
        tools::log->trace("Storing [{: ^6}]: mpos", enum2sv(sinfo.storage_level));
        for(size_t pos = 0; pos < model.get_length(); pos++) { model.get_mpo(pos).save_mpo(h5file, mpo_prefix); }
        h5file.writeAttribute(settings::model::model_size, mpo_prefix, "model_size");
        h5file.writeAttribute(enum2sv(settings::model::model_type), mpo_prefix, "model_type");
    }

    template<typename T>
    void save::data(h5pp::File &h5file, const StorageInfo &sinfo, const T &data, std::string_view data_name, std::string_view prefix) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        auto t_data    = tid::tic_scope("data");
        auto data_path = fmt::format("{}/{}", prefix, data_name);

        // Check if the current entry has already been saved
        auto h5_save_point = tools::common::h5::save::get_last_save_point(h5file, data_path);
        auto save_point    = std::make_pair(sinfo.iter, sinfo.step);
        if(h5_save_point and h5_save_point.value() == save_point) return; // Already saved

        H5D_layout_t layout;
        if(std::is_scalar_v<T>)
            layout = H5D_layout_t::H5D_COMPACT;
        else
            layout = H5D_layout_t::H5D_CHUNKED;
        h5file.writeDataset(data, data_path, layout);
        h5file.writeAttribute(sinfo.iter, data_path, "iter");
        h5file.writeAttribute(sinfo.step, data_path, "step");
        h5file.writeAttribute(sinfo.bond_lim, data_path, "bond_lim");
        h5file.writeAttribute(sinfo.bond_max, data_path, "bond_max");
    }

    template void save::data(h5pp::File &h5file, const StorageInfo &sinfo, const Eigen::Tensor<double, 2> &data, std::string_view data_name,
                             std::string_view prefix);
    template void save::data(h5pp::File &h5file, const StorageInfo &sinfo, const Eigen::Tensor<double, 1> &data, std::string_view data_name,
                             std::string_view prefix);

    template<typename T>
    void save::data(h5pp::File &h5file, const StorageInfo &sinfo, const T &data, std::string_view data_name, CopyPolicy copy_policy) {
        // Setup this save
        auto t_h5     = tid::tic_scope("h5");
        auto t_reason = tid::tic_scope(enum2sv(sinfo.storage_event), tid::level::higher);

        std::string prefix;
        switch(sinfo.storage_event) {
            case StorageEvent::MODEL: {
                prefix = fmt::format("{}/model", sinfo.algo_name);
                break;
            }
            default: {
                prefix = sinfo.get_state_prefix();
                break;
            }
        }

        save::data(h5file, sinfo, data, data_name, prefix);
        tools::common::h5::tmp::copy_from_tmp(h5file, sinfo.iter, sinfo.step, sinfo.storage_event, copy_policy);
    }

    template void save::data(h5pp::File &h5file, const StorageInfo &sinfo, const Eigen::Tensor<std::complex<double>, 2> &data, std::string_view data_name,
                             CopyPolicy copy_policy);

    void save::simulation(h5pp::File &h5file, const TensorsFinite &tensors, const AlgorithmStatus &status, StorageEvent storage_event, CopyPolicy copy_policy) {
        save::simulation(h5file, *tensors.state, *tensors.model, *tensors.edges, status, storage_event, copy_policy);
    }

    void save::simulation(h5pp::File &h5file, const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges, const AlgorithmStatus &status,
                          StorageEvent event, CopyPolicy copy_policy) {
        if(not state.position_is_inward_edge()) return;
        auto sinfo = StorageInfo(status, state.get_name(), event);
        if(sinfo.storage_event == StorageEvent::NONE) return;
        if(sinfo.storage_level == StorageLevel::NONE) return;

        auto t_h5     = tid::tic_scope("h5");
        auto t_sim    = tid::tic_scope("sim");
        auto t_reason = tid::tic_scope(enum2sv(status.event));

        tools::log->debug("Writing to file: Reason [{}] | Level [{}] | state prefix [{}]", enum2sv(sinfo.storage_event), enum2sv(sinfo.storage_level),
                          sinfo.get_state_prefix());

        // The file can be kept open during writes
        h5file.setKeepFileOpened();

        // The main results have now been written. Next we append data to tables
        tools::finite::h5::save::model(h5file, sinfo, model);
        tools::finite::h5::save::state(h5file, sinfo, state);
        tools::common::h5::save::status(h5file, sinfo, status);
        tools::common::h5::save::mem(h5file, sinfo);
        tools::common::h5::save::timer(h5file, sinfo);
        tools::finite::h5::save::measurements(h5file, sinfo, state, model, edges, status);
        tools::finite::h5::save::bond_dimensions(h5file, sinfo, state);
        tools::finite::h5::save::bonds(h5file, sinfo, state);
        //                tools::finite::h5::save::schmidt_values(h5file, sinfo, state);
        tools::finite::h5::save::truncation_errors(h5file, sinfo, state);
        tools::finite::h5::save::entropies_neumann(h5file, sinfo, state);
        tools::finite::h5::save::entropies_renyi(h5file, sinfo, state);
        tools::finite::h5::save::entropies_number(h5file, sinfo, state);
        tools::finite::h5::save::expectations(h5file, sinfo, state);
        tools::finite::h5::save::structure_factors(h5file, sinfo, state);
        tools::finite::h5::save::correlations(h5file, sinfo, state);
        tools::finite::h5::save::kvornings_marker(h5file, sinfo, state);
        tools::finite::h5::save::number_probabilities(h5file, sinfo, state);

        // The file can now be closed
        h5file.setKeepFileClosed();

        // Copy from temporary location to destination depending on given policy
        tools::common::h5::tmp::copy_from_tmp(h5file, sinfo.iter, sinfo.step, sinfo.storage_event, copy_policy);
    }
}