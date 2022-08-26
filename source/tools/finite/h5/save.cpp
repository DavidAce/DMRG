#include "algorithms/AlgorithmStatus.h"
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
#include "tools/finite/h5.h"
#include "tools/finite/measure.h"
#include "tools/finite/ops.h"
#include <complex>
#include <h5pp/h5pp.h>

namespace tools::finite::h5 {

    void tools::finite::h5::save::measurements(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level,
                                               const TensorsFinite &tensors, const AlgorithmStatus &status) {
        save::measurements(h5file, table_prefix, storage_level, *tensors.state, *tensors.model, *tensors.edges, status);
    }

    void save::measurements(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level, const StateFinite &state,
                            const ModelFinite &model, const EdgesFinite &edges, const AlgorithmStatus &status) {
        if(storage_level == StorageLevel::NONE) return;
        std::string table_path = fmt::format("{}/measurements", table_prefix);

        // Check if the current entry has already been saved
        auto h5_save_point = tools::common::h5::save::get_last_save_point(h5file, table_path);
        auto save_point    = std::make_pair(status.iter, status.step);
        if(h5_save_point and h5_save_point.value() == save_point) return; // Already saved

        tools::log->trace("Appending to table: {}", table_path);
        auto t_hdf = tid::tic_scope("measurements", tid::level::extra);
        h5pp_table_measurements_finite::register_table_type();
        if(not h5file.linkExists(table_path)) h5file.createTable(h5pp_table_measurements_finite::h5_type, table_path, "measurements");

        h5pp_table_measurements_finite::table measurement_entry{};
        measurement_entry.step     = static_cast<uint64_t>(status.step);
        measurement_entry.iter     = static_cast<uint64_t>(status.iter);
        measurement_entry.position = static_cast<long>(status.position);
        measurement_entry.length   = static_cast<uint64_t>(tools::finite::measure::length(state));
        measurement_entry.bond_lim = status.bond_lim;
        measurement_entry.bond_max = status.bond_max;
        measurement_entry.bond_mid = static_cast<long>(tools::finite::measure::bond_dimension_midchain(state));
        measurement_entry.norm     = tools::finite::measure::norm(state);
        if(status.algo_type != AlgorithmType::fLBIT) {
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

    void save::correlations(h5pp::File &h5file, std::string_view state_prefix, const StorageLevel &storage_level, const StateFinite &state,
                            const AlgorithmStatus &status) {
        if(storage_level <= StorageLevel::LIGHT) return;
        if(status.algo_type != AlgorithmType::xDMRG) return;
        auto t_hdf = tid::tic_scope("correlations", tid::level::extra);
        tools::finite::measure::correlation_matrix_xyz(state);
        tools::log->trace("Saving correlations to {}", state_prefix);
        /* clang-format off */
        if(state.measurements.correlation_matrix_sx) save::data(h5file, state.measurements.correlation_matrix_sx.value(), "correlation_matrix_sx", state_prefix, storage_level, status);
        if(state.measurements.correlation_matrix_sy) save::data(h5file, state.measurements.correlation_matrix_sy.value(), "correlation_matrix_sy", state_prefix, storage_level, status);
        if(state.measurements.correlation_matrix_sz) save::data(h5file, state.measurements.correlation_matrix_sz.value(), "correlation_matrix_sz", state_prefix, storage_level, status);
        /* clang-format on */
    }
    void save::structure_factors(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level, const StateFinite &state,
                                 const AlgorithmStatus &status) {
        if(storage_level <= StorageLevel::LIGHT) return;
        if(status.algo_type != AlgorithmType::xDMRG) return;
        auto t_hdf = tid::tic_scope("structure_factors", tid::level::extra);
        tools::finite::measure::structure_factors_xyz(state);
        tools::log->trace("Saving structure factors to {}", table_prefix);
        /* clang-format off */
        save::data_as_table(h5file, table_prefix, status, state.measurements.structure_factor_x, "structure_factor_x", "structure factor x", "f");
        save::data_as_table(h5file, table_prefix, status, state.measurements.structure_factor_y, "structure_factor_y", "structure factor y", "f");
        save::data_as_table(h5file, table_prefix, status, state.measurements.structure_factor_z, "structure_factor_z", "structure factor z", "f");
        /* clang-format on */
    }

    void save::kvornings_marker(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level, const StateFinite &state,
                                const AlgorithmStatus &status) {
        if(storage_level <= StorageLevel::LIGHT) return;
        if(status.algo_type != AlgorithmType::xDMRG) return;
        auto t_hdf = tid::tic_scope("kvornings_marker", tid::level::extra);
        tools::finite::measure::kvornings_marker(state);
        tools::log->trace("Saving kvornings marker to {}", table_prefix);
        /* clang-format off */
        save::data_as_table(h5file, table_prefix, status, state.measurements.kvornings_marker, "kvornings_marker", "Kvornings marker", "eigval");
        /* clang-format on */
    }

    template<typename T>
    void save::data_as_table(h5pp::File &h5file, std::string_view table_prefix, const AlgorithmStatus &status, const T *const data, size_t size,
                             std::string_view table_name, std::string_view table_title, std::string_view fieldname) {
        auto table_path = fmt::format("{}/{}", table_prefix, table_name);
        tools::log->trace("Appending to table: {}", table_path);
        // Check if the current entry has already been appended
        auto h5_save_point = tools::common::h5::save::get_last_save_point(h5file, table_path);
        auto save_point    = std::make_pair(status.iter, status.step);
        if(h5_save_point and h5_save_point.value() == save_point) return;

        // Register the table and create if it doesnt exist
        auto h5_type = h5pp_table_data<T>::register_table_type(size, fieldname);
        if(not h5file.linkExists(table_path)) h5file.createTable(h5_type, table_path, table_title);

        // Copy the data into an std::vector<std::byte> stream, which will act as a struct for our table entry
        auto entry = h5pp_table_data<T>::make_entry(status.iter, status.step, status.bond_lim, data, size);
        h5file.appendTableRecords(entry, table_path);
        h5file.writeAttribute(status.iter, table_path, "iter");
        h5file.writeAttribute(status.step, table_path, "step");
        h5file.writeAttribute(status.bond_lim, table_path, "bond_lim");
        h5file.writeAttribute(status.bond_max, table_path, "bond_max");
    }

    void save::bond_dimensions(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level, const StateFinite &state,
                               const AlgorithmStatus &status) {
        if(storage_level == StorageLevel::NONE) return;
        auto t_hdf = tid::tic_scope("bond_dimensions", tid::level::extra);
        data_as_table(h5file, table_prefix, status, tools::finite::measure::bond_dimensions(state), "bond_dims", "Bond Dimensions", "L_");
    }

    void save::schmidt_values(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level, const StateFinite &state,
                              const AlgorithmStatus &status) {
        if(storage_level == StorageLevel::NONE) return;
        auto                     t_hdf      = tid::tic_scope("schmidt_values", tid::level::extra);
        auto                     table_path = fmt::format("{}/{}", table_prefix, "schmidt_values");
        const long               cols       = state.get_length<long>() + 1; // Number of bonds matrices (1d arrays)
        const long               rows       = status.bond_max;              // Max number of schmidt values in each bond (singular values)
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
            h5file.createDataset(table_path, h5pp::util::getH5Type<double>(), H5D_CHUNKED, dims, chnk, std::nullopt, 9);
        }
        h5file.appendToDataset(bonds, table_path, 2);
        h5file.writeAttribute(state.get_position<long>(), table_path, "position");
        h5file.writeAttribute(status.iter, table_path, "iter");
        h5file.writeAttribute(status.step, table_path, "step");
        h5file.writeAttribute(status.bond_lim, table_path, "bond_lim");
        h5file.writeAttribute(status.bond_max, table_path, "bond_max");
    }

    void save::truncation_errors(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level, const StateFinite &state,
                                 const AlgorithmStatus &status) {
        if(storage_level == StorageLevel::NONE) return;
        auto t_hdf = tid::tic_scope("truncation_errors", tid::level::extra);
        data_as_table(h5file, table_prefix, status, tools::finite::measure::truncation_errors(state), "truncation_errors", "Truncation errors", "L_");
    }

    void save::entropies_neumann(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level, const StateFinite &state,
                                 const AlgorithmStatus &status) {
        if(storage_level == StorageLevel::NONE) return;
        auto t_hdf = tid::tic_scope("entropies", tid::level::extra);
        data_as_table(h5file, table_prefix, status, tools::finite::measure::entanglement_entropies(state), "entanglement_entropies", "Entanglement Entropies",
                      "L_");
    }

    void save::number_probabilities(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level, const StateFinite &state,
                                    const AlgorithmStatus &status) {
        if(storage_level == StorageLevel::NONE) return;
        if(not state.measurements.number_probabilities) return;
        auto t_hdf      = tid::tic_scope("probabilities", tid::level::extra);
        auto table_path = fmt::format("{}/{}", table_prefix, "number_probabilities");
        tools::log->trace("Appending to table: {}", table_path);
        // Check if the current entry has already been appended
        auto h5_save_point = tools::common::h5::save::get_last_save_point(h5file, table_path);
        auto save_point    = std::make_pair(status.iter, status.step);
        if(h5_save_point and h5_save_point.value() == save_point) return;

        if(not h5file.linkExists(table_path)) {
            auto                 rows = static_cast<hsize_t>(state.measurements.number_probabilities->dimension(0));
            auto                 cols = static_cast<hsize_t>(state.measurements.number_probabilities->dimension(1));
            std::vector<hsize_t> dims = {rows, cols, 0};
            std::vector<hsize_t> chnk = {rows, cols, 10};
            h5file.createDataset(table_path, h5pp::util::getH5Type<double>(), H5D_CHUNKED, dims, chnk);
        }
        h5file.appendToDataset(state.measurements.number_probabilities.value(), table_path, 2);
        h5file.writeAttribute(status.iter, table_path, "iter");
        h5file.writeAttribute(status.step, table_path, "step");
        h5file.writeAttribute(status.bond_lim, table_path, "bond_lim");
        h5file.writeAttribute(status.bond_max, table_path, "bond_max");
    }

    void save::entropies_renyi(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level, const StateFinite &state,
                               const AlgorithmStatus &status) {
        if(storage_level <= StorageLevel::LIGHT) return;
        auto t_hdf = tid::tic_scope("entropies", tid::level::extra);
        auto inf   = std::numeric_limits<double>::infinity();
        data_as_table(h5file, table_prefix, status, tools::finite::measure::renyi_entropies(state, 2), "renyi_entropies_2", "Renyi Entropy 2", "L_");
        data_as_table(h5file, table_prefix, status, tools::finite::measure::renyi_entropies(state, 3), "renyi_entropies_3", "Renyi Entropy 3", "L_");
        data_as_table(h5file, table_prefix, status, tools::finite::measure::renyi_entropies(state, 4), "renyi_entropies_4", "Renyi Entropy 4", "L_");
        data_as_table(h5file, table_prefix, status, tools::finite::measure::renyi_entropies(state, inf), "renyi_entropies_inf", "Renyi Entropy inf", "L_");
    }

    void save::entropies_number(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level, const StateFinite &state,
                                const AlgorithmStatus &status) {
        if(storage_level == StorageLevel::NONE) return;
        if(status.algo_type != AlgorithmType::fLBIT) return;

        auto t_hdf = tid::tic_scope("entropies", tid::level::extra);
        data_as_table(h5file, table_prefix, status, tools::finite::measure::number_entropies(state), "number_entropies", "Number entropies", "L_");
    }

    void save::expectations(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level, const StateFinite &state,
                            const AlgorithmStatus &status) {
        if(storage_level <= StorageLevel::LIGHT) return;
        if(status.algo_type != AlgorithmType::xDMRG) return;
        auto t_hdf = tid::tic_scope("expectations", tid::level::extra);
        tools::finite::measure::expectation_values_xyz(state);
        tools::log->trace("Saving expectations to {}", table_prefix);
        /* clang-format off */
        save::data_as_table(h5file, table_prefix, status, state.measurements.expectation_values_sx, "expectation_values_sx", "<sigma x>", "L_");
        save::data_as_table(h5file, table_prefix, status, state.measurements.expectation_values_sy, "expectation_values_sy", "<sigma y>", "L_");
        save::data_as_table(h5file, table_prefix, status, state.measurements.expectation_values_sz, "expectation_values_sz", "<sigma z>", "L_");
        /* clang-format on */
    }

    void save::state(h5pp::File &h5file, std::string_view state_prefix, const StorageLevel &storage_level, const StateFinite &state,
                     const AlgorithmStatus &status) {
        if(storage_level <= StorageLevel::LIGHT) return;
        auto t_hdf            = tid::tic_scope("state", tid::level::extra);
        auto dsetname_schmidt = fmt::format("{}/schmidt_midchain", state_prefix);
        auto mps_prefix       = fmt::format("{}/mps", state_prefix);

        // Check if the current entry has already been saved
        auto h5_save_point_dset_schmidt = tools::common::h5::save::get_last_save_point(h5file, dsetname_schmidt);
        auto h5_save_point_mps          = tools::common::h5::save::get_last_save_point(h5file, mps_prefix);
        auto save_point                 = std::make_pair(status.iter, status.step);

        bool skip_dset_schmidt = h5_save_point_dset_schmidt and h5_save_point_dset_schmidt.value() == save_point;
        bool skip_mps          = h5_save_point_mps and h5_save_point_mps.value() == save_point;

        if(not skip_dset_schmidt) {
            /*! Writes down the midchain "Lambda" bond matrix (singular values). */
            tools::log->trace("Storing [{: ^6}]: mid bond matrix", enum2sv(storage_level));
            h5file.writeDataset(state.midchain_bond(), dsetname_schmidt, H5D_CHUNKED);
            h5file.writeAttribute(state.get_truncation_error_midchain(), dsetname_schmidt, "truncation_error");
            h5file.writeAttribute((state.get_length<long>() - 1) / 2, dsetname_schmidt, "position");
            h5file.writeAttribute(status.iter, dsetname_schmidt, "iter");
            h5file.writeAttribute(status.step, dsetname_schmidt, "step");
            h5file.writeAttribute(status.bond_lim, dsetname_schmidt, "bond_lim");
            h5file.writeAttribute(status.bond_max, dsetname_schmidt, "bond_max");
        }

        if(not skip_mps) {
            tools::log->trace("Storing [{: ^6}]: bond matrices", enum2sv(storage_level));
            // There should be one more sites+1 number of L's, because there is also a center bond
            // However L_i always belongs M_i. Stick to this rule!
            // This means that some M_i has two bonds, one L_i to the left, and one L_C to the right.
            for(const auto &mps : state.mps_sites) {
                auto dsetName = fmt::format("{}/L_{}", mps_prefix, mps->get_position<long>());
                h5file.writeDataset(mps->get_L(), dsetName, H5D_CHUNKED);
                h5file.writeAttribute(mps->get_position<long>(), dsetName, "position");
                h5file.writeAttribute(mps->get_L().dimensions(), dsetName, "dimensions");
                h5file.writeAttribute(mps->get_truncation_error(), dsetName, "truncation_error");
                if(mps->isCenter()) {
                    dsetName = fmt::format("{}/L_C", mps_prefix);
                    h5file.writeDataset(mps->get_LC(), dsetName, H5D_CHUNKED);
                    h5file.writeAttribute(mps->get_position<long>(), dsetName, "position");
                    h5file.writeAttribute(mps->get_LC().dimensions(), dsetName, "dimensions");
                    h5file.writeAttribute(mps->get_truncation_error_LC(), dsetName, "truncation_error");
                }
            }
            h5file.writeAttribute(state.get_length(), mps_prefix, "model_size");
            h5file.writeAttribute(state.get_position<long>(), mps_prefix, "position");
            h5file.writeAttribute(state.get_truncation_errors(), mps_prefix, "truncation_errors");
            h5file.writeAttribute(state.get_labels(), mps_prefix, "labels");
            h5file.writeAttribute(status.iter, mps_prefix, "iter");
            h5file.writeAttribute(status.step, mps_prefix, "step");
            h5file.writeAttribute(status.bond_lim, mps_prefix, "bond_lim");
            h5file.writeAttribute(status.bond_max, mps_prefix, "bond_max");
        }

        /*! Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation. */
        if(storage_level < StorageLevel::FULL) { return; }

        if(not skip_mps) {
            tools::log->trace("Storing [{: ^6}]: mps tensors", enum2sv(storage_level));
            for(const auto &mps : state.mps_sites) {
                auto dsetName = fmt::format("{}/M_{}", mps_prefix, mps->get_position<long>());
                h5file.writeDataset(mps->get_M_bare(), dsetName, H5D_CHUNKED); // Important to write bare matrices!!
                h5file.writeAttribute(mps->get_position<long>(), dsetName, "position");
                h5file.writeAttribute(mps->get_M_bare().dimensions(), dsetName, "dimensions");
                h5file.writeAttribute(mps->get_label(), dsetName, "label");
                h5file.writeAttribute(mps->get_unique_id(), dsetName, "unique_id");
            }
        }
    }

    /*! Write down the Hamiltonian model type and site info as attributes */
    void save::model(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level, const ModelFinite &model) {
        if(storage_level == StorageLevel::NONE) return;
        std::string table_path = fmt::format("{}/hamiltonian", table_prefix);
        if(h5file.linkExists(table_path)) return tools::log->debug("The hamiltonian has already been written to [{}]", table_path);

        tools::log->trace("Storing table: [{}]", table_path);
        auto t_hdf = tid::tic_scope("model", tid::level::extra);
        for(auto site = 0ul; site < model.get_length(); site++) model.get_mpo(site).save_hamiltonian(h5file, table_path);
        h5file.writeAttribute(enum2sv(settings::model::model_type), table_path, "model_type");
        h5file.writeAttribute(settings::model::model_size, table_path, "model_size");
    }

    /*! Write all the MPO's with site info in attributes */
    void tools::finite::h5::save::mpo(h5pp::File &h5file, std::string_view model_prefix, const StorageLevel &storage_level, const ModelFinite &model) {
        if(storage_level < StorageLevel::FULL) return;
        std::string mpo_prefix = fmt::format("{}/mpo", model_prefix);
        // We do not expect the MPO's to change. Therefore if they exist, there is nothing else to do here
        if(h5file.linkExists(mpo_prefix)) return tools::log->trace("The MPO's have already been written to [{}]", mpo_prefix);
        auto t_hdf = tid::tic_scope("mpo", tid::level::extra);
        tools::log->trace("Storing [{: ^6}]: mpos", enum2sv(storage_level));
        for(size_t pos = 0; pos < model.get_length(); pos++) { model.get_mpo(pos).save_mpo(h5file, mpo_prefix); }
        h5file.writeAttribute(settings::model::model_size, mpo_prefix, "model_size");
        h5file.writeAttribute(enum2sv(settings::model::model_type), mpo_prefix, "model_type");
    }

    template<typename T>
    void save::data(h5pp::File &h5file, const T &data, std::string_view data_name, std::string_view prefix, const StorageLevel &storage_level,
                    const AlgorithmStatus &status) {
        if(storage_level == StorageLevel::NONE) return;
        auto t_data    = tid::tic_scope("data");
        auto data_path = fmt::format("{}/{}", prefix, data_name);

        // Check if the current entry has already been saved
        auto h5_save_point = tools::common::h5::save::get_last_save_point(h5file, data_path);
        auto save_point    = std::make_pair(status.iter, status.step);
        if(h5_save_point and h5_save_point.value() == save_point) return; // Already saved

        H5D_layout_t layout;
        if(std::is_scalar_v<T>)
            layout = H5D_layout_t::H5D_COMPACT;
        else
            layout = H5D_layout_t::H5D_CHUNKED;
        h5file.writeDataset(data, data_path, layout);
        h5file.writeAttribute(status.iter, data_path, "iter");
        h5file.writeAttribute(status.step, data_path, "step");
        h5file.writeAttribute(status.bond_lim, data_path, "bond_lim");
        h5file.writeAttribute(status.bond_max, data_path, "bond_max");
    }

    template void save::data(h5pp::File &h5file, const Eigen::Tensor<double, 2> &data, std::string_view data_name, std::string_view prefix,
                             const StorageLevel &storage_level, const AlgorithmStatus &status);
    template void save::data(h5pp::File &h5file, const Eigen::Tensor<double, 1> &data, std::string_view data_name, std::string_view prefix,
                             const StorageLevel &storage_level, const AlgorithmStatus &status);

    template<typename T>
    void save::data(h5pp::File &h5file, const T &data, std::string_view data_name, std::string_view state_name, const AlgorithmStatus &status,
                    StorageReason storage_reason, std::optional<CopyPolicy> copy_policy) {
        // Setup this save
        auto                     t_h5     = tid::tic_scope("h5");
        auto                     t_reason = tid::tic_scope(enum2sv(storage_reason), tid::level::extra);
        StorageLevel             storage_level;
        std::string              state_prefix;
        std::string              model_prefix;
        std::string              timer_prefix;
        std::vector<std::string> table_prefxs;
        tools::finite::h5::save::setup_prefix(status, storage_reason, storage_level, state_name, state_prefix, model_prefix, timer_prefix, table_prefxs);

        std::string prefix;
        switch(storage_reason) {
            case StorageReason::MODEL: {
                prefix = model_prefix;
                break;
            }
            default: {
                prefix = state_prefix;
                break;
            }
        }

        save::data(h5file, data, data_name, prefix, storage_level, status);
        tools::common::h5::tmp::copy_from_tmp(status, h5file, storage_reason, copy_policy);
    }

    template void save::data(h5pp::File &h5file, const Eigen::Tensor<std::complex<double>, 2> &data, std::string_view data_name, std::string_view state_name,
                             const AlgorithmStatus &status, StorageReason storage_reason, std::optional<CopyPolicy> copy_policy);

    void save::setup_prefix(const AlgorithmStatus &status, const StorageReason &storage_reason, StorageLevel &storage_level, std::string_view state_name,
                            std::string &state_prefix, std::string &model_prefix, std::string &timer_prefix, std::vector<std::string> &table_prefxs) {
        // Setup this save
        state_prefix  = fmt::format("{}/{}", status.algo_type_sv(), state_name); // May get modified
        model_prefix  = fmt::format("{}/model", status.algo_type_sv());
        timer_prefix  = fmt::format("{}/timers", state_prefix);
        table_prefxs  = {fmt::format("{}/tables", state_prefix)}; // Common tables
        storage_level = StorageLevel::NONE;

        switch(storage_reason) {
            case StorageReason::NONE: return;
            case StorageReason::FINISHED: {
                storage_level = settings::storage::storage_level_finished;
                state_prefix += "/finished";
                break;
            }
            case StorageReason::SAVEPOINT: {
                storage_level = settings::storage::storage_level_savepoint;
                if(settings::storage::savepoint_frequency == 0) storage_level = StorageLevel::NONE;
                if(settings::storage::savepoint_frequency > 0 and num::mod(status.iter, settings::storage::savepoint_frequency) != 0)
                    storage_level = StorageLevel::NONE;

                state_prefix += "/savepoint";
                if(settings::storage::savepoint_keep_newest_only)
                    state_prefix += "/iter_last";
                else
                    state_prefix += fmt::format("/iter_{}", status.iter);
                break;
            }
            case StorageReason::CHECKPOINT: {
                storage_level = settings::storage::storage_level_checkpoint;
                if(settings::storage::checkpoint_frequency == 0) storage_level = StorageLevel::NONE;
                if(settings::storage::checkpoint_frequency > 0 and num::mod(status.iter, settings::storage::checkpoint_frequency) != 0)
                    storage_level = StorageLevel::NONE;

                state_prefix += "/checkpoint";
                if(settings::storage::checkpoint_keep_newest_only)
                    state_prefix += "/iter_last";
                else
                    state_prefix += fmt::format("/iter_{}", status.iter);
                break;
            }
            case StorageReason::BOND_INCREASE: {
                if(settings::strategy::bond_increase_when == UpdateWhen::NEVER) storage_level = StorageLevel::NONE;
                storage_level = settings::storage::storage_level_bond_state;
                state_prefix += "/bond";
                table_prefxs = {state_prefix}; // Does not pollute common tables
                break;
            }
            case StorageReason::TRNC_DECREASE: {
                if(settings::strategy::trnc_decrease_when == UpdateWhen::NEVER) storage_level = StorageLevel::NONE;
                storage_level = settings::storage::storage_level_trnc_state;
                state_prefix += "/trnc";
                table_prefxs = {state_prefix}; // Does not pollute common tables
                break;
            }
            case StorageReason::FES: {
                if(settings::strategy::fes_rate == 0.0) storage_level = StorageLevel::NONE;
                storage_level = settings::storage::storage_level_fes_state;
                state_prefix += "/fes";
                table_prefxs = {state_prefix}; // Does not pollute common tables
                break;
            }
            case StorageReason::PROJ_STATE: {
                storage_level = settings::storage::storage_level_proj_state;
                state_prefix += "/projection";
                table_prefxs = {state_prefix}; // Does not pollute common tables
                break;
            }
            case StorageReason::INIT_STATE: {
                storage_level = settings::storage::storage_level_init_state;
                state_prefix += "/state_init";
                table_prefxs = {state_prefix}; // Does not pollute common tables
                break;
            }
            case StorageReason::EMIN_STATE: {
                storage_level = settings::storage::storage_level_emin_state;
                break;
            }
            case StorageReason::EMAX_STATE: {
                storage_level = settings::storage::storage_level_emax_state;
                break;
            }
            case StorageReason::MODEL: {
                storage_level = settings::storage::storage_level_model;
                break;
            }
        }
    }

    void save::simulation(h5pp::File &h5file, const TensorsFinite &tensors, const AlgorithmStatus &status, const StorageReason &storage_reason,
                          std::optional<CopyPolicy> copy_policy) {
        save::simulation(h5file, *tensors.state, *tensors.model, *tensors.edges, status, storage_reason, copy_policy);
    }

    void save::simulation(h5pp::File &h5file, const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges, const AlgorithmStatus &status,
                          const StorageReason &storage_reason, std::optional<CopyPolicy> copy_policy) {
        auto                     t_h5     = tid::tic_scope("h5");
        auto                     t_sim    = tid::tic_scope("sim");
        auto                     t_reason = tid::tic_scope(enum2sv(storage_reason));
        StorageLevel             storage_level;
        std::string              state_prefix;
        std::string              model_prefix;
        std::string              timer_prefix;
        std::vector<std::string> table_prefxs;
        tools::finite::h5::save::setup_prefix(status, storage_reason, storage_level, state.get_name(), state_prefix, model_prefix, timer_prefix, table_prefxs);

        // Additional constraints
        switch(storage_reason) {
            case StorageReason::NONE: return;
            case StorageReason::SAVEPOINT:
            case StorageReason::CHECKPOINT: {
                if(not state.position_is_inward_edge()) storage_level = StorageLevel::NONE;
                break;
            }
            case StorageReason::PROJ_STATE: {
                if(storage_level != StorageLevel::NONE) {
                    if(not state.position_is_inward_edge()) storage_level = StorageLevel::NONE;
                    auto state_projected = tools::finite::ops::get_projection_to_nearest_sector(state, settings::strategy::target_sector,
                                                                                                svd::config(status.bond_lim, status.trnc_lim));
                    return save::simulation(h5file, state_projected, model, edges, status, storage_reason, copy_policy);
                }
                break;
            }
            case StorageReason::FINISHED: {
                // Either checkpoint or savepoint may already have the same data.
                // We can look for either of those under /common/finished, and make a soft-link
                // If /finished has a higher storage level, it will write what is missing into the soft link
                if(not h5file.linkExists(state_prefix)) {
                    auto duplicate_prefix = common::h5::find::find_duplicate_save(h5file, state_prefix, status);
                    if(duplicate_prefix) {
                        // We have a match! attrName is the path to a group containing finished data already.
                        // We make a soft link to it, and keep adding data if necessary.
                        tools::log->debug("Creating soft link {} -> {}", state_prefix, duplicate_prefix.value());
                        h5file.createSoftLink(duplicate_prefix.value(), state_prefix);
                    }
                }
                break;
            }

            case StorageReason::INIT_STATE:
            case StorageReason::EMIN_STATE:
            case StorageReason::EMAX_STATE:
            case StorageReason::BOND_INCREASE:
            case StorageReason::TRNC_DECREASE:
            case StorageReason::FES:
            case StorageReason::MODEL: break;
        }

        if(storage_level == StorageLevel::NONE) return;
        if(state_prefix.empty()) throw except::runtime_error("State prefix is empty");
        tools::log->debug("Writing to file: Reason [{}] | Level [{}] | state prefix [{}] | model prefix [{}]", enum2sv(storage_reason), enum2sv(storage_level),
                          state_prefix, model_prefix);

        // The file cam be kept open during writes
        h5file.setKeepFileOpened();

        // Start saving tensors and metadata
        if(storage_reason == StorageReason::MODEL) {
            tools::finite::h5::save::model(h5file, model_prefix, storage_level, model);
            tools::finite::h5::save::mpo(h5file, model_prefix, storage_level, model);
        } else {
            switch(storage_reason) {
                case StorageReason::BOND_INCREASE:
                case StorageReason::TRNC_DECREASE:
                case StorageReason::FES: break;
                default: {
                    tools::finite::h5::save::state(h5file, state_prefix, storage_level, state, status);
                    tools::finite::h5::save::correlations(h5file, state_prefix, storage_level, state, status);
                    break;
                }
            }
        }
        tools::common::h5::save::meta(h5file, storage_level, storage_reason, settings::model::model_type, settings::model::model_size, state.get_name(),
                                      state_prefix, model_prefix, table_prefxs, status);

        // The main results have now been written. Next we append data to tables
        for(const auto &table_prefix : table_prefxs) {
            if(storage_reason == StorageReason::MODEL) break;
            tools::common::h5::save::status(h5file, table_prefix, storage_level, status);
            tools::common::h5::save::mem(h5file, table_prefix, storage_level, status);
            tools::common::h5::save::timer(h5file, table_prefix, storage_level, status);
            tools::finite::h5::save::measurements(h5file, table_prefix, storage_level, state, model, edges, status);
            tools::finite::h5::save::bond_dimensions(h5file, table_prefix, storage_level, state, status);
            tools::finite::h5::save::schmidt_values(h5file, table_prefix, storage_level, state, status);
            tools::finite::h5::save::truncation_errors(h5file, table_prefix, storage_level, state, status);
            tools::finite::h5::save::entropies_neumann(h5file, table_prefix, storage_level, state, status);
            tools::finite::h5::save::entropies_renyi(h5file, table_prefix, storage_level, state, status);
            tools::finite::h5::save::entropies_number(h5file, table_prefix, storage_level, state, status);
            tools::finite::h5::save::expectations(h5file, table_prefix, storage_level, state, status);
            tools::finite::h5::save::structure_factors(h5file, table_prefix, storage_level, state, status);
            tools::finite::h5::save::kvornings_marker(h5file, table_prefix, storage_level, state, status);
            tools::finite::h5::save::number_probabilities(h5file, table_prefix, storage_level, state, status);
        }
        h5file.setKeepFileClosed();

        // Copy from temporary location to destination depending on given policy
        tools::common::h5::tmp::copy_from_tmp(status, h5file, storage_reason, copy_policy);
    }
}