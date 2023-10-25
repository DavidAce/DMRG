#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "io/hdf5_types.h"
#include "math/float.h"
#include "math/num.h"
#include "qm/spin.h"
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
        auto t_hdf = tid::tic_scope("measurements", tid::level::highest);

        sinfo.assert_well_defined();
        // Define the table
        std::string table_path = fmt::format("{}/measurements", sinfo.get_state_prefix());

        // Check if the current entry has already been appended
        auto attrs = tools::common::h5::save::get_save_attrs(h5file, table_path);
        if(not attrs.link_exists) h5file.createTable(h5pp_table_measurements_finite::get_h5t(), table_path, "measurements");
        if(attrs == sinfo) return;
        auto offset = tools::common::h5::save::get_table_offset(h5file, table_path, sinfo, attrs);

        // Define the table entry
        tools::log->trace("Appending to table: {}", table_path);
        h5pp_table_measurements_finite::table measurement_entry{};
        measurement_entry.iter     = static_cast<uint64_t>(sinfo.iter);
        measurement_entry.step     = static_cast<uint64_t>(sinfo.step);
        measurement_entry.position = static_cast<long>(sinfo.position);
        measurement_entry.event    = sinfo.storage_event;
        measurement_entry.length   = static_cast<uint64_t>(tools::finite::measure::length(state));
        if(status.algo_type != AlgorithmType::fLBIT) {
            measurement_entry.energy                 = tools::finite::measure::energy(state, model, edges);
            measurement_entry.energy_variance        = tools::finite::measure::energy_variance(state, model, edges);
            measurement_entry.energy_variance_lowest = status.energy_variance_lowest;
        }
        measurement_entry.norm             = tools::finite::measure::norm(state);
        measurement_entry.truncation_error = state.get_truncation_error_midchain();
        measurement_entry.bond_lim         = sinfo.bond_lim;
        measurement_entry.bond_max         = sinfo.bond_max;
        measurement_entry.bond_mid         = static_cast<long>(tools::finite::measure::bond_dimension_midchain(state));

        measurement_entry.entanglement_entropy = tools::finite::measure::entanglement_entropy_midchain(state);
        measurement_entry.renyi_entropy_2      = tools::finite::measure::renyi_entropy_midchain(state, 2);
        measurement_entry.renyi_entropy_3      = tools::finite::measure::renyi_entropy_midchain(state, 3);
        measurement_entry.renyi_entropy_4      = tools::finite::measure::renyi_entropy_midchain(state, 4);
        measurement_entry.renyi_entropy_inf    = tools::finite::measure::renyi_entropy_midchain(state, std::numeric_limits<double>::infinity());
        if(status.algo_type == AlgorithmType::fLBIT) { measurement_entry.number_entropy = tools::finite::measure::number_entropy_midchain(state); }
        measurement_entry.spin_global = tools::finite::measure::spin_components(state);
        measurement_entry.spin_local  = tools::finite::measure::expectation_value_xyz(state);
        if(status.algo_type == AlgorithmType::xDMRG and sinfo.storage_event == StorageEvent::LAST_STATE) {
            measurement_entry.structure_factors = tools::finite::measure::structure_factor_xyz(state);
        }
        measurement_entry.total_time     = status.wall_time;
        measurement_entry.algorithm_time = status.algo_time;
        measurement_entry.physical_time  = status.phys_time;

        tools::log->trace("Writing to table: {} | offset {}", table_path, offset);
        h5file.writeTableRecords(measurement_entry, table_path, offset);
        tools::common::h5::save::set_save_attrs(h5file, table_path, sinfo);
    }

    void save::correlations(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.algo_type != AlgorithmType::xDMRG) return;
        if(sinfo.storage_level <= StorageLevel::LIGHT) return;
        if(sinfo.storage_event != StorageEvent::LAST_STATE) return;
        auto correlation_matrix_xyz = tools::finite::measure::correlation_matrix_xyz(state);
        /* clang-format off */
        save::data(h5file, sinfo, correlation_matrix_xyz[0], "correlation_matrix_sx", sinfo.get_state_prefix());
        save::data(h5file, sinfo, correlation_matrix_xyz[1], "correlation_matrix_sy", sinfo.get_state_prefix());
        save::data(h5file, sinfo, correlation_matrix_xyz[2], "correlation_matrix_sz", sinfo.get_state_prefix());
        /* clang-format on */
    }

    void save::expectations(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.algo_type == AlgorithmType::xDMRG) {
            if(settings::storage::storage_level_tables <= StorageLevel::LIGHT) return; // Midchain is included in measurements table
            if(sinfo.storage_level <= StorageLevel::LIGHT) return;
            if(sinfo.storage_event != StorageEvent::LAST_STATE) return;
            auto t_hdf = tid::tic_scope("spin_local", tid::level::higher);
            tools::log->trace("Saving spin expectation values to {}", sinfo.get_state_prefix());
            auto expectation_values_xyz = tools::finite::measure::expectation_values_xyz(state);
            save::data_as_table(h5file, sinfo, expectation_values_xyz[0], "expectation_values_sx", "<sigma x>", "L_");
            save::data_as_table(h5file, sinfo, expectation_values_xyz[1], "expectation_values_sy", "<sigma y>", "L_");
            save::data_as_table(h5file, sinfo, expectation_values_xyz[2], "expectation_values_sz", "<sigma z>", "L_");
        }
        if(sinfo.algo_type == AlgorithmType::fLBIT) {
            if(settings::storage::storage_level_tables == StorageLevel::NONE) return; // Midchain is included in measurements table
            if(sinfo.storage_level == StorageLevel::NONE) return;
            if(sinfo.storage_event != StorageEvent::ITER_STATE) return;
            auto t_hdf = tid::tic_scope("spin_local", tid::level::higher);
            //            tools::log->trace("Saving spin expectation values to {}", sinfo.get_state_prefix());
            if(not state.measurements.expectation_values_sz.has_value())
                state.measurements.expectation_values_sz = measure::expectation_values(state, qm::spin::half::sz);
            //            save::data_as_table(h5file, sinfo, state.measurements.expectation_values_sz, "expectation_values_sz", "<sigma z>", "L_");
            auto table_path = fmt::format("{}/{}", sinfo.get_state_prefix(), "expectation_values_sz");
            tools::log->trace("Appending to table: {}", table_path);
            // Check if the current entry has already been appended
            auto attrs = tools::common::h5::save::get_save_attrs(h5file, table_path);
            if(attrs == sinfo) return;
            if(not attrs.link_exists) {
                auto                 rows = static_cast<hsize_t>(state.measurements.expectation_values_sz->dimension(0));
                std::vector<hsize_t> dims = {rows, 0};
                std::vector<hsize_t> chnk = {rows, 50};
                h5file.createDataset(table_path, h5pp::type::getH5Type<double>(), H5D_CHUNKED, dims, chnk);
            }
            h5file.appendToDataset(state.measurements.expectation_values_sz.value(), table_path, 1);
            tools::common::h5::save::set_save_attrs(h5file, table_path, sinfo);
        }
    }

    void save::kvornings_marker(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.algo_type != AlgorithmType::xDMRG) return;
        if(sinfo.storage_level <= StorageLevel::LIGHT) return;
        if(sinfo.storage_event != StorageEvent::LAST_STATE) return;
        auto kvornings_marker = tools::finite::measure::kvornings_marker(state);
        save::data_as_table(h5file, sinfo, kvornings_marker, "kvornings_marker", "Kvornings marker", "eigval");
    }

    template<typename T>
    void save::data_as_table(h5pp::File &h5file, const StorageInfo &sinfo, const T *const data, size_t size, std::string_view table_name,
                             std::string_view table_title, std::string_view fieldname) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        auto t_hdf = tid::tic_scope(table_name, tid::level::highest);
        sinfo.assert_well_defined();

        // Define the table
        auto table_path = fmt::format("{}/{}", sinfo.get_state_prefix(), table_name);
        auto h5_type    = h5pp_table_data<T>::register_table_type(h5pp::type::getH5Type<T>(), size, fieldname);

        // Check if the current entry has already been appended
        auto attrs = tools::common::h5::save::get_save_attrs(h5file, table_path);
        if(not attrs.link_exists) h5file.createTable(h5_type, table_path, table_title);
        if(attrs == sinfo) return;
        auto offset = tools::common::h5::save::get_table_offset(h5file, table_path, sinfo, attrs);

        tools::log->trace("Writing to table: {} | offset {}", table_path, offset);
        // Copy the data into an std::vector<std::byte> stream, which will act as a struct for our table entry
        auto entry = h5pp_table_data<T>::make_entry(sinfo.iter, sinfo.step, sinfo.position, sinfo.storage_event, sinfo.bond_lim, data, size);
        h5file.writeTableRecords(entry, table_path, offset);
        tools::common::h5::save::set_save_attrs(h5file, table_path, sinfo);
    }

    template<typename T>
    void save::data_as_table_vla(h5pp::File &h5file, const StorageInfo &sinfo, const std::vector<T> &data, const h5pp::hid::h5t &h5elem_t,
                                 std::string_view table_name, std::string_view table_title, std::string_view fieldname) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        auto t_hdf = tid::tic_scope(table_name);
        sinfo.assert_well_defined();

        // Define the table
        auto table_path = fmt::format("{}/{}", sinfo.get_state_prefix(), table_name);
        auto h5_type    = h5pp_table_data<T>::register_table_type(h5elem_t, data.size(), fieldname);

        // Check if the current entry has already been appended
        auto attrs = tools::common::h5::save::get_save_attrs(h5file, table_path);
        if(not attrs.link_exists) h5file.createTable(h5_type, table_path, table_title);
        if(attrs == sinfo) return;
        auto offset = tools::common::h5::save::get_table_offset(h5file, table_path, sinfo, attrs);

        tools::log->trace("Appending to table: {}", table_path);

        // Copy the data into an std::vector<std::byte> stream, which will act as a struct for our table entry
        tools::log->trace("Writing to table: {} | offset {}", table_path, offset);
        auto entry = h5pp_table_data<T>::make_entry(sinfo.iter, sinfo.step, sinfo.position, sinfo.storage_event, sinfo.bond_lim, data.data(), data.size());
        h5file.writeTableRecords(entry, table_path, offset);
        tools::common::h5::save::set_save_attrs(h5file, table_path, sinfo);
    }

    void save::bond_dimensions(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        if(settings::storage::storage_level_tables <= StorageLevel::LIGHT) return; // Midchain is included in measurements table
        data_as_table(h5file, sinfo, tools::finite::measure::bond_dimensions(state), "bond_dims", "Bond Dimensions", "L_");
    }

    void save::truncation_errors(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        if(settings::storage::storage_level_tables <= StorageLevel::LIGHT) return; // Midchain is included in measurements table
        auto t_hdf = tid::tic_scope("truncation_errors", tid::level::higher);
        data_as_table(h5file, sinfo, tools::finite::measure::truncation_errors(state), "truncation_errors", "Truncation errors", "L_");
    }

    void save::entropies_neumann(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        if(settings::storage::storage_level_tables <= StorageLevel::LIGHT) return; // Midchain is included in measurements table
        auto t_hdf = tid::tic_scope("entropies", tid::level::higher);
        data_as_table(h5file, sinfo, tools::finite::measure::entanglement_entropies(state), "entanglement_entropies", "Entanglement Entropies", "L_");
    }

    void save::number_probabilities(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.algo_type != AlgorithmType::fLBIT) return;
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event != StorageEvent::ITER_STATE) return;
        if(not state.measurements.number_probabilities) return;
        auto t_hdf      = tid::tic_scope("number_probabilities", tid::level::higher);
        auto table_path = fmt::format("{}/{}", sinfo.get_state_prefix(), "number_probabilities");
        tools::log->trace("Appending to table: {}", table_path);
        // Check if the current entry has already been appended
        auto attrs = tools::common::h5::save::get_save_attrs(h5file, table_path);
        if(attrs == sinfo) return;
        if(not attrs.link_exists) {
            auto                 rows = static_cast<hsize_t>(state.measurements.number_probabilities->dimension(0));
            auto                 cols = static_cast<hsize_t>(state.measurements.number_probabilities->dimension(1));
            std::vector<hsize_t> dims = {rows, cols, 0};
            std::vector<hsize_t> chnk = {rows, cols, 50};
            h5file.createDataset(table_path, h5pp::type::getH5Type<double>(), H5D_CHUNKED, dims, chnk, std::nullopt, 2);
            h5file.writeAttribute("n_count, site, time", table_path,  "index");
            h5file.writeAttribute("Probability of finding n_count particles to the left of a site at a time index", table_path,  "description");
        }
        // Do not append if the iteration number is smaller than the dataset iter dimension
        auto num_entries = h5file.getDatasetDimensions(table_path).back();
        if(sinfo.iter >= num_entries * settings::storage::storage_interval) {
            h5file.appendToDataset(state.measurements.number_probabilities.value(), table_path, 2);
            tools::common::h5::save::set_save_attrs(h5file, table_path, sinfo);
        }
    }

    void save::entropies_renyi(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level <= StorageLevel::LIGHT) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        if(settings::storage::storage_level_tables <= StorageLevel::LIGHT) return; // Midchain is included in measurements table
        auto t_hdf = tid::tic_scope("entropies", tid::level::higher);
        auto inf   = std::numeric_limits<double>::infinity();
        data_as_table(h5file, sinfo, tools::finite::measure::renyi_entropies(state, 2), "renyi_entropies_2", "Midchain Renyi Entropy 2", "L_");
        data_as_table(h5file, sinfo, tools::finite::measure::renyi_entropies(state, 3), "renyi_entropies_3", "Midchain Renyi Entropy 3", "L_");
        data_as_table(h5file, sinfo, tools::finite::measure::renyi_entropies(state, 4), "renyi_entropies_4", "Midchain Renyi Entropy 4", "L_");
        data_as_table(h5file, sinfo, tools::finite::measure::renyi_entropies(state, inf), "renyi_entropies_inf", "Midchain Renyi Entropy inf", "L_");
    }

    void save::entropies_number(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level <= StorageLevel::LIGHT) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        if(sinfo.algo_type != AlgorithmType::fLBIT) return;
        if(settings::storage::storage_level_tables <= StorageLevel::LIGHT) return; // Midchain is included in measurements table
        auto t_hdf = tid::tic_scope("entropies", tid::level::higher);
        data_as_table(h5file, sinfo, tools::finite::measure::number_entropies(state), "number_entropies", "Number entropies", "L_");
    }

    void save::bonds(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level == StorageLevel::NONE) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        if(settings::storage::storage_level_tables <= StorageLevel::LIGHT) return; // Midchain is included in measurements table
        auto t_hdf = tid::tic_scope("bonds", tid::level::higher);
        // Transform from cplx to real to save space
        auto        h5real      = h5pp::type::getH5Type<real>();
        auto        bonds_real  = std::vector<h5pp::varr_t<real>>();
        std::string column_name = "L_";
        std::string table_title = "Bonds (singular values)";
        bonds_real.reserve(state.get_length<size_t>() + 1);
        for(const auto &mps : state.mps_sites) {
            bonds_real.emplace_back(Eigen::Tensor<real, 1>(mps->get_L().real()));
            if(mps->isCenter()) { bonds_real.emplace_back(Eigen::Tensor<real, 1>(mps->get_LC().real())); }
        }
        data_as_table_vla(h5file, sinfo, bonds_real, h5real, "bonds", table_title, column_name);
    }

    void save::state(h5pp::File &h5file, const StorageInfo &sinfo, const StateFinite &state) {
        if(sinfo.storage_level < StorageLevel::FULL) return;
        if(sinfo.storage_event <= StorageEvent::MODEL) return;
        auto t_hdf      = tid::tic_scope("state", tid::level::higher);
        auto mps_prefix = sinfo.get_mps_prefix();

        // Check if the current entry has already been saved
        auto attrs = tools::common::h5::save::get_save_attrs(h5file, mps_prefix);
        if(attrs == sinfo) return;

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
        h5file.writeAttribute(sinfo.algo_type, mps_prefix, "algo_type", std::nullopt, h5_enum_algo_type::get_h5t());
        tools::common::h5::save::set_save_attrs(h5file, mps_prefix, sinfo);
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
        auto attrs = tools::common::h5::save::get_save_attrs(h5file, data_path);
        if(attrs == sinfo) return;

        H5D_layout_t layout;
        if(std::is_scalar_v<T>)
            layout = H5D_layout_t::H5D_COMPACT;
        else
            layout = H5D_layout_t::H5D_CHUNKED;
        h5file.writeDataset(data, data_path, layout);
        tools::common::h5::save::set_save_attrs(h5file, data_path, sinfo);
    }

    template void save::data(h5pp::File &h5file, const StorageInfo &sinfo, const Eigen::Tensor<double, 2> &data, std::string_view data_name,
                             std::string_view prefix);
    template void save::data(h5pp::File &h5file, const StorageInfo &sinfo, const Eigen::Tensor<double, 1> &data, std::string_view data_name,
                             std::string_view prefix);

    template<typename T>
    void save::data(h5pp::File &h5file, const StorageInfo &sinfo, const T &data, std::string_view data_name, CopyPolicy copy_policy) {
        // Setup this save
        auto t_h5    = tid::tic_scope("h5");
        auto t_event = tid::tic_scope(enum2sv(sinfo.storage_event), tid::level::highest);

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

        auto t_h5    = tid::tic_scope("h5");
        auto t_event = tid::tic_scope(enum2sv(status.event), tid::level::highest);

        tools::log->debug("Writing to file: Reason [{}] | Level [{}] | state prefix [{}]", enum2sv(sinfo.storage_event), enum2sv(sinfo.storage_level),
                          sinfo.get_state_prefix());

        // The file can be kept open during writes
        h5file.setKeepFileOpened();

        // The main results have now been written. Next we append data to tables
        tools::finite::h5::save::model(h5file, sinfo, model);
        tools::finite::h5::save::state(h5file, sinfo, state);
        tools::common::h5::save::initial_state_attrs(h5file, sinfo); // Save the initial state type and pattern (rather than the MPS itself)
        tools::common::h5::save::status(h5file, sinfo, status);
        tools::common::h5::save::mem(h5file, sinfo);
        tools::common::h5::save::timer(h5file, sinfo);
        tools::finite::h5::save::measurements(h5file, sinfo, state, model, edges, status);
        tools::finite::h5::save::bond_dimensions(h5file, sinfo, state);
        tools::finite::h5::save::bonds(h5file, sinfo, state);
        tools::finite::h5::save::truncation_errors(h5file, sinfo, state);
        tools::finite::h5::save::entropies_neumann(h5file, sinfo, state);
        tools::finite::h5::save::entropies_renyi(h5file, sinfo, state);
        tools::finite::h5::save::entropies_number(h5file, sinfo, state);
        tools::finite::h5::save::number_probabilities(h5file, sinfo, state);
        tools::finite::h5::save::expectations(h5file, sinfo, state);
        tools::finite::h5::save::correlations(h5file, sinfo, state);
        tools::finite::h5::save::kvornings_marker(h5file, sinfo, state);

        // The file can now be closed
        h5file.setKeepFileClosed();

        // Copy from temporary location to destination depending on given policy
        tools::common::h5::tmp::copy_from_tmp(h5file, sinfo.iter, sinfo.step, sinfo.storage_event, copy_policy);
    }
}