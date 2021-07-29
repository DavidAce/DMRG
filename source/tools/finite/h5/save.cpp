#include <algorithms/AlgorithmStatus.h>
#include <complex>
#include <config/settings.h>
#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <math/num.h>
#include <regex>
#include <tensors/model/ModelFinite.h>
#include <tensors/site/mpo/MpoSite.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateFinite.h>
#include <tensors/TensorsFinite.h>
#include <tid/tid.h>
#include <tools/common/h5.h>
#include <tools/common/log.h>
#include <tools/finite/h5.h>
#include <tools/finite/measure.h>
#include <tools/finite/ops.h>

namespace tools::finite::h5 {

    int save::decide_layout(std::string_view prefix_path) {
        return H5D_CHUNKED; // Let everything be chunked a while. When resuming, rewriting into savepoint/iter_? can lead datasets of different sizes
        std::string str(prefix_path);
        std::regex  rx(R"(savepoint/iter_[0-9])"); // Declare the regex with a raw string literal
        std::smatch m;
        if(regex_search(str, m, rx))
            return H5D_CONTIGUOUS;
        else
            return H5D_CHUNKED;
    }

    void save::bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile,
                                  const std::vector<std::string> &links) {
        if(save_log.empty()) {
            try {
                for(auto &link : links) {
                    if(h5ppFile.linkExists(link)) {
                        auto step      = h5ppFile.readAttribute<uint64_t>("step", link);
                        auto iter      = h5ppFile.readAttribute<uint64_t>("iteration", link);
                        save_log[link] = std::make_pair(iter, step);
                    }
                }
            } catch(const std::exception &ex) { tools::log->warn("Could not bootstrap save_log: {}", ex.what()); }
        }
    }

    void save::bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile,
                                  const std::string &link) {
        if(save_log.empty()) {
            try {
                if(h5ppFile.linkExists(link)) {
                    auto step      = h5ppFile.readAttribute<uint64_t>("step", link);
                    auto iter      = h5ppFile.readAttribute<uint64_t>("iteration", link);
                    save_log[link] = std::make_pair(iter, step);
                }
            } catch(const std::exception &ex) { tools::log->warn("Could not bootstrap save_log: {}", ex.what()); }
        }
    }

    void tools::finite::h5::save::measurements(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level,
                                               const TensorsFinite &tensors, const AlgorithmStatus &status, AlgorithmType algo_type) {
        save::measurements(h5ppFile, table_prefix, storage_level, *tensors.state, *tensors.model, *tensors.edges, status, algo_type);
    }

    void save::measurements(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level, const StateFinite &state,
                            const ModelFinite &model, const EdgesFinite &edges, const AlgorithmStatus &status, AlgorithmType algo_type) {
        if(storage_level == StorageLevel::NONE) return;
        // Check if the current entry has already been appended
        std::string                                                           table_path = fmt::format("{}/measurements", table_prefix);
        static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
        bootstrap_save_log(save_log, h5ppFile, table_path);
        auto save_point = std::make_pair(status.iter, status.step);
        if(save_log[table_path] == save_point) return;

        log->trace("Appending to table: {}", table_path);
        auto t_hdf = tid::tic_scope("measurements");
        h5pp_table_measurements_finite::register_table_type();
        if(not h5ppFile.linkExists(table_path)) h5ppFile.createTable(h5pp_table_measurements_finite::h5_type, table_path, "measurements");

        h5pp_table_measurements_finite::table measurement_entry{};
        measurement_entry.step                          = static_cast<uint64_t>(status.step);
        measurement_entry.iter                          = static_cast<uint64_t>(status.iter);
        measurement_entry.position                      = static_cast<long>(status.position);
        measurement_entry.length                        = static_cast<uint64_t>(tools::finite::measure::length(state));
        measurement_entry.bond_dimension_midchain       = static_cast<long>(tools::finite::measure::bond_dimension_midchain(state));
        measurement_entry.bond_dimension_current        = static_cast<long>(tools::finite::measure::bond_dimension_current(state));
        measurement_entry.bond_dimension_limit          = status.chi_lim;
        measurement_entry.bond_dimension_maximum        = status.chi_lim_max;
        measurement_entry.entanglement_entropy_midchain = tools::finite::measure::entanglement_entropy_midchain(state);
        measurement_entry.entanglement_entropy_current  = tools::finite::measure::entanglement_entropy_current(state);
        measurement_entry.number_entropy_midchain       = tools::finite::measure::number_entropy_midchain(state);
        measurement_entry.number_entropy_current        = tools::finite::measure::number_entropy_current(state);
        measurement_entry.norm                          = tools::finite::measure::norm(state);
        measurement_entry.energy                        = tools::finite::measure::energy(state, model, edges);
        measurement_entry.energy_per_site               = tools::finite::measure::energy_per_site(state, model, edges);
        if(algo_type != AlgorithmType::fLBIT) {
            measurement_entry.energy_variance                 = tools::finite::measure::energy_variance(state, model, edges);
            measurement_entry.energy_variance_per_site        = tools::finite::measure::energy_variance_per_site(state, model, edges);
            measurement_entry.energy_variance_lowest          = status.energy_variance_lowest;
            measurement_entry.energy_variance_per_site_lowest = status.energy_variance_lowest / state.get_length<double>();
        }
        measurement_entry.spin_components  = tools::finite::measure::spin_components(state);
        measurement_entry.truncation_error = state.get_truncation_error_midchain();
        measurement_entry.total_time       = status.wall_time;
        measurement_entry.algorithm_time   = status.algo_time;
        measurement_entry.physical_time    = status.phys_time;
        auto t_app                         = tid::tic_scope("append");
        h5ppFile.appendTableRecords(measurement_entry, table_path);
        h5ppFile.writeAttribute(status.iter, "iteration", table_path);
        h5ppFile.writeAttribute(status.step, "step", table_path);
        save_log[table_path] = save_point;
    }

    void save::state(h5pp::File &h5ppFile, const std::string &state_prefix, const StorageLevel &storage_level, const StateFinite &state,
                     const AlgorithmStatus &status) {
        if(storage_level == StorageLevel::NONE) return;
        auto t_hdf = tid::tic_scope("state");

        // Checks if the current entry has already been saved
        // If it is empty because we are resuming, check if there is a log entry on file already
        static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
        bootstrap_save_log(save_log, h5ppFile, {state_prefix + "/schmidt_midchain", state_prefix + "/mps"});

        auto save_point = std::make_pair(status.iter, status.step);

        auto layout = static_cast<H5D_layout_t>(decide_layout(state_prefix));

        std::string dsetName = state_prefix + "/schmidt_midchain";
        if(save_log[dsetName] != save_point) {
            /*! Writes down the center "Lambda" bond matrix (singular values). */
            tools::log->trace("Storing [{: ^6}]: mid bond matrix", enum2str(storage_level));
            h5ppFile.writeDataset(state.midchain_bond(), dsetName, layout);
            h5ppFile.writeAttribute(state.get_truncation_error_midchain(), "truncation_error", dsetName);
            h5ppFile.writeAttribute((state.get_length<long>() - 1) / 2, "position", dsetName);
            h5ppFile.writeAttribute(status.iter, "iteration", dsetName);
            h5ppFile.writeAttribute(status.step, "step", dsetName);
            h5ppFile.writeAttribute(status.chi_lim, "chi_lim", dsetName);
            h5ppFile.writeAttribute(status.chi_lim_max, "chi_lim_max", dsetName);
            save_log[dsetName] = save_point;
        }

        if(storage_level < StorageLevel::NORMAL) return;

        std::string mps_prefix = state_prefix + "/mps";
        if(save_log[mps_prefix] != save_point) {
            tools::log->trace("Storing [{: ^6}]: bond matrices", enum2str(storage_level));
            // There should be one more sites+1 number of L's, because there is also a center bond
            // However L_i always belongs M_i. Stick to this rule!
            // This means that some M_i has two bonds, one L_i to the left, and one L_C to the right.
            for(const auto &mps : state.mps_sites) {
                dsetName = fmt::format("{}/L_{}", mps_prefix, mps->get_position<long>());
                if(save_log[dsetName] == save_point) continue;
                h5ppFile.writeDataset(mps->get_L(), dsetName, layout);
                h5ppFile.writeAttribute(mps->get_position<long>(), "position", dsetName);
                h5ppFile.writeAttribute(mps->get_L().dimensions(), "dimensions", dsetName);
                h5ppFile.writeAttribute(mps->get_truncation_error(), "truncation_error", dsetName);
                if(mps->isCenter()) {
                    dsetName = mps_prefix + "/L_C";
                    h5ppFile.writeDataset(mps->get_LC(), dsetName, layout);
                    h5ppFile.writeAttribute(mps->get_position<long>(), "position", dsetName);
                    h5ppFile.writeAttribute(mps->get_LC().dimensions(), "dimensions", dsetName);
                    h5ppFile.writeAttribute(mps->get_truncation_error_LC(), "truncation_error", dsetName);
                }
                save_log[dsetName] = save_point;
            }
            h5ppFile.writeAttribute(state.get_length(), "model_size", mps_prefix);
            h5ppFile.writeAttribute(state.get_position<long>(), "position", mps_prefix);
            h5ppFile.writeAttribute(state.get_truncation_errors(), "truncation_errors", mps_prefix);
            h5ppFile.writeAttribute(state.get_labels(), "labels", mps_prefix);
            h5ppFile.writeAttribute(status.iter, "iteration", mps_prefix);
            h5ppFile.writeAttribute(status.step, "step", mps_prefix);
        }

        /*! Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation. */
        if(storage_level < StorageLevel::FULL) {
            save_log[mps_prefix] = save_point;
            return;
        }

        if(save_log[mps_prefix] != save_point) {
            tools::log->trace("Storing [{: ^6}]: mps tensors", enum2str(storage_level));
            for(const auto &mps : state.mps_sites) {
                dsetName = fmt::format("{}/M_{}", mps_prefix, mps->get_position<long>());
                if(save_log[dsetName] == save_point) continue;
                h5ppFile.writeDataset(mps->get_M_bare(), dsetName, layout); // Important to write bare matrices!!
                h5ppFile.writeAttribute(mps->get_position<long>(), "position", dsetName);
                h5ppFile.writeAttribute(mps->get_M_bare().dimensions(), "dimensions", dsetName);
                h5ppFile.writeAttribute(mps->get_label(), "label", dsetName);
                h5ppFile.writeAttribute(mps->get_unique_id(), "unique_id", dsetName);
                save_log[dsetName] = save_point;
            }
            save_log[mps_prefix] = save_point;
        }
    }

    /*! Write down the Hamiltonian model type and site info as attributes */
    void save::model(h5pp::File &h5ppFile, const std::string &table_prefix, const StorageLevel &storage_level, const ModelFinite &model) {
        if(storage_level == StorageLevel::NONE) return;
        std::string table_path = fmt::format("{}/hamiltonian", table_prefix);
        if(h5ppFile.linkExists(table_path)) return tools::log->debug("The hamiltonian has already been written to [{}]", table_path);

        tools::log->trace("Storing table: [{}]", table_path);
        auto t_hdf = tid::tic_scope("model");
        for(auto site = 0ul; site < model.get_length(); site++) model.get_mpo(site).save_hamiltonian(h5ppFile, table_path);
        h5ppFile.writeAttribute(enum2str(settings::model::model_type), "model_type", table_path);
        h5ppFile.writeAttribute(settings::model::model_size, "model_size", table_path);
    }

    /*! Write all the MPO's with site info in attributes */
    void tools::finite::h5::save::mpo(h5pp::File &h5ppFile, const std::string &model_prefix, const StorageLevel &storage_level, const ModelFinite &model) {
        if(storage_level < StorageLevel::FULL) return;
        std::string mpo_prefix = fmt::format("{}/mpo", model_prefix);
        // We do not expect the MPO's to change. Therefore if they exist, there is nothing else to do here
        if(h5ppFile.linkExists(mpo_prefix)) return tools::log->trace("The MPO's have already been written to [{}]", mpo_prefix);
        auto t_hdf = tid::tic_scope("mpo");
        tools::log->trace("Storing [{: ^6}]: mpos", enum2str(storage_level));
        for(size_t pos = 0; pos < model.get_length(); pos++) { model.get_mpo(pos).save_mpo(h5ppFile, mpo_prefix); }
        h5ppFile.writeAttribute(settings::model::model_size, "model_size", mpo_prefix);
        h5ppFile.writeAttribute(enum2str(settings::model::model_type), "model_type", mpo_prefix);
    }

    /*! Write down measurements that can't fit in a table */
    void save::entgm(h5pp::File &h5ppFile, const std::string &state_prefix, const StorageLevel &storage_level, const StateFinite &state,
                     const AlgorithmStatus &status) {
        if(storage_level < StorageLevel::NORMAL) return;
        state.do_all_measurements();
        auto t_hdf = tid::tic_scope("entgm");

        tools::log->trace("Storing [{: ^6}]: bond dimensions", enum2str(storage_level));
        h5ppFile.writeDataset(tools::finite::measure::bond_dimensions(state), state_prefix + "/bond_dimensions");
        h5ppFile.writeAttribute(status.chi_lim, "chi_lim", state_prefix + "/bond_dimensions");
        h5ppFile.writeAttribute(status.chi_lim_max, "chi_lim_max", state_prefix + "/bond_dimensions");

        tools::log->trace("Storing [{: ^6}]: entanglement entropies", enum2str(storage_level));
        h5ppFile.writeDataset(tools::finite::measure::entanglement_entropies(state), state_prefix + "/entanglement_entropies");
        h5ppFile.writeDataset(tools::finite::measure::renyi_entropies(state, 2), state_prefix + "/renyi_2");
        h5ppFile.writeDataset(tools::finite::measure::renyi_entropies(state, 3), state_prefix + "/renyi_3");
        h5ppFile.writeDataset(tools::finite::measure::renyi_entropies(state, 4), state_prefix + "/renyi_4");
        h5ppFile.writeDataset(tools::finite::measure::renyi_entropies(state, 100), state_prefix + "/renyi_100");
        h5ppFile.writeDataset(tools::finite::measure::truncation_errors(state), state_prefix + "/truncation_errors");
        if(not tools::finite::measure::number_entropies(state).empty())
            h5ppFile.writeDataset(tools::finite::measure::number_entropies(state), state_prefix + "/number_entropies");
    }

    template<typename T>
    void save::data(h5pp::File &h5ppFile, const T &data, const std::string &data_name, const std::string &state_name, const AlgorithmStatus &status,
                    StorageReason storage_reason, std::optional<CopyPolicy> copy_policy) {
        // Setup this save
        auto                     t_save = tid::tic_scope(fmt::format("h5.{}", enum2str(storage_reason)));
        StorageLevel             storage_level;
        std::string              state_prefix;
        std::string              model_prefix;
        std::vector<std::string> table_prefxs;
        tools::finite::h5::save::setup_prefix(status, storage_reason, storage_level, state_name, state_prefix, model_prefix, table_prefxs);

        std::string data_path;
        switch(storage_reason) {
            case StorageReason::MODEL: {
                data_path = fmt::format("{}/{}", model_prefix, data_name);
                break;
            }
            default: {
                data_path = fmt::format("{}/{}", state_prefix, data_name);
                break;
            }
        }
        // Checks if the current entry has already been saved
        // If it is empty because we are resuming, check if there is a log entry on file already
        static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
        bootstrap_save_log(save_log, h5ppFile, {data_path});
        auto save_point = std::make_pair(status.iter, status.step);
        if(save_log[data_path] != save_point or status.step == 0) {
            auto layout = static_cast<H5D_layout_t>(decide_layout(data_path));
            h5ppFile.writeDataset(data, data_path, layout);
            h5ppFile.writeAttribute(status.iter, "iteration", data_path);
            h5ppFile.writeAttribute(status.step, "step", data_path);
            save_log[data_path] = save_point;
            tools::common::h5::tmp::copy_from_tmp(status, h5ppFile, storage_reason, copy_policy);
        }
    }

    template void save::data(h5pp::File &h5ppFile, const Eigen::Tensor<std::complex<double>, 2> &data, const std::string &data_name,
                             const std::string &state_name, const AlgorithmStatus &status, StorageReason storage_reason, std::optional<CopyPolicy> copy_policy);

    void save::setup_prefix(const AlgorithmStatus &status, const StorageReason &storage_reason, StorageLevel &storage_level, const std::string &state_name,
                            std::string &state_prefix, std::string &model_prefix, std::vector<std::string> &table_prefxs) {
        // Setup this save
        state_prefix  = fmt::format("{}/{}", status.algo_type_sv(), state_name); // May get modified
        model_prefix  = fmt::format("{}/model", status.algo_type_sv());
        table_prefxs  = {fmt::format("{}/tables", state_prefix)}; // Common tables
        storage_level = StorageLevel::NONE;

        switch(storage_reason) {
            case StorageReason::FINISHED: {
                if(status.algorithm_has_succeeded)
                    storage_level = settings::output::storage_level_good_state;
                else
                    storage_level = settings::output::storage_level_fail_state;
                state_prefix += "/finished";
                table_prefxs.emplace_back(state_prefix); // Appends to its own table as well as the common ones
                break;
            }
            case StorageReason::SAVEPOINT: {
                storage_level = settings::output::storage_level_savepoint;
                if(status.algo_stop == AlgorithmStop::NONE) {
                    if(num::mod(status.iter, settings::output::savepoint_frequency) != 0) storage_level = StorageLevel::NONE;
                }
                state_prefix += "/savepoint";
                if(settings::output::savepoint_keep_newest_only)
                    state_prefix += "/iter_last";
                else
                    state_prefix += fmt::format("/iter_{}", status.iter);
                table_prefxs = {state_prefix}; // Does not pollute common tables
                break;
            }
            case StorageReason::CHECKPOINT: {
                storage_level = settings::output::storage_level_checkpoint;
                if(status.algo_stop == AlgorithmStop::NONE) {
                    if(num::mod(status.iter, settings::output::checkpoint_frequency) != 0) storage_level = StorageLevel::NONE;
                }
                state_prefix += "/checkpoint";
                if(settings::output::checkpoint_keep_newest_only)
                    state_prefix += "/iter_last";
                else
                    state_prefix += fmt::format("/iter_{}", status.iter);
                table_prefxs.emplace_back(state_prefix); // Appends to its own table as well as the common ones
                break;
            }
            case StorageReason::CHI_UPDATE: {
                storage_level = settings::output::storage_level_checkpoint;
                if(not settings::chi_lim_grow(status.algo_type)) storage_level = StorageLevel::NONE;
                // If we have updated chi we may want to write a projection too
                state_prefix += fmt::format("/checkpoint/chi_{}", status.chi_lim);
                table_prefxs = {state_prefix}; // Does not pollute common tables
                break;
            }
            case StorageReason::PROJ_STATE: {
                storage_level = settings::output::storage_level_proj_state;
                state_prefix += "/projection";
                table_prefxs = {state_prefix}; // Does not pollute common tables
                break;
            }
            case StorageReason::INIT_STATE: {
                storage_level = settings::output::storage_level_init_state;
                state_prefix += "/state_init";
                table_prefxs = {state_prefix}; // Does not pollute common tables
                break;
            }
            case StorageReason::EMIN_STATE: {
                storage_level = settings::output::storage_level_emin_state;
                break;
            }
            case StorageReason::EMAX_STATE: {
                storage_level = settings::output::storage_level_emax_state;
                break;
            }
            case StorageReason::MODEL: {
                storage_level = settings::output::storage_level_model;
                break;
            }
        }
    }

    void save::simulation(h5pp::File &h5ppfile, const TensorsFinite &tensors, const AlgorithmStatus &status, const StorageReason &storage_reason,
                          std::optional<CopyPolicy> copy_policy) {
        save::simulation(h5ppfile, *tensors.state, *tensors.model, *tensors.edges, status, storage_reason, copy_policy);
    }

    void save::simulation(h5pp::File &h5ppfile, const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges, const AlgorithmStatus &status,
                          const StorageReason &storage_reason, std::optional<CopyPolicy> copy_policy) {
        auto                     t_save = tid::tic_scope(fmt::format("h5.{}", enum2str(storage_reason)));
        StorageLevel             storage_level;
        std::string              state_prefix;
        std::string              model_prefix;
        std::vector<std::string> table_prefxs;
        tools::finite::h5::save::setup_prefix(status, storage_reason, storage_level, state.get_name(), state_prefix, model_prefix, table_prefxs);
        // Additional constraints
        switch(storage_reason) {
            case StorageReason::FINISHED: break;
            case StorageReason::SAVEPOINT:
            case StorageReason::CHECKPOINT: {
                if(not state.position_is_inward_edge()) storage_level = StorageLevel::NONE;
                break;
            }
            case StorageReason::CHI_UPDATE: {
                if(not settings::chi_lim_grow(status.algo_type)) storage_level = StorageLevel::NONE;
                break;
            }
            case StorageReason::PROJ_STATE: {
                auto abs_spin_component = std::abs(tools::finite::measure::spin_component(state, settings::strategy::target_sector));
                if(std::abs(abs_spin_component - 1.0) > 1e-6) {
                    auto state_projected =
                        tools::finite::ops::get_normalized_projection_to_nearest_sector(state, settings::strategy::target_sector, status.chi_lim);
                    return save::simulation(h5ppfile, state_projected, model, edges, status, storage_reason, copy_policy);
                }
                break;
            }
            case StorageReason::INIT_STATE:
            case StorageReason::EMIN_STATE:
            case StorageReason::EMAX_STATE:
            case StorageReason::MODEL: break;
        }

        if(storage_level == StorageLevel::NONE) return;
        if(state_prefix.empty()) throw std::runtime_error("State prefix is empty");
        tools::log->info("Writing to file: Reason [{}] | Level [{}] | state prefix [{}] | model prefix [{}]", enum2str(storage_reason), enum2str(storage_level),
                         state_prefix, model_prefix);

        // The file cam be kept open during writes
        h5ppfile.setKeepFileOpened();

        // Start saving tensors and metadata
        if(storage_reason == StorageReason::MODEL) {
            tools::finite::h5::save::model(h5ppfile, model_prefix, storage_level, model);
            tools::finite::h5::save::mpo(h5ppfile, model_prefix, storage_level, model);
        } else {
            tools::finite::h5::save::state(h5ppfile, state_prefix, storage_level, state, status);
            tools::finite::h5::save::entgm(h5ppfile, state_prefix, storage_level, state, status);
            tools::common::h5::save::timer(h5ppfile, state_prefix, storage_level, status);
        }

        tools::common::h5::save::meta(h5ppfile, storage_level, storage_reason, settings::model::model_type, settings::model::model_size, status.algo_type,
                                      state.get_name(), state_prefix, model_prefix, status);

        // The main results have now been written. Next we append data to tables
        for(const auto &table_prefix : table_prefxs) {
            if(storage_reason == StorageReason::MODEL) break;
            tools::common::h5::save::status(h5ppfile, table_prefix, storage_level, status);
            tools::common::h5::save::mem(h5ppfile, table_prefix, storage_level, status);
            tools::finite::h5::save::measurements(h5ppfile, table_prefix, storage_level, state, model, edges, status, status.algo_type);
        }
        h5ppfile.setKeepFileClosed();

        // Copy from temporary location to destination depending on given policy
        tools::common::h5::tmp::copy_from_tmp(status, h5ppfile, storage_reason, copy_policy);
    }
}