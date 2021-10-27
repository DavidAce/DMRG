#include "../h5.h"
#include <algorithms/AlgorithmStatus.h>
#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <regex>
#include <tensors/edges/EdgesInfinite.h>
#include <tensors/model/ModelInfinite.h>
#include <tensors/site/mpo/MpoSite.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateInfinite.h>
#include <tensors/TensorsInfinite.h>
#include <tid/tid.h>
#include <tools/infinite/measure.h>

namespace tools::infinite::h5::save {
    void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5file,
                            const std::vector<std::string> &links) {
        if(save_log.empty()) {
            try {
                for(auto &link : links) {
                    if(h5file.linkExists(link)) {
                        auto step      = h5file.readAttribute<uint64_t>("step", link);
                        auto iter      = h5file.readAttribute<uint64_t>("iteration", link);
                        save_log[link] = std::make_pair(iter, step);
                    }
                }
            } catch(const std::exception &ex) { tools::log->warn("Could not bootstrap save_log: {}", ex.what()); }
        }
    }
}

int tools::infinite::h5::save::decide_layout(std::string_view prefix_path) {
    std::string str(prefix_path);
    std::regex  rx(R"(point/iter_[0-9])"); // Declare the regex with a raw string literal
    std::smatch m;
    if(regex_search(str, m, rx))
        return H5D_CONTIGUOUS;
    else
        return H5D_CHUNKED;
}

void tools::infinite::h5::save::state(h5pp::File &h5file, std::string_view state_prefix, const StorageLevel &storage_level, const StateInfinite &state,
                                      const AlgorithmStatus &status) {
    if(storage_level == StorageLevel::NONE) return;

    // Checks if the current entry has already been saved
    // If it is empty because we are resuming, check if there is a log entry on file already
    auto                                                                  tic = tid::tic_token("state");
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    auto                                                                  dset_schmidt = fmt::format("{}/schmidt_midchain", state_prefix);
    auto                                                                  mps_prefix   = fmt::format("{}/mps", state_prefix);

    bootstrap_save_log(save_log, h5file, {dset_schmidt, mps_prefix});
    auto save_point = std::make_pair(status.iter, status.step);
    auto layout     = static_cast<H5D_layout_t>(decide_layout(state_prefix));

    if(save_log[dset_schmidt] != save_point) {
        tools::log->trace("Storing [{: ^6}]: mid bond matrix", enum2sv(storage_level));
        h5file.writeDataset(state.LC(), dset_schmidt, layout);
        h5file.writeAttribute(state.get_truncation_error(), "truncation_error", dset_schmidt);
        h5file.writeAttribute(status.chi_lim, "chi_lim", dset_schmidt);
        h5file.writeAttribute(status.chi_lim_max, "chi_lim_max", dset_schmidt);
        save_log[dset_schmidt] = save_point;
    }
    if(storage_level < StorageLevel::NORMAL) return;
    if(save_log[mps_prefix] != save_point) {
        tools::log->trace("Storing [{: ^6}]: bond matrices", enum2sv(storage_level));
        auto dsetName = mps_prefix + "/L_A";
        if(save_log[dsetName] != save_point) {
            h5file.writeDataset(state.LA(), dsetName);
            h5file.writeAttribute(state.LA().dimensions(), "dimensions", dsetName);
            h5file.writeAttribute(state.get_truncation_error(), "truncation_error", dsetName);
            h5file.writeAttribute(status.chi_lim, "chi_lim", dsetName);
            h5file.writeAttribute(status.chi_lim_max, "chi_lim_max", dsetName);
            save_log[dsetName] = save_point;
        }
        dsetName = mps_prefix + "/L_B";
        if(save_log[dsetName] != save_point) {
            h5file.writeDataset(state.LB(), dsetName);
            h5file.writeAttribute(state.LB().dimensions(), "dimensions", dsetName);
            h5file.writeAttribute(state.get_truncation_error(), "truncation_error", dsetName);
            h5file.writeAttribute(status.chi_lim, "chi_lim", dsetName);
            h5file.writeAttribute(status.chi_lim_max, "chi_lim_max", dsetName);
            save_log[dsetName] = save_point;
        }
        dsetName = mps_prefix + "/L_C";
        if(save_log[dsetName] != save_point) {
            h5file.writeDataset(state.LC(), dsetName);
            h5file.writeAttribute(state.LC().dimensions(), "dimensions", dsetName);
            h5file.writeAttribute(state.get_truncation_error(), "truncation_error", dsetName);
            h5file.writeAttribute(status.chi_lim, "chi_lim", dsetName);
            h5file.writeAttribute(status.chi_lim_max, "chi_lim_max", dsetName);
            save_log[dsetName] = save_point;
        }
        h5file.writeAttribute(status.iter, "iteration", mps_prefix);
        h5file.writeAttribute(status.step, "step", mps_prefix);
    }
    /*! Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation. */
    if(storage_level < StorageLevel::FULL) {
        save_log[mps_prefix] = save_point;
        return;
    }

    if(save_log[mps_prefix] != save_point) {
        h5file.writeDataset(state.A_bare(), mps_prefix + "/M_A");
        h5file.writeAttribute(state.A_bare().dimensions(), "dimensions", mps_prefix + "/M_A");
        h5file.writeDataset(state.B(), mps_prefix + "/M_B");
        h5file.writeAttribute(state.B().dimensions(), "dimensions", mps_prefix + "/M_B");
        save_log[mps_prefix] = save_point;
    }
}

void tools::infinite::h5::save::edges(h5pp::File &h5file, std::string_view edges_prefix, const StorageLevel &storage_level, const EdgesInfinite &edges) {
    if(storage_level < StorageLevel::NORMAL) return;
    auto        tic = tid::tic_token("edges");
    const auto &ene = edges.get_ene_blk();
    const auto &var = edges.get_var_blk();
    h5file.writeDataset(ene.L, fmt::format("{}/eneL", edges_prefix));
    h5file.writeDataset(ene.R, fmt::format("{}/eneR", edges_prefix));
    h5file.writeDataset(var.L, fmt::format("{}/varL", edges_prefix));
    h5file.writeDataset(var.R, fmt::format("{}/varR", edges_prefix));
}

/*! Write down the Hamiltonian model type and site info as attributes */
void tools::infinite::h5::save::model(h5pp::File &h5file, std::string_view model_prefix, const StorageLevel &storage_level, const ModelInfinite &model) {
    if(storage_level < StorageLevel::LIGHT) return;
    tools::log->trace("Writing Hamiltonian model");
    auto table_path = fmt::format("{}/hamiltonian", model_prefix);
    if(h5file.linkExists(table_path)) return tools::log->debug("The hamiltonian has already been written to [{}]", table_path);
    tools::log->trace("Storing table: [{}]", table_path);
    auto t_ham = tid::tic_token("save_hamiltonian");
    model.get_mpo_siteA().save_hamiltonian(h5file, table_path);
    model.get_mpo_siteB().save_hamiltonian(h5file, table_path);
    h5file.writeAttribute(enum2sv(settings::model::model_type), "model_type", table_path);
    h5file.writeAttribute(settings::model::model_size, "model_size", table_path);
}

void tools::infinite::h5::save::mpo(h5pp::File &h5file, std::string_view model_prefix, const StorageLevel &storage_level, const ModelInfinite &model) {
    if(storage_level < StorageLevel::FULL) return;
    // We do not expect the MPO's to change. Therefore if they exist, there is nothing else to do here
    if(h5file.linkExists(model_prefix)) return tools::log->trace("The model has already been written to [{}]", model_prefix);
    tools::log->trace("Storing [{: ^6}]: mpo tensors", enum2sv(storage_level));
    auto tic        = tid::tic_token("mpo");
    auto mpo_prefix = fmt::format("{}/mpo", model_prefix);
    model.get_mpo_siteA().save_mpo(h5file, mpo_prefix);
    model.get_mpo_siteB().save_mpo(h5file, mpo_prefix);

    h5file.writeAttribute(2, "model_size", mpo_prefix);
    h5file.writeAttribute(enum2sv(settings::model::model_type), "model_type", mpo_prefix);

    h5file.writeAttribute(2, "model_size", model_prefix);
    h5file.writeAttribute(enum2sv(settings::model::model_type), "model_type", model_prefix);
}

void tools::infinite::h5::save::measurements(h5pp::File &h5file, std::string_view table_prefix, const StorageLevel &storage_level,
                                             const TensorsInfinite &tensors, const AlgorithmStatus &status) {
    if(storage_level == StorageLevel::NONE) return;
    auto table_path = fmt::format("{}/measurements", table_prefix);
    // Check if the current entry has already been appended
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    auto                                                                  save_point = std::make_pair(status.iter, status.step);
    if(save_log[table_path] == save_point) return;

    log->trace("Appending to table: {}", table_path);
    h5pp_table_measurements_infinite::register_table_type();
    if(not h5file.linkExists(table_path)) h5file.createTable(h5pp_table_measurements_infinite::h5_type, table_path, "measurements");

    h5pp_table_measurements_infinite::table measurement_entry{};
    const auto                             &state = *tensors.state;

    measurement_entry.step                   = static_cast<uint64_t>(status.step);
    measurement_entry.iter                   = static_cast<uint64_t>(status.iter);
    measurement_entry.position               = static_cast<int64_t>(status.position);
    measurement_entry.length                 = static_cast<uint64_t>(tools::infinite::measure::length(tensors));
    measurement_entry.bond_dimension         = tools::infinite::measure::bond_dimension(state);
    measurement_entry.bond_dimension_limit   = status.chi_lim;
    measurement_entry.bond_dimension_maximum = status.chi_lim_max;
    measurement_entry.entanglement_entropy   = tools::infinite::measure::entanglement_entropy(state);
    measurement_entry.norm                   = tools::infinite::measure::norm(*tensors.state);
    if(std::abs(status.delta_t) == 0) {
        // MPO calculations do not make sense for iTEBD, and we know delta_t != 0 on iTEBD.
        measurement_entry.energy_mpo                   = tools::infinite::measure::energy_mpo(tensors);
        measurement_entry.energy_per_site_mpo          = tools::infinite::measure::energy_per_site_mpo(tensors);
        measurement_entry.energy_variance_mpo          = tools::infinite::measure::energy_variance_mpo(tensors);
        measurement_entry.energy_variance_per_site_mpo = tools::infinite::measure::energy_variance_per_site_mpo(tensors);
    }
    measurement_entry.energy_per_site_ham          = tools::infinite::measure::energy_per_site_ham(tensors);
    measurement_entry.energy_per_site_mom          = tools::infinite::measure::energy_per_site_mom(tensors);
    measurement_entry.energy_variance_per_site_ham = tools::infinite::measure::energy_variance_per_site_ham(tensors);
    measurement_entry.energy_variance_per_site_mom = tools::infinite::measure::energy_variance_per_site_mom(tensors);
    measurement_entry.truncation_error             = tools::infinite::measure::truncation_error(state);
    measurement_entry.wall_time                    = status.wall_time;
    measurement_entry.phys_time                    = status.phys_time;
    measurement_entry.time_step                    = status.delta_t;
    auto t_app                                     = tid::tic_scope("append_table");
    h5file.appendTableRecords(measurement_entry, table_path);
    save_log[table_path] = save_point;
}
