//
// Created by david on 2019-03-09.
//

#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <h5pp/h5pp.h>
#include <regex>
#include <tensors/model/class_model_finite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/io.h>
#include <tools/finite/measure.h>
using Scalar = std::complex<double>;

namespace tools::finite::io::h5dset{
    void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> & save_log, const h5pp::File &h5ppFile, const std::vector<std::string> & links){
        if(save_log.empty()){
            try {
                for(auto &link : links) {
                    if(h5ppFile.linkExists(link)) {
                        auto step      = h5ppFile.readAttribute<uint64_t>("step", link);
                        auto iter      = h5ppFile.readAttribute<uint64_t>("iteration", link);
                        save_log[link] = std::make_pair(iter, step);
                    }
                }
            }catch(const std::exception & ex){
                tools::log->warn("Could not bootstrap save_log: {}", ex.what());
            }
        }
    }
}



int tools::finite::io::h5dset::decide_layout(std::string_view prefix_path) {
    return H5D_CHUNKED; // Let everything be chunked a while. When resuming, rewriting into checkpoint/iter_? can lead datasets of different sizes
    std::string str(prefix_path);
    std::regex  rx(R"(checkpoint/iter_[0-9])"); // Declare the regex with a raw string literal
    std::smatch m;
    if(regex_search(str, m, rx)) return H5D_CONTIGUOUS;
    else
        return H5D_CHUNKED;
}

void tools::finite::io::h5dset::save_state(h5pp::File &h5ppFile, const std::string &state_prefix, const StorageLevel &storage_level,
                                           const class_state_finite &state, const class_algorithm_status &status) {
    if(storage_level == StorageLevel::NONE) return;

    // Checks if the current entry has already been saved
    // If it is empty because we are resuming, check if there is a log entry on file already
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    bootstrap_save_log(save_log,h5ppFile,{state_prefix + "/schmidt_midchain", state_prefix + "/mps"});

    auto save_point = std::make_pair(status.iter, status.step);

    auto layout = static_cast<H5D_layout_t>(decide_layout(state_prefix));

    std::string dsetName = state_prefix + "/schmidt_midchain";
    if(save_log[dsetName] != save_point) {
        /*! Writes down the center "Lambda" bond matrix (singular values). */
        tools::log->trace("Storing [{: ^6}]: mid bond matrix", enum2str(storage_level));
        tools::common::profile::get_default_prof()["t_hdf"]->tic();
        h5ppFile.writeDataset(state.midchain_bond(), dsetName, layout);
        h5ppFile.writeAttribute(state.get_truncation_error_midchain(), "truncation_error", dsetName);
        h5ppFile.writeAttribute((state.get_length() - 1) / 2, "position", dsetName);
        h5ppFile.writeAttribute(status.iter, "iteration", dsetName);
        h5ppFile.writeAttribute(status.step, "step", dsetName);
        h5ppFile.writeAttribute(status.chi_lim, "chi_lim", dsetName);
        h5ppFile.writeAttribute(status.chi_lim_max, "chi_lim_max", dsetName);
        tools::common::profile::get_default_prof()["t_hdf"]->toc();
        save_log[dsetName] = save_point;
    }

    if(storage_level < StorageLevel::NORMAL) return;

    std::string mps_prefix = state_prefix + "/mps";
    if(save_log[mps_prefix] != save_point) {
        tools::log->trace("Storing [{: ^6}]: bond matrices", enum2str(storage_level));
        tools::common::profile::get_default_prof()["t_hdf"]->tic();
        // There should be one more sites+1 number of L's, because there is also a center bond
        // However L_i always belongs M_i. Stick to this rule!
        // This means that some M_i has two bonds, one L_i to the left, and one L_C to the right.
        for(size_t pos = 0; pos < state.get_length(); pos++) {
            dsetName = fmt::format("{}/L_{}", mps_prefix, pos);
            if(save_log[dsetName] == save_point) continue;
            h5ppFile.writeDataset(state.get_mps_site(pos).get_L(), dsetName, layout);
            h5ppFile.writeAttribute(pos, "position", dsetName);
            h5ppFile.writeAttribute(state.get_mps_site(pos).get_L().dimensions(), "dimensions", dsetName);
            h5ppFile.writeAttribute(state.get_mps_site(pos).get_truncation_error(), "truncation_error", dsetName);
            if(state.get_mps_site(pos).isCenter()) {
                dsetName = mps_prefix + "/L_C";
                h5ppFile.writeDataset(state.get_mps_site(pos).get_LC(), dsetName, layout);
                h5ppFile.writeAttribute(pos, "position", dsetName);
                h5ppFile.writeAttribute(state.get_mps_site(pos).get_LC().dimensions(), "dimensions", dsetName);
                h5ppFile.writeAttribute(state.get_mps_site(pos).get_truncation_error_LC(), "truncation_error", dsetName);
            }
            save_log[dsetName] = save_point;
        }
        h5ppFile.writeAttribute(state.get_length(), "model_size", mps_prefix);
        h5ppFile.writeAttribute(state.get_position(), "position", mps_prefix);
        h5ppFile.writeAttribute(state.get_truncation_errors(), "truncation_errors", mps_prefix);
        h5ppFile.writeAttribute(state.get_labels(), "labels", mps_prefix);
        h5ppFile.writeAttribute(status.iter, "iteration", mps_prefix);
        h5ppFile.writeAttribute(status.step, "step", mps_prefix);
        tools::common::profile::get_default_prof()["t_hdf"]->toc();
    }

    /*! Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation. */
    if(storage_level < StorageLevel::FULL) {
        save_log[mps_prefix] = save_point;
        return;
    }

    if(save_log[mps_prefix] != save_point) {
        tools::log->trace("Storing [{: ^6}]: mps tensors", enum2str(storage_level));
        tools::common::profile::get_default_prof()["t_hdf"]->tic();
        for(size_t pos = 0; pos < state.get_length(); pos++) {
            dsetName = fmt::format("{}/M_{}", mps_prefix, pos);
            if(save_log[dsetName] == save_point) continue;
            h5ppFile.writeDataset(state.get_mps_site(pos).get_M_bare(), dsetName, layout); // Important to write bare matrices!!
            h5ppFile.writeAttribute(pos, "position", dsetName);
            h5ppFile.writeAttribute(state.get_mps_site(pos).get_M_bare().dimensions(), "dimensions", dsetName);
            h5ppFile.writeAttribute(state.get_mps_site(pos).get_label(), "label", dsetName);
            save_log[dsetName] = save_point;
        }
        tools::common::profile::get_default_prof()["t_hdf"]->toc();
        save_log[mps_prefix] = save_point;
    }
}

/*! Write all the MPO's with site info in attributes */
void tools::finite::io::h5dset::save_model(h5pp::File &h5ppFile, const std::string &model_prefix, const StorageLevel &storage_level,
                                           const class_model_finite &model) {
    if(storage_level < StorageLevel::FULL) return;
    // We do not expect the MPO's to change. Therefore if they exist, there is nothing else to do here
    if(h5ppFile.linkExists(model_prefix)) return tools::log->trace("The MPO's have already been written to [{}]", model_prefix);
    tools::log->trace("Storing [{: ^6}]: mpo tensors", enum2str(storage_level));
    tools::common::profile::get_default_prof()["t_hdf"]->tic();
    for(size_t pos = 0; pos < model.get_length(); pos++) { model.get_mpo(pos).save_mpo(h5ppFile, model_prefix); }
    h5ppFile.writeAttribute(settings::model::model_size, "model_size", model_prefix);
    h5ppFile.writeAttribute(enum2str(settings::model::model_type), "model_type", model_prefix);
    tools::common::profile::get_default_prof()["t_hdf"]->toc();
}

/*! Write down measurements that can't fit in a table */
void tools::finite::io::h5dset::save_entgm(h5pp::File &h5ppFile, const std::string &state_prefix, const StorageLevel &storage_level,
                                           const class_state_finite &state, const class_algorithm_status &status) {
    if(storage_level < StorageLevel::NORMAL) return;
    state.do_all_measurements();
    tools::common::profile::get_default_prof()["t_hdf"]->tic();

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

    tools::log->trace("Storing [{: ^6}]: truncation errors", enum2str(storage_level));
    h5ppFile.writeDataset(state.get_truncation_errors(), state_prefix + "/truncation_errors");
    tools::common::profile::get_default_prof()["t_hdf"]->toc();
}
