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

int tools::finite::io::h5dset::decide_layout(std::string_view prefix_path) {
    std::string str(prefix_path);
    std::regex  rx(R"(checkpoint/iter_[0-9])"); // Declare the regex with a raw string literal
    std::smatch m;
    if(regex_search(str, m, rx)) return H5D_CONTIGUOUS;
    else
        return H5D_CHUNKED;
}

void tools::finite::io::h5dset::save_state(h5pp::File &h5ppFile, const std::string &state_prefix, const StorageLevel &storage_level,
                                           const class_state_finite &state) {
    if(storage_level == StorageLevel::NONE) return;
    // Checks if the current entry has already been saved
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    auto                                                                  save_point = std::make_pair(state.get_iteration(), state.get_step());

    auto layout = static_cast<H5D_layout_t>(decide_layout(state_prefix));

    /*! Writes down the center "Lambda" bond matrix (singular values). */
    std::string dsetName = state_prefix + "/schmidt_midchain";
    if(save_log[dsetName] != save_point) {
        tools::log->trace("Storing [{: ^6}]: mid bond matrix", enum2str(storage_level));
        tools::common::profile::get_default_prof()["t_hdf"]->tic();
        h5ppFile.writeDataset(state.midchain_bond(), dsetName, layout);
        h5ppFile.writeAttribute(state.get_truncation_error_midchain(), "truncation_error", dsetName);
        h5ppFile.writeAttribute((state.get_length() - 1) / 2, "position", dsetName);
        h5ppFile.writeAttribute(state.get_iteration(), "iteration", dsetName);
        h5ppFile.writeAttribute(state.get_step(), "step", dsetName);
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
            if(state.get_mps_site(pos).isCenter()) {
                dsetName = mps_prefix + "/L_C";
                h5ppFile.writeDataset(state.get_mps_site(pos).get_LC(), dsetName, layout);
                h5ppFile.writeAttribute(pos, "position", dsetName);
                h5ppFile.writeAttribute(state.get_mps_site(pos).get_LC().dimensions(), "dimensions", dsetName);
            }
            save_log[dsetName] = save_point;
        }
        h5ppFile.writeAttribute(state.get_length(), "model_size", mps_prefix);
        h5ppFile.writeAttribute(state.get_position(), "position", mps_prefix);
        h5ppFile.writeAttribute(state.get_iteration(), "iteration", mps_prefix);
        h5ppFile.writeAttribute(state.get_step(), "step", mps_prefix);
        h5ppFile.writeAttribute(state.get_truncation_errors(), "truncation_errors", mps_prefix);
        tools::common::profile::get_default_prof()["t_hdf"]->toc();
    }

    /*! Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation. */
    if(storage_level < StorageLevel::FULL) {
        save_log[mps_prefix] = save_point;
        return;
    }

    if(save_log[mps_prefix] != save_point){
        tools::log->trace("Storing [{: ^6}]: mps tensors", enum2str(storage_level));
        tools::common::profile::get_default_prof()["t_hdf"]->tic();
        for(size_t pos = 0; pos < state.get_length(); pos++) {
            dsetName = fmt::format("{}/M_{}", mps_prefix, pos);
            if(save_log[dsetName] == save_point) continue;
            h5ppFile.writeDataset(state.get_mps_site(pos).get_M_bare(), dsetName, layout); // Important to write bare matrices!!
            h5ppFile.writeAttribute(pos, "position", dsetName);
            h5ppFile.writeAttribute(state.get_mps_site(pos).get_M_bare().dimensions(), "dimensions", dsetName);
            save_log[dsetName] = save_point;
        }
        tools::common::profile::get_default_prof()["t_hdf"]->toc();
        save_log[mps_prefix] = save_point;
    }
}

/*! Write all the MPO's with site info in attributes */
void tools::finite::io::h5dset::save_model(h5pp::File &h5ppFile, const std::string &mpo_path, const StorageLevel &storage_level,
                                           const class_model_finite &model) {
    if(storage_level < StorageLevel::FULL) return;
    // We do not expect the MPO's to change. Therefore if they exist, there is nothing else to do here
    if(h5ppFile.linkExists(mpo_path)) return tools::log->trace("The MPO's have already been written to [{}]", mpo_path);
    tools::log->trace("Storing [{: ^6}]: mpo tensors", enum2str(storage_level));
    tools::common::profile::get_default_prof()["t_hdf"]->tic();
    for(size_t pos = 0; pos < model.get_length(); pos++) { model.get_mpo(pos).save_mpo(h5ppFile, mpo_path); }
    h5ppFile.writeAttribute(settings::model::model_size, "model_size", mpo_path);
    h5ppFile.writeAttribute(enum2str(settings::model::model_type), "model_type", mpo_path);
    tools::common::profile::get_default_prof()["t_hdf"]->toc();
}

/*! Write down measurements that can't fit in a table */
void tools::finite::io::h5dset::save_bonds(h5pp::File &h5ppFile, const std::string &state_prefix, const StorageLevel &storage_level,
                                           const class_state_finite &state) {
    if(storage_level < StorageLevel::NORMAL) return;
    state.do_all_measurements();
    tools::common::profile::get_default_prof()["t_hdf"]->tic();
    tools::log->trace("Storing [{: ^6}]: bond dims", enum2str(storage_level));
    h5ppFile.writeDataset(tools::finite::measure::bond_dimensions(state), state_prefix + "/bond_dimensions");
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
