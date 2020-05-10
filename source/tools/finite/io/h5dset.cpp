//
// Created by david on 2019-03-09.
//

#include <h5pp/h5pp.h>
#include <model/class_model_factory.h>
#include <regex>
#include <simulation/class_simulation_status.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/io.h>
#include <tools/finite/measure.h>
using Scalar = std::complex<double>;

H5D_layout_t decide_layout(std::string_view prefix_path) {
    std::string str(prefix_path);
    std::regex  rx(R"(journal/iter_[0-9])"); // Declare the regex with a raw string literal
    std::smatch m;
    if(regex_search(str, m, rx))
        return H5D_CONTIGUOUS;
    else
        return H5D_CHUNKED;
}

void tools::finite::io::h5dset::write_mps(h5pp::File &h5ppFile, const std::string &state_prefix, const StorageLevel &storage_level, const class_state_finite &state) {
    if(storage_level == StorageLevel::NONE) return;

    /*! Writes down midchain and/or all the "Lambda" bond matrices (singular value matrices), so we can obtain the entanglement spectrum easily. */
    tools::log->trace("Storing [{: ^6}]: mid bond matrix", enum2str(storage_level));
    tools::common::profile::t_hdf->tic();
    H5D_layout_t layout = decide_layout(state_prefix);
    /*! Writes down the center "Lambda" bond matrix (singular values). */
    std::string dsetName = state_prefix + "/schmidt_midchain";
    h5ppFile.writeDataset(state.midchain_bond(), dsetName, layout);
    h5ppFile.writeAttribute(state.get_truncation_error_midchain(), "truncation_error", dsetName);
    h5ppFile.writeAttribute((state.get_length() - 1) / 2, "position", dsetName);
    h5ppFile.writeAttribute(state.get_iteration(), "iteration", dsetName);
    h5ppFile.writeAttribute(state.get_step(), "step", dsetName);
    h5ppFile.writeAttribute(state.get_chi_lim(), "chi_lim", dsetName);
    h5ppFile.writeAttribute(state.get_chi_max(), "chi_max", dsetName);
    tools::common::profile::t_hdf->toc();

    if(storage_level < StorageLevel::NORMAL) return;
    tools::log->trace("Storing [{: ^6}]: all bond matrices", enum2str(storage_level));
    tools::common::profile::t_hdf->tic();
    // There should be one more sites+1 number of L's, because there is also a center bond
    // However L_i always belongs M_i. Stick to this rule!
    // This means that some M_i has two bonds, one L_i to the left, and one L_C to the right.
    std::string mps_path = state_prefix + "/mps";
    for(size_t pos = 0; pos < state.get_length(); pos++) {
        dsetName = mps_path + "/L_" + std::to_string(pos);
        h5ppFile.writeDataset(state.get_MPS(pos).get_L(), dsetName, layout);
        h5ppFile.writeAttribute(pos, "position", dsetName);
        h5ppFile.writeAttribute(state.get_MPS(pos).get_L().dimensions(), "dimensions", dsetName);
        if(state.get_MPS(pos).isCenter()) {
            dsetName = mps_path + "/L_C";
            h5ppFile.writeDataset(state.get_MPS(pos).get_LC(), dsetName, layout);
            h5ppFile.writeAttribute(pos, "position", dsetName);
            h5ppFile.writeAttribute(state.get_MPS(pos).get_LC().dimensions(), "dimensions", dsetName);
        }
    }
    h5ppFile.writeAttribute(state.get_length(), "model_size", mps_path);
    h5ppFile.writeAttribute(state.get_position(), "position", mps_path);
    h5ppFile.writeAttribute(state.get_iteration(), "iteration", mps_path);
    h5ppFile.writeAttribute(state.get_step(), "step", mps_path);
    h5ppFile.writeAttribute(state.get_chi_lim(), "chi_lim", mps_path);
    h5ppFile.writeAttribute(state.get_chi_max(), "chi_max", mps_path);
    h5ppFile.writeAttribute(state.get_truncation_errors(), "truncation_errors", mps_path);
    tools::common::profile::t_hdf->toc();

    /*! Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation. */
    if(storage_level < StorageLevel::FULL) return;
    tools::log->trace("Storing [{: ^6}]: mps tensors", enum2str(storage_level));
    tools::common::profile::t_hdf->tic();
    for(size_t pos = 0; pos < state.get_length(); pos++) {
        dsetName = mps_path + "/M_" + std::to_string(pos);
        h5ppFile.writeDataset(state.get_MPS(pos).get_M_bare(), dsetName, layout); // Important to write bare matrices!!
        h5ppFile.writeAttribute(pos, "position", dsetName);
        h5ppFile.writeAttribute(state.get_MPS(pos).get_M_bare().dimensions(), "dimensions", dsetName);
    }
    tools::common::profile::t_hdf->toc();
}

/*! Write all the MPO's with site info in attributes */
void tools::finite::io::h5dset::write_mpo(h5pp::File &h5ppFile, const std::string &model_prefix, const StorageLevel &storage_level, const class_state_finite &state) {
    if(storage_level < StorageLevel::FULL) return;
    // We do not expect the MPO's to change. Therefore if they exist, there is nothing else to do here
    if(h5ppFile.linkExists(model_prefix + "/mpo")) return tools::log->trace("The MPO's have already been written to [{}]", model_prefix + "/mpo");


    tools::log->trace("Storing [{: ^6}]: mpo tensors", enum2str(storage_level));
    tools::common::profile::t_hdf->tic();
    for(size_t pos = 0; pos < state.get_length(); pos++) {
        state.get_MPO(pos).write_mpo(h5ppFile, model_prefix);
    }
    h5ppFile.writeAttribute(settings::model::model_size, "model_size", model_prefix + "/mpo");
    h5ppFile.writeAttribute(settings::model::model_type, "model_type", model_prefix + "/mpo");
    tools::common::profile::t_hdf->toc();
}

/*! Write down measurements that can't fit in a table */
void tools::finite::io::h5dset::write_measurements(h5pp::File &h5ppFile, const std::string &prefix, const StorageLevel &storage_level,
                                                   const class_state_finite &state) {
    if(storage_level < StorageLevel::NORMAL) return;
    state.do_all_measurements();
    tools::common::profile::t_hdf->tic();
    tools::log->trace("Storing [{: ^6}]: bond dims", enum2str(storage_level));
    h5ppFile.writeDataset(tools::finite::measure::bond_dimensions(state), prefix + "/bond_dimensions");
    tools::log->trace("Storing [{: ^6}]: entanglement entropies", enum2str(storage_level));
    h5ppFile.writeDataset(tools::finite::measure::entanglement_entropies(state), prefix + "/entanglement_entropies");
    tools::log->trace("Storing [{: ^6}]: truncation errors", enum2str(storage_level));
    h5ppFile.writeDataset(state.get_truncation_errors(), prefix + "/truncation_errors");
    tools::common::profile::t_hdf->toc();
}
