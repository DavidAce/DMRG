//
// Created by david on 2019-03-09.
//

#include <tools/finite/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <h5pp/h5pp.h>
#include <general/nmspc_quantum_mechanics.h>
#include <model/class_model_factory.h>
#include <simulation/class_simulation_status.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_finite.h>
#include <regex>
using Scalar = std::complex<double>;

H5D_layout_t decide_layout(std::string_view prefix_path) {
    std::string str(prefix_path);
    std::regex rx(R"(journal/iter_[0-9])"); // Declare the regex with a raw string literal
    std::smatch m;
    if (regex_search(str, m, rx))
        return H5D_CONTIGUOUS;
    else
        return H5D_CHUNKED;
}

void tools::finite::io::h5dset::write_all(h5pp::File & h5ppFile, const std::string &prefix, const StorageLevel & storage_level, const class_state_finite & state) {
    if(storage_level == StorageLevel::NONE) return;

    tools::log->trace("Writing {} state to prefix: {}", enum2str(storage_level), prefix);
    write_mps(h5ppFile, prefix, storage_level, state);
    write_mpo(h5ppFile, prefix, storage_level, state);
}

void tools::finite::io::h5dset::write_mps(h5pp::File & h5ppFile, const std::string &prefix, const StorageLevel & storage_level, const class_state_finite & state){
    if(storage_level < StorageLevel::LIGHT) return;
    /*! Writes down midchain and/or all the "Lambda" bond matrices (singular value matrices), so we can obtain the entanglement spectrum easily. */
    tools::log->trace("Writing mid bond matrix");
    tools::common::profile::t_hdf->tic();
    H5D_layout_t layout = decide_layout(prefix);
    /*! Writes down the center "Lambda" bond matrix (singular values). */
    std::string dsetName = prefix + "/schmidt_midchain";
    h5ppFile.writeDataset(state.midchain_bond(), dsetName , layout);
    h5ppFile.writeAttribute(state.get_truncation_error_midchain(), "truncation_error" , dsetName);
    h5ppFile.writeAttribute((state.get_length()-1) / 2, "position" , dsetName);
    h5ppFile.writeAttribute(state.get_iteration(), "iteration", dsetName);
    h5ppFile.writeAttribute(state.get_step(), "step", dsetName);
    h5ppFile.writeAttribute(state.get_chi_lim(),"chi_lim",dsetName);
    h5ppFile.writeAttribute(state.get_chi_max(),"chi_max",dsetName);
    tools::common::profile::t_hdf->toc();


    if(storage_level < StorageLevel::NORMAL) return;
    tools::log->trace("Writing all bond matrices");
    tools::common::profile::t_hdf->tic();
    // There should be one more sites+1 number of L's, because there is also a center bond
    // However L_i always belongs M_i. Stick to this rule!
    // This means that some M_i has two bonds, one L_i to the left, and one L_C to the right.
    for(size_t pos = 0; pos < state.get_length(); pos++) {
        dsetName = prefix + "/mps/L_" + std::to_string(pos);
        h5ppFile.writeDataset(state.get_MPS(pos).get_L(),dsetName , layout);
        h5ppFile.writeAttribute(pos, "position" , dsetName);
        h5ppFile.writeAttribute( state.get_MPS(pos).get_L().dimensions(), "dimensions", dsetName);
        if(state.get_MPS(pos).isCenter()){
            dsetName = prefix + "/mps/L_C";
            h5ppFile.writeDataset(state.get_MPS(pos).get_LC(), dsetName, layout);
            h5ppFile.writeAttribute(pos, "position" , dsetName);
            h5ppFile.writeAttribute(state.get_MPS(pos).get_LC().dimensions(), "dimensions" , dsetName);
        }
    }
    h5ppFile.writeAttribute(state.get_length(), "sites" , prefix + "/mps");
    h5ppFile.writeAttribute(state.get_position(), "position", prefix + "/mps");
    h5ppFile.writeAttribute(state.get_iteration(), "iteration", prefix + "/mps");
    h5ppFile.writeAttribute(state.get_step(), "step", prefix + "/mps");
    h5ppFile.writeAttribute(state.get_chi_lim(),"chi_lim",prefix + "/mps");
    h5ppFile.writeAttribute(state.get_chi_max(),"chi_max",prefix + "/mps");
    h5ppFile.writeAttribute(state.get_truncation_errors(),"truncation_errors",prefix + "/mps");
    tools::common::profile::t_hdf->toc();


    /*! Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation. */
    if(storage_level < StorageLevel::FULL) return;
    tools::log->trace("Writing MPS tensors");
    tools::common::profile::t_hdf->tic();
    for(size_t pos = 0; pos < state.get_length(); pos++) {
        dsetName = prefix + "/mps/M_" + std::to_string(pos);
        h5ppFile.writeDataset(state.get_MPS(pos).get_M(),dsetName, layout);
        h5ppFile.writeAttribute(pos, "position" , dsetName);
        h5ppFile.writeAttribute(state.get_MPS(pos).get_M().dimensions(), "dimensions" , dsetName);
    }
    tools::common::profile::t_hdf->toc();
}

/*! Write all the MPO's with site info in attributes */
void tools::finite::io::h5dset::write_mpo(h5pp::File & h5ppFile, const std::string &prefix, const StorageLevel & storage_level, const class_state_finite & state) {
    if(storage_level < StorageLevel::FULL) return;
    tools::log->trace("Writing MPO tensors");
    tools::common::profile::t_hdf->tic();
    H5D_layout_t layout = decide_layout(prefix);
    for(auto site = 0ul; site < state.get_length(); site++) {
        std::string dataset_name = prefix + "/mpo/H_" + std::to_string(site);
        h5ppFile.writeDataset(state.get_MPO(site).MPO(), dataset_name, layout);
        // Write MPO properties as attributes
        for(auto &params : state.get_MPO(site).get_parameters()) {
            if(params.second.type() == typeid(double)) h5ppFile.writeAttribute(std::any_cast<double>(params.second), params.first, dataset_name);
            if(params.second.type() == typeid(size_t)) h5ppFile.writeAttribute(std::any_cast<size_t>(params.second), params.first, dataset_name);
            if(params.second.type() == typeid(int)) h5ppFile.writeAttribute(std::any_cast<int>(params.second), params.first, dataset_name);
            if(params.second.type() == typeid(bool)) h5ppFile.writeAttribute(std::any_cast<bool>(params.second), params.first, dataset_name);
            if(params.second.type() == typeid(std::string)) h5ppFile.writeAttribute(std::any_cast<std::string>(params.second), params.first, dataset_name);
        }
    }
    h5ppFile.writeAttribute(state.get_length(), "sites" , prefix + "/mpo");
    h5ppFile.writeAttribute(state.get_position(), "position", prefix + "/mpo");
    h5ppFile.writeAttribute(state.get_iteration(), "iteration", prefix + "/mpo");
    h5ppFile.writeAttribute(state.get_step(), "step", prefix + "/mpo");
    tools::common::profile::t_hdf->toc();
}

/*! Write down measurements that can't fit in a table */
void tools::finite::io::h5dset::write_array_measurements(h5pp::File & h5ppFile, const std::string &prefix, const StorageLevel & storage_level, const class_state_finite & state) {
    if(storage_level < StorageLevel::NORMAL) return;
    state.do_all_measurements();
    tools::log->trace("Writing array measurements");
    tools::common::profile::t_hdf->tic();
    h5ppFile.writeDataset(tools::finite::measure::bond_dimensions(state), prefix + "/bond_dimensions");
    h5ppFile.writeDataset(tools::finite::measure::entanglement_entropies(state), prefix + "/entanglement_entropies");
    h5ppFile.writeDataset(state.get_truncation_errors(), prefix + "/truncation_errors");
    tools::common::profile::t_hdf->toc();
}
