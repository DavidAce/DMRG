//
// Created by david on 2019-03-09.
//

#include <tools/finite/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <h5pp/h5pp.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <model/class_model_factory.h>
#include <simulation/class_simulation_status.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_finite.h>
#include <state/class_state_infinite.h>
#include <stdexcept>

using Scalar = std::complex<double>;

H5D_layout_t decide_layout(std::string_view prefix_path) {
    if(prefix_path.find("/log") != std::string::npos) // A logged dataset is not supposed to be chunked, just written once.
        return H5D_CONTIGUOUS;
    else
        return H5D_CHUNKED;
}

void tools::finite::io::h5dset::write_all_state(h5pp::File & h5ppFile, const std::string &path, const StorageLevel & storage_level, const class_state_finite & state) {
    tools::log->trace("Writing state to path: {}...", path + "/state");
    write_bond_matrix(h5ppFile,path, storage_level,state);
    write_bond_matrices(h5ppFile,path, storage_level,state);
    h5ppFile.writeDataset(state.get_length(), path + "/state/sites");
    h5ppFile.writeDataset(state.get_position(), path + "/state/position");
    h5ppFile.writeDataset(state.get_truncation_error_midchain(), path + "/state/truncation_error_midchain");
    h5ppFile.writeDataset(state.get_truncation_error(), path + "/state/truncation_error_current");
    write_full_mps(h5ppFile,path, storage_level,state);
    write_full_mpo(h5ppFile,path, storage_level,state);
    tools::log->trace("Writing state to datasets in path: {}... OK", path + "/state");
}

/*! Writes down all the "Lambda" bond matrices (singular value matrices), so we can obtain the entanglement spectrum easily. */
void tools::finite::io::h5dset::write_bond_matrices(h5pp::File & h5ppFile, const std::string &path, const StorageLevel & storage_level, const class_state_finite & state){
    tools::common::profile::t_hdf->tic();
    H5D_layout_t layout = decide_layout(path);
    size_t count = 0; // There should be one more sites+1 number of L's, because there is also a center bond
    for(size_t i = 0; i < state.get_length(); i++) {
        h5ppFile.writeDataset(state.get_MPS(i).get_L(), path + "/state/mps/L_" + std::to_string(count++), layout);
        if(state.get_MPS(i).isCenter())
            h5ppFile.writeDataset(state.get_MPS(i).get_LC(), path + "/state/mps/L_" + std::to_string(count++), layout);
    }
    h5ppFile.writeDataset(state.get_truncation_errors(), path + "/state/truncation_errors", layout);
    tools::common::profile::t_hdf->toc();
}

/*! Writes down the center "Lambda" bond matrix (singular values), so we can obtain the entanglement spectrum easily. */
void tools::finite::io::h5dset::write_bond_matrix(h5pp::File & h5ppFile, const std::string &path, const StorageLevel & storage_level, const class_state_finite & state){
    tools::common::profile::t_hdf->tic();
    H5D_layout_t layout = decide_layout(path);
    h5ppFile.writeDataset(state.midchain_bond(), path + "/state/mps/L_C", layout);
    tools::common::profile::t_hdf->toc();
}


/*! Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation. */
void tools::finite::io::h5dset::write_full_mps(h5pp::File & h5ppFile, const std::string &path, const StorageLevel & storage_level, const class_state_finite & state){
    tools::common::profile::t_hdf->tic();
    H5D_layout_t layout = decide_layout(path);
    std::string mpstype = "A";
    for(size_t i = 0; i < state.get_length(); i++) {
        h5ppFile.writeDataset(state.get_MPS(i).get_M(), path + "/state/mps/" + mpstype + "_" + std::to_string(i), layout);
        if(state.get_MPS(i).isCenter())
            mpstype = "B";
    }
    tools::common::profile::t_hdf->toc();
}

/*! Write all the MPO's with site info in attributes */
void tools::finite::io::h5dset::write_full_mpo(h5pp::File & h5ppFile, const std::string &path, const StorageLevel & storage_level, const class_state_finite & state) {
    if(settings::output::storage_level_results == StorageLevel::NONE) return;
    tools::common::profile::t_hdf->tic();
    H5D_layout_t layout = decide_layout(path);
    for(auto site = 0ul; site < state.get_length(); site++) {
        std::string dataset_name = path + "/state/mpo/H_" + std::to_string(site);
        h5ppFile.writeDataset(state.get_MPO(site).MPO(), dataset_name, layout);
        // Write MPO properties as attributes
        h5ppFile.writeDataset(site, dataset_name);
        for(auto &params : state.get_MPO(site).get_parameters()) {
            if(params.second.type() == typeid(double)) h5ppFile.writeAttribute(std::any_cast<double>(params.second), params.first, dataset_name);
            if(params.second.type() == typeid(size_t)) h5ppFile.writeAttribute(std::any_cast<size_t>(params.second), params.first, dataset_name);
            if(params.second.type() == typeid(int)) h5ppFile.writeAttribute(std::any_cast<int>(params.second), params.first, dataset_name);
            if(params.second.type() == typeid(bool)) h5ppFile.writeAttribute(std::any_cast<bool>(params.second), params.first, dataset_name);
            if(params.second.type() == typeid(std::string)) h5ppFile.writeAttribute(std::any_cast<std::string>(params.second), params.first, dataset_name);
        }
    }
    tools::common::profile::t_hdf->toc();
}

/*! Write down the Hamiltonian model type and site info as attributes */
void tools::finite::io::h5dset::write_model(h5pp::File & h5ppFile, const std::string &path, const StorageLevel & storage_level, const class_state_finite & state) {
    if(settings::output::storage_level_results == StorageLevel::NONE) return;
    tools::common::profile::t_hdf->tic();
    h5ppFile.writeDataset(settings::model::model_type, path + "/model/model_type");
    std::string table_name = path + "/model/Hamiltonian";
    for(auto site = 0ul; site < state.get_length(); site++)
        state.get_MPO(site).write_parameters(h5ppFile,table_name);
    tools::common::profile::t_hdf->toc();
}

/*! Write down measurements that can't fit in a table */
void tools::finite::io::h5dset::write_array_measurements(h5pp::File & h5ppFile, const std::string &path, const StorageLevel & storage_level, const class_state_finite & state) {
    if(settings::output::storage_level_results == StorageLevel::NONE) return;
    state.do_all_measurements();
    tools::log->trace("Writing all measurements...");
    tools::common::profile::t_hdf->tic();
    h5ppFile.writeDataset(state.measurements.bond_dimensions.value(), path + "/bond_dimensions");
    h5ppFile.writeDataset(state.measurements.entanglement_entropies.value(), path + "/entanglement_entropies");
    h5ppFile.writeDataset(state.measurements.spin_components.value(), path + "/spin_components");
    tools::common::profile::t_hdf->toc();
    tools::log->trace("Writing all measurements... OK");
}
