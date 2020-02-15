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

void tools::finite::io::h5dset::write_all_state(const class_state_finite &state, h5pp::File &h5ppFile, const std::string &prefix_path) {
    tools::log->trace("Writing state to datasets in path: {}...", prefix_path + "/state");
    tools::common::profile::t_hdf->tic();

    if(settings::output::storage_level >= StorageLevel::LIGHT) {
        write_bond_matrix(state, h5ppFile, prefix_path);
        h5ppFile.writeDataset(state.get_length(), prefix_path + "/state/sites");
        h5ppFile.writeDataset(state.get_position(), prefix_path + "/state/position");
    }
    if(settings::output::storage_level >= StorageLevel::NORMAL) {
        write_bond_matrices(state, h5ppFile, prefix_path);
    }
    if(settings::output::storage_level >= StorageLevel::FULL) {
        write_full_mps(state, h5ppFile, prefix_path);
        write_full_mpo(state, h5ppFile, prefix_path);
    }

    tools::common::profile::t_hdf->toc();
    tools::log->trace("Writing state to datasets in path: {}... OK", prefix_path + "/state");
}

void tools::finite::io::h5dset::write_bond_matrices(const class_state_finite &state, h5pp::File &h5ppFile, const std::string &prefix_path)
/*! Writes down all the "Lambda" bond matrices (singular value matrices), so we can obtain the entanglement spectrum easily.
 */
{
    H5D_layout_t layout = decide_layout(prefix_path);
    auto         middle = (size_t)(state.get_length() / 2);
    for(size_t i = 0; i < state.get_length(); i++) {
        h5ppFile.writeDataset(state.get_MPS(i).get_L(), prefix_path + "/state/mps/L_" + std::to_string(i), layout);
        if(i == middle) {
            h5ppFile.writeDataset(state.midchain_bond(), prefix_path + "/state/mps/L_C", layout);
        }
    }
    h5ppFile.writeDataset(state.get_truncation_errors(), prefix_path + "/state/truncation_errors", layout);
}

void tools::finite::io::h5dset::write_bond_matrix(const class_state_finite &state, h5pp::File &h5ppFile, const std::string &prefix_path)
/*! Writes down the center "Lambda" bond matrix (singular values), so we can obtain the entanglement spectrum easily.
 */
{
    H5D_layout_t layout = decide_layout(prefix_path);
    auto         middle = (size_t)(state.get_length() / 2);
    h5ppFile.writeDataset(state.midchain_bond(), prefix_path + "/state/mps/L_C", layout);
    h5ppFile.writeDataset(state.get_truncation_error(middle), prefix_path + "/state/truncation_error", layout);
}

void tools::finite::io::h5dset::write_full_mps(const class_state_finite &state, h5pp::File &h5ppFile, const std::string &prefix_path)
/*!
 * Writes down the full MPS in "L-G-L-G- LC -G-L-G-L" notation.
 *
 */
{
    H5D_layout_t layout = decide_layout(prefix_path);
    write_bond_matrices(state, h5ppFile, prefix_path);
    for(size_t i = 0; i < state.get_length(); i++) {
        h5ppFile.writeDataset(state.get_MPS(i).get_M(), prefix_path + "/state/mps/M_" + std::to_string(i), layout);
    }
}

void tools::finite::io::h5dset::write_full_mpo(const class_state_finite &state, h5pp::File &h5ppFile, const std::string &prefix_path) {
    // Write all the MPO's with site info in attributes
    if(settings::output::storage_level == StorageLevel::NONE) return;
    H5D_layout_t layout = decide_layout(prefix_path);
    for(auto site = 0ul; site < state.get_length(); site++) {
        std::string dataset_name = prefix_path + "/state/mpo/H_" + std::to_string(site);
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


}

void tools::finite::io::h5dset::write_model(const class_state_finite &state, h5pp::File &h5ppFile, const std::string &prefix_path) {
    // Write down the Hamiltonian model type and site info as attributes
    if(settings::output::storage_level == StorageLevel::NONE) return;
    h5ppFile.writeDataset(settings::model::model_type, prefix_path + "/model/model_type");
    for(auto site = 0ul; site < state.get_length(); site++) {
        // Write MPO properties as attributes
        std::string dataset_name = prefix_path + "/model/hamiltonian_" + std::to_string(site);
        h5ppFile.writeDataset(site, dataset_name);
        for(auto &params : state.get_MPO(site).get_parameters()) {
            if(params.second.type() == typeid(double)) h5ppFile.writeAttribute(std::any_cast<double>(params.second), params.first, dataset_name);
            if(params.second.type() == typeid(size_t)) h5ppFile.writeAttribute(std::any_cast<size_t>(params.second), params.first, dataset_name);
            if(params.second.type() == typeid(int)) h5ppFile.writeAttribute(std::any_cast<int>(params.second), params.first, dataset_name);
            if(params.second.type() == typeid(bool)) h5ppFile.writeAttribute(std::any_cast<bool>(params.second), params.first, dataset_name);
            if(params.second.type() == typeid(std::string)) h5ppFile.writeAttribute(std::any_cast<std::string>(params.second), params.first, dataset_name);
        }
    }
}

void tools::finite::io::h5dset::write_array_measurements(const class_state_finite &state, h5pp::File &h5ppFile, const std::string &prefix_path) {
    state.do_all_measurements();
    tools::log->trace("Writing all measurements...");
    tools::common::profile::t_hdf->tic();
    h5ppFile.writeDataset(state.measurements.bond_dimensions.value(), prefix_path + "/bond_dimensions");
    h5ppFile.writeDataset(state.measurements.entanglement_entropies.value(), prefix_path + "/entanglement_entropies");
    h5ppFile.writeDataset(state.measurements.spin_components.value(), prefix_path + "/spin_components");
    tools::common::profile::t_hdf->toc();
    tools::log->trace("Writing all measurements... OK");
}
