//
// Created by david on 2019-03-09.
//

#include <h5pp/h5pp.h>
#include <regex>
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/prof.h>
#include <tools/infinite/io.h>

int tools::infinite::io::h5dset::decide_layout(std::string_view prefix_path) {
    std::string str(prefix_path);
    std::regex  rx(R"(checkpoint/iter_[0-9])"); // Declare the regex with a raw string literal
    std::smatch m;
    if(regex_search(str, m, rx))
        return H5D_CONTIGUOUS;
    else
        return H5D_CHUNKED;
}


void tools::infinite::io::h5dset::write_state(h5pp::File &h5ppFile, const std::string &state_prefix,const StorageLevel & storage_level, const class_state_infinite &state){
    if(storage_level == StorageLevel::NONE) return;
    tools::log->trace("Storing [{: ^6}]: mid bond matrix", enum2str(storage_level));
    tools::common::profile::t_hdf->tic();
    auto layout = static_cast<H5D_layout_t>(decide_layout(state_prefix));
    std::string dsetName = state_prefix + "/schmidt_midchain";
    h5ppFile.writeDataset(state.LC(), dsetName, layout);
    h5ppFile.writeAttribute(state.get_truncation_error(), "truncation_error", dsetName);
    h5ppFile.writeAttribute(state.get_chi_lim(), "chi_lim", dsetName);
    h5ppFile.writeAttribute(state.get_chi_max(), "cfg_chi_lim_max", dsetName);
    tools::common::profile::t_hdf->toc();

    if(storage_level < StorageLevel::NORMAL) return;

    tools::log->trace("Storing [{: ^6}]: bond matrices", enum2str(storage_level));
    tools::common::profile::t_hdf->tic();
    h5ppFile.writeDataset(state.LC() , state_prefix + "/mps/L_C");
    h5ppFile.writeAttribute(state.get_truncation_error(), "truncation_error", state_prefix + "/mps/L_C");
    h5ppFile.writeAttribute(state.LC().dimensions(), "dimensions", state_prefix + "/mps/L_C");
    h5ppFile.writeDataset(state.LA() , state_prefix + "/mps/L_A");
    h5ppFile.writeDataset(state.LB() , state_prefix + "/mps/L_B");
    h5ppFile.writeAttribute(state.get_chi_lim(), "chi_lim", state_prefix + "/mps/L_C");
    h5ppFile.writeAttribute(state.get_chi_max(), "cfg_chi_lim_max", state_prefix + "/mps/L_C");
    tools::common::profile::t_hdf->toc();

    if(storage_level < StorageLevel::FULL) return;

    tools::common::profile::t_hdf->tic();
    h5ppFile.writeDataset(state.A_bare(), state_prefix + "/mps/M_A");
    h5ppFile.writeAttribute(state.A_bare().dimensions(), "dimensions", state_prefix + "/mps/M_A");
    h5ppFile.writeDataset(state.B() , state_prefix + "/mps/M_B");
    h5ppFile.writeAttribute(state.B().dimensions(), "dimensions", state_prefix + "/mps/M_B");
    tools::common::profile::t_hdf->toc();
}

void tools::infinite::io::h5dset::write_model(h5pp::File &h5ppFile, const std::string & model_prefix, const StorageLevel & storage_level, const class_model_infinite &model){
    if(storage_level < StorageLevel::FULL) return;
    if(h5ppFile.linkExists(model_prefix + "/mpo")) return tools::log->trace("The model has already been written to [{}]", model_prefix + "/mpo");

    tools::common::profile::t_hdf->tic();
    model.get_mpo_siteA().write_mpo(h5ppFile,model_prefix);
    model.get_mpo_siteB().write_mpo(h5ppFile,model_prefix);
    h5ppFile.writeAttribute(2, "model_size", model_prefix + "/mpo");
    h5ppFile.writeAttribute(enum2str(settings::model::model_type), "model_type", model_prefix + "/mpo");
    tools::common::profile::t_hdf->toc();
}

void tools::infinite::io::h5dset::write_edges(h5pp::File &h5ppFile, const std::string &edges_prefix, const StorageLevel & storage_level, const class_edges_infinite &edges){
    if(storage_level < StorageLevel::NORMAL) return;
    tools::common::profile::t_hdf->tic();
    const auto & ene = edges.get_ene_blk();
    const auto & var = edges.get_var_blk();
    h5ppFile.writeDataset(ene.L, edges_prefix + "/eneL");
    h5ppFile.writeDataset(ene.R, edges_prefix + "/eneR");
    h5ppFile.writeDataset(var.L, edges_prefix + "/varL");
    h5ppFile.writeDataset(var.R, edges_prefix + "/varR");
    tools::common::profile::t_hdf->toc();
}


//
//void tools::infinite::io::h5dset::write_hamiltonian_params(h5pp::File &h5pp_file, const std::string & sim_name, const StorageLevel & storage_level, const class_state_infinite &state){
//    tools::common::profile::t_hdf->tic();
//    h5pp_file.writeDataset(enum2str(settings::model::model_type), sim_name + "/model/model_type");
//    auto paramsA = state.HA->get_parameters();
//    auto paramsB = state.HB->get_parameters();
//        // Write MPO properties as attributes
//    std::string dataset_name = sim_name + "/model/hamiltonian_A";
//    h5pp_file.writeDataset("A", dataset_name);
//    for(auto &params : state.HA->get_parameters()) {
//        if(params.second.type() == typeid(double)) h5pp_file.writeAttribute(std::any_cast<double>(params.second), params.first, dataset_name);
//        if(params.second.type() == typeid(size_t)) h5pp_file.writeAttribute(std::any_cast<size_t>(params.second), params.first, dataset_name);
//        if(params.second.type() == typeid(int)) h5pp_file.writeAttribute(std::any_cast<int>(params.second), params.first, dataset_name);
//        if(params.second.type() == typeid(bool)) h5pp_file.writeAttribute(std::any_cast<bool>(params.second), params.first, dataset_name);
//        if(params.second.type() == typeid(std::string)) h5pp_file.writeAttribute(std::any_cast<std::string>(params.second), params.first, dataset_name);
//    }
//    dataset_name = sim_name + "/model/hamiltonian_B";
//    h5pp_file.writeDataset("B", dataset_name);
//    for(auto &params : state.HB->get_parameters()) {
//        if(params.second.type() == typeid(double)) h5pp_file.writeAttribute(std::any_cast<double>(params.second), params.first, dataset_name);
//        if(params.second.type() == typeid(size_t)) h5pp_file.writeAttribute(std::any_cast<size_t>(params.second), params.first, dataset_name);
//        if(params.second.type() == typeid(int)) h5pp_file.writeAttribute(std::any_cast<int>(params.second), params.first, dataset_name);
//        if(params.second.type() == typeid(bool)) h5pp_file.writeAttribute(std::any_cast<bool>(params.second), params.first, dataset_name);
//        if(params.second.type() == typeid(std::string)) h5pp_file.writeAttribute(std::any_cast<std::string>(params.second), params.first, dataset_name);
//    }
//    tools::common::profile::t_hdf->toc();
//}
//
//void tools::infinite::io::h5dset::write_all_measurements  (h5pp::File &h5pp_file, const std::string & sim_name, const StorageLevel & storage_level, const class_simulation_status & status, const class_state_infinite &state){
//    state.do_all_measurements();
//    tools::common::profile::t_hdf->tic();
//    h5pp_file.writeDataset(state.measurements.length.value()                      , sim_name + "/measurements/2site/length");
//    h5pp_file.writeDataset(state.measurements.bond_dimension.value()              , sim_name + "/measurements/2site/bond_dimension");
//    h5pp_file.writeDataset(state.measurements.norm.value()                        , sim_name + "/measurements/2site/norm");
//    h5pp_file.writeDataset(state.measurements.truncation_error.value()            , sim_name + "/measurements/2site/truncation_error");
//    h5pp_file.writeDataset(state.measurements.energy_mpo.value()                  , sim_name + "/measurements/2site/energy");
//    h5pp_file.writeDataset(state.measurements.energy_per_site_mpo.value()         , sim_name + "/measurements/2site/energy_per_site");
//    h5pp_file.writeDataset(state.measurements.energy_per_site_ham.value()         , sim_name + "/measurements/2site/energy_per_site_mom");
//    h5pp_file.writeDataset(state.measurements.energy_per_site_mom.value()         , sim_name + "/measurements/2site/energy_per_site_mom");
//    h5pp_file.writeDataset(state.measurements.energy_variance_per_site_mpo.value(), sim_name + "/measurements/2site/energy_variance_per_site");
//    h5pp_file.writeDataset(state.measurements.energy_variance_per_site_ham.value(), sim_name + "/measurements/2site/energy_variance_per_site_ham");
//    h5pp_file.writeDataset(state.measurements.energy_variance_per_site_mom.value(), sim_name + "/measurements/2site/energy_variance_per_site_mom");
//    h5pp_file.writeDataset(state.measurements.entanglement_entropy.value(), sim_name + "/measurements/2site/entanglement_entropy_midchain");
//    tools::common::profile::t_hdf->toc();
//}

