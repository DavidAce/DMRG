//
// Created by david on 2019-03-09.
//

#include <algorithms/class_algorithm_status.h>
#include <h5pp/h5pp.h>
#include <regex>
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/prof.h>
#include <tools/infinite/io.h>

namespace tools::infinite::io::h5dset {
    void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile,
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
}

int tools::infinite::io::h5dset::decide_layout(std::string_view prefix_path) {
    std::string str(prefix_path);
    std::regex  rx(R"(checkpoint/iter_[0-9])"); // Declare the regex with a raw string literal
    std::smatch m;
    if(regex_search(str, m, rx)) return H5D_CONTIGUOUS;
    else
        return H5D_CHUNKED;
}

void tools::infinite::io::h5dset::save_state(h5pp::File &h5ppFile, const std::string &state_prefix, const StorageLevel &storage_level,
                                             const class_state_infinite &state, const class_algorithm_status &status) {
    if(storage_level == StorageLevel::NONE) return;

    // Checks if the current entry has already been saved
    // If it is empty because we are resuming, check if there is a log entry on file already
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    bootstrap_save_log(save_log, h5ppFile, {state_prefix + "/schmidt_midchain", state_prefix + "/mps"});
    auto save_point = std::make_pair(status.iter, status.step);
    auto        layout   = static_cast<H5D_layout_t>(decide_layout(state_prefix));

    std::string dsetName = state_prefix + "/schmidt_midchain";
    if(save_log[dsetName] != save_point) {
        tools::log->trace("Storing [{: ^6}]: mid bond matrix", enum2str(storage_level));
        tools::common::profile::get_default_prof()["t_hdf"]->tic();
        h5ppFile.writeDataset(state.LC(), dsetName, layout);
        h5ppFile.writeAttribute(state.get_truncation_error(), "truncation_error", dsetName);
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
        dsetName = mps_prefix + "/L_A";
        if(save_log[dsetName] != save_point){
            h5ppFile.writeDataset(state.LA(), dsetName);
            h5ppFile.writeAttribute(state.LA().dimensions(), "dimensions", dsetName);
            h5ppFile.writeAttribute(state.get_truncation_error(), "truncation_error", dsetName);
            h5ppFile.writeAttribute(status.chi_lim, "chi_lim", dsetName);
            h5ppFile.writeAttribute(status.chi_lim_max, "chi_lim_max", dsetName);
            save_log[dsetName] = save_point;
        }
        dsetName = mps_prefix + "/L_B";
        if(save_log[dsetName] != save_point){
            h5ppFile.writeDataset(state.LB(), dsetName);
            h5ppFile.writeAttribute(state.LB().dimensions(), "dimensions", dsetName);
            h5ppFile.writeAttribute(state.get_truncation_error(), "truncation_error", dsetName);
            h5ppFile.writeAttribute(status.chi_lim, "chi_lim", dsetName);
            h5ppFile.writeAttribute(status.chi_lim_max, "chi_lim_max", dsetName);
            save_log[dsetName] = save_point;
        }
        dsetName = mps_prefix + "/L_C";
        if(save_log[dsetName] != save_point){
            h5ppFile.writeDataset(state.LC(), dsetName);
            h5ppFile.writeAttribute(state.LC().dimensions(), "dimensions", dsetName);
            h5ppFile.writeAttribute(state.get_truncation_error(), "truncation_error", dsetName);
            h5ppFile.writeAttribute(status.chi_lim, "chi_lim", dsetName);
            h5ppFile.writeAttribute(status.chi_lim_max, "chi_lim_max", dsetName);
            save_log[dsetName] = save_point;
        }
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
        tools::common::profile::get_default_prof()["t_hdf"]->tic();
        h5ppFile.writeDataset(state.A_bare(), mps_prefix + "/M_A");
        h5ppFile.writeAttribute(state.A_bare().dimensions(), "dimensions", mps_prefix + "/M_A");
        h5ppFile.writeDataset(state.B(), mps_prefix + "/M_B");
        h5ppFile.writeAttribute(state.B().dimensions(), "dimensions", mps_prefix + "/M_B");
        tools::common::profile::get_default_prof()["t_hdf"]->toc();
        save_log[mps_prefix] = save_point;
    }


}

void tools::infinite::io::h5dset::save_model(h5pp::File &h5ppFile, const std::string &mpo_path, const StorageLevel &storage_level,
                                             const class_model_infinite &model) {
    if(storage_level < StorageLevel::FULL) return;
    // We do not expect the MPO's to change. Therefore if they exist, there is nothing else to do here
    if(h5ppFile.linkExists(mpo_path)) return tools::log->trace("The model has already been written to [{}]", mpo_path);
    tools::log->trace("Storing [{: ^6}]: mpo tensors", enum2str(storage_level));
    tools::common::profile::get_default_prof()["t_hdf"]->tic();
    model.get_mpo_siteA().save_mpo(h5ppFile, mpo_path);
    model.get_mpo_siteB().save_mpo(h5ppFile, mpo_path);
    h5ppFile.writeAttribute(2, "model_size", mpo_path);
    h5ppFile.writeAttribute(enum2str(settings::model::model_type), "model_type", mpo_path);
    tools::common::profile::get_default_prof()["t_hdf"]->toc();
}

void tools::infinite::io::h5dset::save_edges(h5pp::File &h5ppFile, const std::string &edges_prefix, const StorageLevel &storage_level,
                                             const class_edges_infinite &edges) {
    if(storage_level < StorageLevel::NORMAL) return;
    tools::common::profile::get_default_prof()["t_hdf"]->tic();
    const auto &ene = edges.get_ene_blk();
    const auto &var = edges.get_var_blk();
    h5ppFile.writeDataset(ene.L, edges_prefix + "/eneL");
    h5ppFile.writeDataset(ene.R, edges_prefix + "/eneR");
    h5ppFile.writeDataset(var.L, edges_prefix + "/varL");
    h5ppFile.writeDataset(var.R, edges_prefix + "/varR");
    tools::common::profile::get_default_prof()["t_hdf"]->toc();
}

//
// void tools::infinite::io::h5dset::write_hamiltonian_params(h5pp::File &h5pp_file, const std::string & sim_name, const StorageLevel & storage_level, const
// class_state_infinite &state){
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
// void tools::infinite::io::h5dset::write_all_measurements  (h5pp::File &h5pp_file, const std::string & sim_name, const StorageLevel & storage_level, const
// class_simulation_status & status, const class_state_infinite &state){
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
