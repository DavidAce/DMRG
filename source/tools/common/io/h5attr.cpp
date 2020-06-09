
#include <h5pp/h5pp.h>
#include <algorithms/class_algorithm_status.h>
#include <config/enums.h>
#include <string>
#include <tools/common/io.h>
#include <tools/common/log.h>
void tools::common::io::h5attr::write_meta(h5pp::File &h5ppFile, const std::string &algo_name, const std::string & state_name, const std::string &state_prefix, const std::string &model_prefix,
                                           ModelType model_type, const StorageLevel &storage_level, const class_algorithm_status &status) {

    // Write metadata into /common so that we can find state paths later.
    // In particular, one dataset should contain the root paths to all states



    std::string storage_level_str(enum2str(storage_level));
    std::string model_type_str   (enum2str(model_type));
    std::string state_root        = algo_name + "/" + state_name;
    std::string ham_prefix        = model_prefix + "/Hamiltonian";
    std::string mpo_prefix        = model_prefix + "/mpo";
    std::string mps_prefix        = state_prefix + "/mps";

    // Register this state prefix as the lastest written if its a checkpoint or a final result
    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", state_root, storage_level_str, algo_name);
    h5ppFile.writeAttribute(status.algorithm_has_finished,state_root, algo_name);

    // Add the storage level of all states as an attribute of this simulation group.
    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", state_prefix, storage_level_str, algo_name);
    h5ppFile.writeAttribute(storage_level_str, state_prefix, algo_name);

    // Write storage level into an attribute of the current state that is being written to
    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", "storage_level", storage_level_str, state_prefix);
    h5ppFile.writeAttribute(storage_level_str, "storage_level", state_prefix);

    // Collect the storage level of all states that have been written on any simulation
    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", "common/storage_level", state_prefix, storage_level_str);
    h5ppFile.writeDataset("The attributes on this dataset are the storage level of each state", "common/storage_level");
    h5ppFile.writeAttribute(storage_level_str, state_prefix, "common/storage_level");

    // Write the iteration and step number of this state (useful when resuming to compare state freshness)
    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", "iteration", status.iter, state_prefix);
    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", "step", status.step, state_prefix);
    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", "position", status.position, state_prefix);
    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", "model_type", model_type_str, state_prefix);
    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", "ham_prefix", ham_prefix, state_prefix);
    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", "mpo_prefix", mpo_prefix, state_prefix);
    h5ppFile.writeAttribute(status.iter, "iteration", state_prefix);
    h5ppFile.writeAttribute(status.step, "step", state_prefix);
    h5ppFile.writeAttribute(status.position, "position", state_prefix);
    h5ppFile.writeAttribute(model_type_str, "model_type", state_prefix);
    h5ppFile.writeAttribute(ham_prefix, "ham_prefix", state_prefix);
    h5ppFile.writeAttribute(mpo_prefix, "mpo_prefix", state_prefix);

    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", state_prefix, status.iter, "common/state_roots");
    h5ppFile.writeDataset("The attributes on this dataset contains the root paths to all states", "common/state_roots");
    h5ppFile.writeAttribute(status.algorithm_has_finished, state_root, "common/state_roots");


    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", state_prefix, status.iter, "common/iteration");
    h5ppFile.writeDataset("The attributes on this dataset are the current step of each state", "common/iteration");
    h5ppFile.writeAttribute(status.iter, state_prefix, "common/iteration");

    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", state_prefix, status.step, "common/step");
    h5ppFile.writeDataset("The attributes on this dataset are the current step of each state", "common/step");
    h5ppFile.writeAttribute(status.step, state_prefix, "common/step");

    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", state_prefix, ham_prefix, "common/ham_prefix");
    h5ppFile.writeDataset("The attributes on this dataset are paths to the Hamiltonian Table", "common/ham_prefix");
    h5ppFile.writeAttribute(ham_prefix, state_prefix, "common/ham_prefix");

    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", state_prefix, mpo_prefix, "common/mpo_prefix");
    h5ppFile.writeDataset("The attributes on this dataset are paths to MPOs", "common/mpo_prefix");
    h5ppFile.writeAttribute(mpo_prefix, state_prefix, "common/mpo_prefix");

    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", state_prefix, mps_prefix, "common/mps_prefix");
    h5ppFile.writeDataset("The attributes on this dataset are paths to MPSs", "common/mps_prefix");
    h5ppFile.writeAttribute(mps_prefix, state_prefix, "common/mps_prefix");

    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", state_prefix, status.algorithm_has_succeeded, "common/success");
    h5ppFile.writeDataset("The attributes on this dataset tell if the simulation succeeded", "common/success");
    h5ppFile.writeAttribute(status.algorithm_has_succeeded, state_prefix, "common/success");

    tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", state_prefix, status.algorithm_has_finished, "common/finish");
    h5ppFile.writeDataset("The attributes on this dataset tell if the simulation finished", "common/finish");
    h5ppFile.writeAttribute(status.algorithm_has_finished, state_prefix, "common/finish");
}
