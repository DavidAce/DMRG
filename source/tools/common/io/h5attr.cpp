
#include <h5pp/h5pp.h>
#include <simulation/enums.h>
#include <string>
#include <tools/common/io.h>
void tools::common::io::h5attr::write_prefix_meta(h5pp::File &h5ppFile, const std::string &prefix, const std::string &sim_name, const std::string &sim_tag, const std::string & model_type,
                                                  const StorageLevel &storage_level, size_t iteration, size_t step, size_t position) {
    std::string storage_level_str = enum2str(storage_level);
    std::string model_path = sim_name + "/model/Hamiltonian";
    // Add the storage level of all prefixes as an attribute of this simulation group.
    h5ppFile.writeAttribute(storage_level_str, prefix, sim_name);

    // Write storage level into an attribute of the current prefix that is being written to
    h5ppFile.writeAttribute(storage_level_str, "storage_level", prefix);

    // Collect the storage level of all states that have been written on any simulation
    std::string help_storage = "The attributes on this dataset specify the storage level of each written state";
    h5ppFile.writeDataset(help_storage, "common/storage_level");
    h5ppFile.writeAttribute(storage_level_str, prefix, "common/storage_level");

    // Write the iteration and step number of this prefix (useful when resuming to compare state freshness)
    h5ppFile.writeAttribute(iteration, "iteration", prefix);
    h5ppFile.writeAttribute(step, "step", prefix);
    h5ppFile.writeAttribute(position, "position", prefix);
    h5ppFile.writeAttribute(model_type, "model_type", prefix);
    h5ppFile.writeAttribute(model_path, "model_path", prefix);

    // Collect the iteration and step number of all states that have been written on any simulation
    std::string help_iter = "The attributes on this dataset specify the iteration number reached on each state";
    h5ppFile.writeDataset(help_iter, "common/iteration");
    h5ppFile.writeAttribute(iteration, prefix, "common/iteration");

    std::string help_step = "The attributes on this dataset specify the step number reached on each state";
    h5ppFile.writeDataset(help_step, "common/step");
    h5ppFile.writeAttribute(step, prefix, "common/step");

    std::string help_position = "The attributes on this dataset specify the current position on the chain for each state";
    h5ppFile.writeDataset(help_step, "common/position");
    h5ppFile.writeAttribute(position, prefix, "common/position");

    std::string help_model_type = "The attributes on this dataset specify the model type used on the simulation";
    h5ppFile.writeDataset(help_model_type, "common/model_type");
    h5ppFile.writeAttribute(model_type, prefix, "common/model_type");

    std::string help_model_path = "The attributes on this dataset specify path to the model Hamiltonian table used on the simulation";
    h5ppFile.writeDataset(help_model_path, "common/model_path");
    h5ppFile.writeAttribute(model_path, prefix, "common/model_path");
}
