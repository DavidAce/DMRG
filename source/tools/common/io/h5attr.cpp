
#include <algorithms/class_algorithm_status.h>
#include <config/enums.h>
#include <h5pp/h5pp.h>
#include <string>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

namespace tools::common::io::h5attr {
    template<typename AttrType>
    void save_attr(h5pp::File &file, const AttrType &attrData, const std::string &attrName, const std::string &dsetPath, const std::string &dsetText) {
        tools::log->trace("Attribute -- {: <40} = {: <40} on dset: {}", attrName, attrData, dsetPath);
        if(not file.linkExists(dsetPath)) file.writeDataset(dsetText, dsetPath);
        file.writeAttribute(attrData, attrName, dsetPath);
    }
}

namespace tools::common::io::h5attr {
    void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile,
                            const std::string &state_prefix) {
        if(save_log.empty()) {
            try {
                uint64_t step = 0;
                uint64_t iter = 0;
                if(h5ppFile.linkExists("common/step")) step = h5ppFile.readAttribute<uint64_t>(state_prefix, "common/step");
                if(h5ppFile.linkExists("common/iteration")) iter = h5ppFile.readAttribute<uint64_t>(state_prefix, "common/iteration");
                save_log[state_prefix] = std::make_pair(iter, step);

            } catch(const std::exception &ex) { tools::log->warn("Could not bootstrap save_log for {}: {}", state_prefix, ex.what()); }
        }
    }
}

void tools::common::io::h5attr::save_meta(h5pp::File &h5ppFile, const StorageLevel &storage_level, const StorageReason &storage_reason,
                                          const ModelType &model_type, size_t model_size, const AlgorithmType &algo_type, const std::string &state_name,
                                          const std::string &state_prefix, const std::string &model_prefix, const class_algorithm_status &status) {
    // Write metadata into /common so that we can find state paths later.
    // Under /common we have the following dummy datasets defining attributes which map state_prefix to useful information:
    //
    // -- finished       | true or 1 if state_prefix is finished
    // -- storage_level  | state_prefix -> storage_level
    // -- storage_reason | state_prefix -> storage_reason
    // -- state_root     | state_prefix -> state_root
    // -- hamiltonian    | state_prefix -> path to hamiltonian table
    // -- mps_prefix     | state_prefix -> mps_prefix
    // -- mpo_prefix     | state_prefix -> mpo_prefix
    // -- model_type     | state_prefix -> model_type
    // -- model_size     | state_prefix -> model_size
    // -- algo_type      | state_prefix -> algo_type
    // -- state_name     | state_prefix -> state_name
    // -- iteration      | state_prefix -> iteration
    // -- step           | state_prefix -> step
    // -- position       | state_prefix -> position of the mps

    // Checks if the current entries have already been written
    auto t_hdf = tools::common::profile::get_default_prof()["t_hdf"]->tic_token();

    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    bootstrap_save_log(save_log, h5ppFile, state_prefix);

    auto save_point = std::make_pair(status.iter, status.step);
    if(save_log[state_prefix] == save_point) return;


    std::string storage_level_str(enum2str(storage_level));
    std::string storage_reason_str(enum2str(storage_reason));
    std::string model_name(enum2str(model_type));
    std::string algo_name(enum2str(algo_type));
    std::string state_root  = algo_name + '/' + state_name;
    std::string hamiltonian = model_prefix + "/hamiltonian";
    std::string mpo_prefix  = model_prefix + "/mpo";
    std::string mps_prefix  = state_prefix + "/mps";

    save_attr(h5ppFile, status.algorithm_has_finished, state_prefix, "common/finished", "Maps state_prefix -> finished");
    save_attr(h5ppFile, storage_level_str, state_prefix, "common/storage_level", "Maps state_prefix -> storage_level");
    save_attr(h5ppFile, storage_reason_str, state_prefix, "common/storage_reason", "Maps state_prefix -> storage_reason");
    save_attr(h5ppFile, state_root, state_prefix, "common/state_root", "Maps state_prefix -> state_root");
    save_attr(h5ppFile, model_prefix, state_prefix, "common/model_prefix", "Maps state_prefix -> model_prefix");
    save_attr(h5ppFile, hamiltonian, state_prefix, "common/hamiltonian", "Maps state_prefix -> hamiltonian");
    save_attr(h5ppFile, mpo_prefix, state_prefix, "common/mpo_prefix", "Maps state_prefix -> mpo_prefix");
    save_attr(h5ppFile, mps_prefix, state_prefix, "common/mps_prefix", "Maps state_prefix -> mps_prefix");
    save_attr(h5ppFile, model_name, state_prefix, "common/model_type", "Maps state_prefix -> model_type");
    save_attr(h5ppFile, model_size, state_prefix, "common/model_size", "Maps state_prefix -> model_size");
    save_attr(h5ppFile, algo_name, state_prefix, "common/algo_type", "Maps state_prefix -> algo_type");
    save_attr(h5ppFile, state_name, state_prefix, "common/state_name", "Maps state_prefix -> state_name");
    save_attr(h5ppFile, status.iter, state_prefix, "common/iteration", "Maps state_prefix -> iteration");
    save_attr(h5ppFile, status.step, state_prefix, "common/step", "Maps state_prefix -> step");
    save_attr(h5ppFile, status.position, state_prefix, "common/position", "Maps state_prefix -> position");
    save_log[state_prefix] = save_point;
}
