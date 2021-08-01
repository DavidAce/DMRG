#pragma once
#include <string>
class StateFinite;
class ModelFinite;
class EdgesFinite;
class TensorsFinite;
class AlgorithmStatus;
enum class StorageLevel;
enum class StorageReason;
enum class CopyPolicy;
enum class AlgorithmType;
namespace h5pp {
    class File;
}

namespace tools::finite::h5 {
    struct SavePack {
        h5pp::File                                                  &file;
        std::string                                                  target_prefix;
        std::optional<std::reference_wrapper<const StateFinite>>     state;
        std::optional<std::reference_wrapper<const ModelFinite>>     model;
        std::optional<std::reference_wrapper<const EdgesFinite>>     edges;
        std::optional<std::reference_wrapper<const AlgorithmStatus>> status;
        StorageLevel                                                 storage_level;
        StorageReason                                                storage_reason;
        std::optional<CopyPolicy>                                    copy_policy;
        SavePack(h5pp::File &file) : file(file) {}
    };
    struct LoadPack {
        const h5pp::File                                      &file;
        std::string                                            target_prefix;
        std::optional<std::reference_wrapper<StateFinite>>     state;
        std::optional<std::reference_wrapper<ModelFinite>>     model;
        std::optional<std::reference_wrapper<EdgesFinite>>     edges;
        std::optional<std::reference_wrapper<AlgorithmStatus>> status;
        StorageLevel                                           storage_level;
        StorageReason                                          storage_reason;
        LoadPack(const h5pp::File &file) : file(file) {}
    };

    /* clang-format off */
    namespace save {
        void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile, const std::vector<std::string_view> &links);
//        void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5ppFile, std::string_view link);

        extern int decide_layout(std::string_view prefix_path);
        template<typename T>
        extern void data      (h5pp::File & h5ppFile, const T &data, std::string_view data_name, std::string_view state_name, const AlgorithmStatus &status, StorageReason storage_reason, std::optional<CopyPolicy> copy_policy);
        extern void state     (h5pp::File & h5ppFile, std::string_view  state_prefix, const StorageLevel & storage_level, const StateFinite & state, const AlgorithmStatus & status);
        extern void model     (h5pp::File & h5ppFile, std::string_view  model_prefix, const StorageLevel & storage_level, const ModelFinite & model);
        extern void mpo       (h5pp::File & h5ppFile, std::string_view  model_prefix, const StorageLevel & storage_level, const ModelFinite & model);

        extern void entgm  (h5pp::File & h5ppFile, std::string_view  state_prefix, const StorageLevel & storage_level, const StateFinite & state, const AlgorithmStatus & status);
        extern void measurements   (h5pp::File & h5ppFile, std::string_view  table_prefix, const StorageLevel & storage_level, const TensorsFinite & tensors, const AlgorithmStatus & status, AlgorithmType algo_type);
        extern void measurements   (h5pp::File & h5ppFile, std::string_view  table_prefix, const StorageLevel & storage_level, const StateFinite & state, const ModelFinite & model, const EdgesFinite & edges, const AlgorithmStatus & status, AlgorithmType algo_type);
        extern void setup_prefix(const AlgorithmStatus &status, const StorageReason &storage_reason, StorageLevel &storage_level,
                                 std::string_view state_name, std::string &state_prefix, std::string &model_prefix, std::vector<std::string> &table_prefxs);
        extern void simulation(h5pp::File &h5ppfile,
                                 const TensorsFinite & tensors,
                                 const AlgorithmStatus &status,
                                 const StorageReason &storage_reason,
                                 std::optional<CopyPolicy> copy_policy);
        extern void simulation(h5pp::File &h5ppfile,
                                 const StateFinite & state,
                                 const ModelFinite & model,
                                 const EdgesFinite & edges,
                                 const AlgorithmStatus &status,
                                 const StorageReason &storage_reason,
                                 std::optional<CopyPolicy> copy_policy);

    }
    namespace load {
        extern void simulation (const h5pp::File & h5ppFile, std::string_view  state_prefix, TensorsFinite & tensors, AlgorithmStatus & status, AlgorithmType algo_type);
        extern void state   (const h5pp::File & h5ppFile, std::string_view  state_prefix, StateFinite & state, const AlgorithmStatus & status);
        extern void model   (const h5pp::File & h5ppFile, std::string_view  state_prefix, ModelFinite & model);
        extern void validate (const h5pp::File & h5ppFile, std::string_view  state_prefix, TensorsFinite & tensors, AlgorithmType algo_type);
    }

    namespace tmp{
        extern void copy(const AlgorithmStatus &status, const h5pp::File &h5ppfile, StorageReason storage_reason, std::optional<CopyPolicy> copy_policy);

    }

//    namespace load{

//        extern std::vector<SimulationTask>
//            getTaskList(const h5pp::File &h5ppFile, std::string_view sim_name, std::string_view state_prefix, const TensorsFinite & tensors, const class_algorithm_status & status);
//    }
    /* clang-format on */
}
