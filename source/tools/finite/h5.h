#pragma once
#include "general/sfinae.h"
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

struct StorageMeta;

namespace tools::finite::h5 {
    /* clang-format off */
    namespace save {
//        void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5file, const std::vector<std::string_view> &links);
//        void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5file, std::string_view link);
        template<typename T>
        void data_as_table(h5pp::File &h5file, std::string_view table_prefix, const AlgorithmStatus &status,const T * const data, size_t size, std::string_view table_name, std::string_view table_title, std::string_view fieldname);
        template<typename T>
        void data_as_table(h5pp::File &h5file, std::string_view table_prefix, const AlgorithmStatus &status,const T & data, std::string_view table_name, std::string_view table_title, std::string_view fieldname){
            if constexpr(sfinae::is_specialization_v<T, std::optional>){
                if(data.has_value()) data_as_table(h5file, table_prefix, status, data.value(), table_name, table_title, fieldname);
            }
            else if constexpr (sfinae::has_data_v<T> and sfinae::has_size_v<T>) data_as_table(h5file, table_prefix, status, data.data(), static_cast<size_t>(data.size()), table_name, table_title, fieldname);
            else if constexpr (std::is_arithmetic_v<T>) data_as_table(h5file, table_prefix, status, &data, 1, table_name, table_title, fieldname);
            else static_assert(sfinae::invalid_type_v<T> and "Datatype must have .data() and .size() (or be std::optional of such)");
        }

        extern int decide_layout(std::string_view prefix_path);
        template<typename T>
        extern void data      (h5pp::File & h5file, const T &data, std::string_view data_name, std::string_view state_name, const AlgorithmStatus &status, StorageReason storage_reason, std::optional<CopyPolicy> copy_policy);
        template<typename T>
        extern void data      (h5pp::File & h5file, const T &data, std::string_view data_name, std::string_view state_prefix, const StorageLevel & storage_level, const AlgorithmStatus & status);
        extern void state     (h5pp::File & h5file, std::string_view  state_prefix, const StorageLevel & storage_level, const StateFinite & state, const AlgorithmStatus & status);
        extern void model     (h5pp::File & h5file, std::string_view  model_prefix, const StorageLevel & storage_level, const ModelFinite & model);
        extern void mpo       (h5pp::File & h5file, std::string_view  model_prefix, const StorageLevel & storage_level, const ModelFinite & model);

        extern void correlations (h5pp::File & h5file, std::string_view  state_prefix, const StorageLevel & storage_level, const StateFinite & state, const AlgorithmStatus & status);

        extern void measurements        (h5pp::File & h5file, std::string_view  table_prefix, const StorageLevel & storage_level, const TensorsFinite & tensors, const AlgorithmStatus & status);
        extern void measurements        (h5pp::File & h5file, std::string_view  table_prefix, const StorageLevel & storage_level, const StateFinite & state, const ModelFinite & model, const EdgesFinite & edges, const AlgorithmStatus & status);
        extern void bond_dimensions     (h5pp::File & h5file, std::string_view  table_prefix, const StorageLevel & storage_level, const StateFinite & state, const AlgorithmStatus & status);
        extern void truncation_errors   (h5pp::File & h5file, std::string_view  table_prefix, const StorageLevel & storage_level, const StateFinite & state, const AlgorithmStatus & status);
        extern void entropies_neumann   (h5pp::File & h5file, std::string_view  table_prefix, const StorageLevel & storage_level, const StateFinite & state, const AlgorithmStatus & status);
        extern void entropies_renyi     (h5pp::File & h5file, std::string_view  table_prefix, const StorageLevel & storage_level, const StateFinite & state, const AlgorithmStatus & status);
        extern void entropies_number    (h5pp::File & h5file, std::string_view  table_prefix, const StorageLevel & storage_level, const StateFinite & state, const AlgorithmStatus & status);
        extern void expectations        (h5pp::File & h5file, std::string_view  table_prefix, const StorageLevel & storage_level, const StateFinite & state, const AlgorithmStatus & status);
        extern void structure_factors   (h5pp::File & h5file, std::string_view  table_prefix, const StorageLevel & storage_level, const StateFinite & state, const AlgorithmStatus & status);

        extern void setup_prefix(const AlgorithmStatus &status, const StorageReason &storage_reason, StorageLevel &storage_level,
                                 std::string_view state_name, std::string &state_prefix, std::string &model_prefix, std::string & timer_prefix, std::vector<std::string> &table_prefxs);
        extern void simulation(h5pp::File &h5file,
                                 const TensorsFinite & tensors,
                                 const AlgorithmStatus &status,
                                 const StorageReason &storage_reason,
                                 std::optional<CopyPolicy> copy_policy);
        extern void simulation(h5pp::File &h5file,
                                 const StateFinite & state,
                                 const ModelFinite & model,
                                 const EdgesFinite & edges,
                                 const AlgorithmStatus &status,
                                 const StorageReason &storage_reason,
                                 std::optional<CopyPolicy> copy_policy);

    }
    namespace load {
        extern void simulation (const h5pp::File & h5file, std::string_view  state_prefix, TensorsFinite & tensors, AlgorithmStatus & status, AlgorithmType algo_type);
        extern void state   (const h5pp::File & h5file, std::string_view  state_prefix, StateFinite & state, const AlgorithmStatus & status);
        extern void model   (const h5pp::File & h5file, std::string_view  state_prefix, ModelFinite & model);
        extern void validate (const h5pp::File & h5file, std::string_view  state_prefix, TensorsFinite & tensors, AlgorithmType algo_type);
    }

    namespace tmp{
        extern void copy(const AlgorithmStatus &status, const h5pp::File &h5file, StorageReason storage_reason, std::optional<CopyPolicy> copy_policy);

    }

//    namespace load{

//        extern std::vector<SimulationTask>
//            getTaskList(const h5pp::File &h5file, std::string_view sim_name, std::string_view state_prefix, const TensorsFinite & tensors, const class_algorithm_status & status);
//    }
    /* clang-format on */
}
