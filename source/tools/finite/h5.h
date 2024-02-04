#pragma once
#include "general/sfinae.h"
#include <string>
class StateFinite;
class ModelFinite;
class EdgesFinite;
class TensorsFinite;
class AlgorithmStatus;
enum class StorageLevel;
enum class StorageEvent;
enum class StoragePolicy;
enum class CopyPolicy;
enum class AlgorithmType;
namespace h5pp {
    class File;
    namespace hid {
        class h5t;
    }
}
namespace tools::common::h5 {
    namespace save {
        bool should_save(const StorageInfo &sinfo, StoragePolicy policy);
    }
    struct MpsInfo;
}

struct StorageInfo;

namespace tools::finite::h5 {

    /* clang-format off */
    namespace save {
//        void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5file, const std::vector<std::string_view> &links);
//        void bootstrap_save_log(std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> &save_log, const h5pp::File &h5file, std::string_view link);
        template<typename T>
        void data_as_table(h5pp::File &h5file, const StorageInfo & sinfo, const T * const data, size_t size, std::string_view table_name, std::string_view table_title, std::string_view fieldname);
        template<typename T>
        void data_as_table(h5pp::File &h5file, const StorageInfo & sinfo, const T & data, std::string_view table_name, std::string_view table_title, std::string_view fieldname){
            if constexpr(sfinae::is_std_optional_v<T>){
                if(data.has_value()) data_as_table(h5file, sinfo, data.value(), table_name, table_title, fieldname);
            }
            else if constexpr (sfinae::has_data_v<T> and sfinae::has_size_v<T>) data_as_table(h5file, sinfo, data.data(), static_cast<size_t>(data.size()), table_name, table_title, fieldname);
            else if constexpr (std::is_arithmetic_v<T>) data_as_table(h5file, sinfo, &data, 1, table_name, table_title, fieldname);
            else static_assert(sfinae::invalid_type_v<T> and "Datatype must have .data() and .size() (or be std::optional of such)");
        }
        template<typename T>
        void data_as_table_vla(h5pp::File &h5file, const  StorageInfo & sinfo, const std::vector<T> & data, const h5pp::hid::h5t & h5elem_t, std::string_view table_name, std::string_view table_title, std::string_view fieldname);
        extern int decide_layout(std::string_view prefix_path);
        template<typename T>
        extern void data      (h5pp::File & h5file, const StorageInfo & sinfo,  const T &data, std::string_view data_name, CopyPolicy copy_policy);
        template<typename T>
        extern void data      (h5pp::File & h5file, const StorageInfo & sinfo, const T &data, std::string_view data_name, std::string_view prefix);
        extern void bonds     (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);
        extern void state     (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);
        extern void model     (h5pp::File & h5file, const StorageInfo & sinfo, const ModelFinite & model);
        extern void mpo       (h5pp::File & h5file, const StorageInfo & sinfo, const ModelFinite & model);

        extern void correlations (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);

        extern void measurements                    (h5pp::File & h5file, const StorageInfo & sinfo, const TensorsFinite & tensors, const AlgorithmStatus & status);
        extern void measurements                    (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state, const ModelFinite & model, const EdgesFinite & edges, const AlgorithmStatus & status);
        extern void bond_dimensions                 (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);
        extern void schmidt_values                  (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);
        extern void truncation_errors               (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);
        extern void entanglement_entropies          (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);
        extern void subsystem_entanglement_entropies(h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);
        extern void renyi_entropies                 (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);
        extern void number_entropies                (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);
        extern void expectations                    (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);
        extern void structure_factors               (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);
        extern void kvornings_marker                (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);
        extern void number_probabilities            (h5pp::File & h5file, const StorageInfo & sinfo, const StateFinite & state);

        [[nodiscard]] extern StorageInfo get_storage_info(const StateFinite & state, const AlgorithmStatus &status);

        extern void simulation(h5pp::File &h5file,
                               const TensorsFinite & tensors,
                               const AlgorithmStatus &status,
                               CopyPolicy copy_policy);
        extern void simulation(h5pp::File &h5file,
                               const StateFinite & state,
                               const ModelFinite & model,
                               const EdgesFinite & edges,
                               const AlgorithmStatus &status,
                               CopyPolicy copy_policy);

    }
    namespace find{
        extern std::string find_unitary_circuit(const h5pp::File &h5file, AlgorithmType algo_type, std::string_view name,
                                                                  size_t iter);
    }
    namespace load {
        using MpsInfo = tools::common::h5::MpsInfo;
        extern void simulation (const h5pp::File & h5file, std::string_view  state_prefix, TensorsFinite & tensors, AlgorithmStatus & status, AlgorithmType algo_type);
        extern void state   (const h5pp::File & h5file, std::string_view  state_prefix, StateFinite & state, MpsInfo & mpsinfo);
        extern void model   (const h5pp::File & h5file, AlgorithmType algo_type, ModelFinite & model);
        extern void validate (const h5pp::File & h5file, std::string_view  state_prefix, TensorsFinite & tensors, AlgorithmStatus & status, AlgorithmType algo_type);
    }

    namespace tmp{
        extern void copy(const AlgorithmStatus &status, const h5pp::File &h5file, StorageEvent storage_reason, std::optional<CopyPolicy> copy_policy);

    }

//    namespace load{

//        extern std::vector<SimulationTask>
//            getTaskList(const h5pp::File &h5file, std::string_view sim_name, std::string_view state_prefix, const TensorsFinite & tensors, const class_algorithm_status & status);
//    }
    /* clang-format on */
}
