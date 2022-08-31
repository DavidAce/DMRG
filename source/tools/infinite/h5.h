#pragma once

#include <string>
#include <vector>
class TensorsInfinite;
class StateInfinite;
class ModelInfinite;
class EdgesInfinite;
class AlgorithmStatus;
enum class StorageLevel;
enum class SimulationTask;
namespace h5pp {
    class File;
}

namespace tools::infinite::h5 {
    /* clang-format off */
    namespace save{
        extern int decide_layout(std::string_view prefix_path);
        extern void bonds            (h5pp::File & h5file, std::string_view  state_prefix, const StorageLevel & storage_level, const StateInfinite & state, const AlgorithmStatus &status);
        extern void state            (h5pp::File & h5file, std::string_view  state_prefix, const StorageLevel & storage_level, const StateInfinite & state, const AlgorithmStatus &status);
        extern void edges            (h5pp::File & h5file, std::string_view  edges_prefix, const StorageLevel & storage_level, const EdgesInfinite & edges);
        extern void model            (h5pp::File & h5file, std::string_view  model_prefix, const StorageLevel & storage_level, const ModelInfinite & model);
        extern void mpo              (h5pp::File & h5file, std::string_view  model_prefix, const StorageLevel & storage_level, const ModelInfinite & model);
        extern void measurements     (h5pp::File & h5file, std::string_view  table_prefix, const StorageLevel & storage_level, const TensorsInfinite & tensors, const AlgorithmStatus & status);
    }

    namespace load{
        extern void tensors (const h5pp::File & h5file, std::string_view  state_prefix, TensorsInfinite & state, AlgorithmStatus & status);
        extern void state (const h5pp::File & h5file, std::string_view  state_prefix, StateInfinite & state, const AlgorithmStatus & status);
        extern void model (const h5pp::File & h5file, std::string_view  state_prefix, ModelInfinite & state, const AlgorithmStatus & status);
        extern void validate (const h5pp::File & h5file, std::string_view  state_prefix, const TensorsInfinite & state, const AlgorithmStatus & status);
        extern std::vector<SimulationTask>
            getTaskList(const h5pp::File &h5file, std::string_view state_prefix, const StateInfinite & state, const AlgorithmStatus & status);

    }
    /* clang-format on */
}
