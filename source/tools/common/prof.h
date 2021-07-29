#pragma once

#include <string_view>
class AlgorithmStatus;
enum class AlgorithmType;

namespace tools::common::profile {
    // Profiling

    //    extern void print_profiling_all();
    extern void print_profiling(const AlgorithmStatus &status);
    extern void print_profiling();
    //    extern void print_profiling_delta();
    //    extern void print_profiling_laps(std::optional<AlgorithmType> algo_type = std::nullopt);

    extern void print_mem_usage();

    extern double mem_usage_in_mb(std::string_view name);
    extern double mem_rss_in_mb();
    extern double mem_hwm_in_mb();
    extern double mem_vm_in_mb();

}