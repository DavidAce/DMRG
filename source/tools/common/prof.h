#pragma once

#include <memory>
#include <optional>
#include <vector>
#include <general/class_tic_toc.h>

class class_tic_toc;
enum class AlgorithmType;
namespace tools::common::profile {
    // Profiling
    inline std::unique_ptr<class_tic_toc> t_tot;          // + Total Time

    namespace internal {
        // Implement a custom ordered map
        // In the map, we want keep the order in which they were appended. This behavior does not exist in stl ordered maps.
        template<typename KeyT, typename ValT>
        class insert_ordered_map {
            private:
            std::vector<std::pair<KeyT, ValT>> data;
            public:
            using iterator = typename std::vector<std::pair<KeyT, ValT>>::iterator;
            using const_iterator = typename std::vector<std::pair<KeyT, ValT>>::const_iterator;
            ValT &operator[](const KeyT &key);
            void append(const KeyT & key, ValT val);
            [[nodiscard]] iterator begin();
            [[nodiscard]] iterator end();
            [[nodiscard]] const_iterator begin() const;
            [[nodiscard]] const_iterator end() const;
            [[nodiscard]] iterator find(const KeyT & key);
            [[nodiscard]] const_iterator find(const KeyT & key) const;
        };
        using MapTicTocUnique   = insert_ordered_map<std::string, std::unique_ptr<class_tic_toc>>;
        using AlgoProf          = insert_ordered_map<AlgorithmType, MapTicTocUnique>;
        inline std::optional<AlgorithmType> default_algo_type = std::nullopt;
    }

    inline internal::AlgoProf prof;
    extern internal::MapTicTocUnique & get_default_prof();
    extern void set_default_prof(AlgorithmType algo_type);
    AlgorithmType get_current_algo_type();

    extern void print_profiling_all();
    extern void print_profiling(std::optional<AlgorithmType> algo_type = std::nullopt);
    extern void print_profiling_delta();
    extern void print_profiling_laps(std::optional<AlgorithmType> algo_type = std::nullopt);
    extern void init_profiling();
    //    extern void reset_profiling();
    extern void reset_profiling(std::optional<AlgorithmType> algo_type = std::nullopt, const std::vector<std::string> & excl = {});
//    extern void reset_for_run_algorithm(std::optional<AlgorithmType> algo_type = std::nullopt, const std::vector<std::string> & excl = {"t_pre","t_pos"});

    extern void print_mem_usage();

    extern double mem_usage_in_mb(std::string_view name);
    extern double mem_rss_in_mb();
    extern double mem_hwm_in_mb();
    extern double mem_vm_in_mb();

}