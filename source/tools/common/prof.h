#pragma once

#include <general/class_tic_toc.h>
#include <map>
#include <memory>
#include <optional>
#include <vector>
#include <tools/common/fmt.h>

class class_tic_toc;
enum class AlgorithmType;
namespace tools::common::profile {
    // Profiling
    inline std::unique_ptr<class_tic_toc> t_tot;          // + Total Time
    namespace internal {
        template<typename KeyT, typename ValT>
        class insert_ordered_map {
            private:
            std::vector<std::pair<KeyT, ValT>> data;

            public:
            const ValT &operator[](const KeyT &key) const {
                auto it = find(key);
                if(it == data.end()){
                    if constexpr(std::is_convertible_v<KeyT,std::string>)
                        throw std::runtime_error(fmt::format("Invalid key: {}",key));
                    else throw std::runtime_error("Invalid key");
                }
                return it->second;
            }
            ValT &operator[](const KeyT &key) {
                auto it = find(key);
                if(it == data.end()) {
                    data.emplace_back(std::make_pair(key, ValT()));
                    return data.back().second;
                } else
                    return it->second;
            }

            [[nodiscard]] auto begin() { return data.begin(); }
            [[nodiscard]] auto end() { return data.end(); }
            [[nodiscard]] auto begin() const { return data.begin(); }
            [[nodiscard]] auto end() const { return data.end(); }
            [[nodiscard]] auto find(const KeyT & key){
                return std::find_if(data.begin(), data.end(), [&key](const auto &element) { return element.first == key; });
            }
            [[nodiscard]] auto find(const KeyT & key) const {
                return std::find_if(data.begin(), data.end(), [&key](const auto &element) { return element.first == key; });
            }
        };
        using MapTicTocUnique   = insert_ordered_map<std::string, std::unique_ptr<class_tic_toc>>;
        using MapTicTocShared   = insert_ordered_map<std::string, std::shared_ptr<class_tic_toc>>;
        using AlgoProf    = insert_ordered_map<AlgorithmType, MapTicTocUnique>;
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