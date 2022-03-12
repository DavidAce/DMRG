#pragma once

#include <array>
#include <optional>
#include <vector>

class StateFinite;
class ModelFinite;
class EdgesFinite;

namespace tools::finite::multisite {
    extern long                get_problem_size(const StateFinite &state, std::optional<std::vector<size_t>> active_sites);
    extern std::array<long, 3> get_dimensions(const StateFinite &state, std::optional<std::vector<size_t>> sites = std::nullopt);
    extern std::array<long, 4> get_dimensions(const ModelFinite &model, std::optional<std::vector<size_t>> sites = std::nullopt);
    extern std::array<long, 4> get_dimensions_squared(const ModelFinite &model, std::optional<std::vector<size_t>> active_sites = std::nullopt);
    extern std::vector<size_t> generate_site_list(StateFinite &state, long threshold, size_t max_sites, size_t min_sites = 1, const std::string &tag = "");
    extern std::vector<size_t> generate_truncated_site_list(StateFinite &state, long threshold, long bond_limit, size_t max_sites, size_t min_sites = 1);
    extern std::vector<size_t> activate_sites(long threshold, size_t max_sites, size_t min_sites = 1);
    extern std::vector<size_t> activate_truncated_sites(long threshold, long bond_limit, size_t max_sites, size_t min_sites = 1);
    extern std::array<long, 3> active_dimensions();
    extern long                active_problem_size();
}
