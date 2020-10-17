#pragma once

#include <general/nmspc_tensor_extra.h>
#include <optional>

class class_state_finite;
class class_model_finite;
class class_edges_finite;

namespace tools::finite::multisite {
    extern Eigen::DSizes<long, 3> get_dimensions(const class_state_finite &state, std::optional<std::vector<size_t>> active_sites = std::nullopt);
    extern Eigen::DSizes<long, 4> get_dimensions(const class_model_finite &model, std::optional<std::vector<size_t>> active_sites = std::nullopt);
    extern Eigen::DSizes<long, 4> get_dimensions_squared(const class_model_finite &model, std::optional<std::vector<size_t>> active_sites = std::nullopt);
    extern long                   get_problem_size(const class_state_finite &state, std::optional<std::vector<size_t>> active_sites);
    extern std::vector<size_t>    generate_site_list(class_state_finite &state, long threshold, size_t max_sites, size_t min_sites = 2);
    extern std::vector<size_t>    generate_truncated_site_list(class_state_finite &state, long threshold, long chi_lim, size_t max_sites, size_t min_sites = 2);
    extern std::vector<size_t>    activate_sites(long threshold, size_t max_sites, size_t min_sites = 2);
    extern std::vector<size_t>    activate_truncated_sites(long threshold, long chi_lim, size_t max_sites, size_t min_sites = 2);
    extern Eigen::DSizes<long, 3> active_dimensions();
    extern long                   active_problem_size();
}
