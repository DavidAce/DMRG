#pragma once

#include <general/nmspc_tensor_omp.h>
#include <list>
#include <optional>

class class_state_finite;
class class_model_finite;
class class_edges_finite;

namespace tools::finite::multisite {
    extern Eigen::DSizes<long, 3> get_dimensions(const class_state_finite &state, std::optional<std::list<size_t>> active_sites = std::nullopt);
    extern Eigen::DSizes<long, 4> get_dimensions(const class_model_finite &model, std::optional<std::list<size_t>> active_sites = std::nullopt);
    extern size_t                 get_problem_size(const class_state_finite &state, std::optional<std::list<size_t>> active_sites);
    extern std::list<size_t>      generate_site_list(class_state_finite &state, size_t threshold, size_t max_sites, size_t min_sites = 2);
    extern std::list<size_t> generate_truncated_site_list(class_state_finite &state, size_t threshold, size_t chi_lim, size_t max_sites, size_t min_sites = 2);
    extern std::list<size_t> activate_sites(size_t threshold, size_t max_sites, size_t min_sites = 2);
    extern std::list<size_t> activate_truncated_sites(long threshold, size_t chi_lim, size_t max_sites, size_t min_sites = 2);
    extern Eigen::DSizes<long, 3> active_dimensions();
    extern size_t                 active_problem_size();
}
