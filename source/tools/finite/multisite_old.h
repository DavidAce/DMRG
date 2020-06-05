#pragma once

#include <general/nmspc_tensor_omp.h>
#include <list>

class class_state_finite;

namespace tools::finite::multisite{
    extern Eigen::DSizes<long,3> get_dimensions  (const class_state_finite &state, const std::vector<size_t> &active_sites);
    extern size_t                get_problem_size(const class_state_finite &state, const std::vector<size_t> &list_of_sites);
    extern std::vector<size_t>     generate_site_list(class_state_finite &state, const size_t threshold, const size_t max_sites, const size_t min_sites = 2);
    extern std::vector<size_t>     generate_truncated_site_list(class_state_finite &state, const size_t threshold, const size_t chi_lim, const size_t max_sites, const size_t min_sites = 2);


    std::vector<size_t>      activate_sites(const size_t threshold, const size_t max_sites, const size_t min_sites = 2);
    std::vector<size_t>      activate_truncated_sites(const long threshold, const size_t chi_lim, const size_t max_sites, const size_t min_sites = 2);
    Eigen::DSizes<long, 3> active_dimensions();
    size_t                 active_problem_size();

}
