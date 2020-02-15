#pragma once

#include <general/nmspc_tensor_omp.h>
#include <list>

class class_state_finite;

namespace tools::finite::multisite{
    extern Eigen::DSizes<long,3> get_dimensions  (const class_state_finite &state, const std::list<size_t> &list_of_sites);
    extern size_t                get_problem_size(const class_state_finite &state, const std::list<size_t> &list_of_sites);
    extern std::list<size_t>     generate_site_list(class_state_finite &state, const size_t threshold, const size_t max_sites, const size_t min_sites = 2);
    extern std::list<size_t>     generate_truncated_site_list(class_state_finite &state, const size_t threshold, const size_t chi_lim, const size_t max_sites, const size_t min_sites = 2);
}
