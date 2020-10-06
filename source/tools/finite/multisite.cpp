//
// Created by david on 2019-06-24.
//

#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/multisite.h>

Eigen::DSizes<long, 3> tools::finite::multisite::get_dimensions(const class_state_finite &state, std::optional<std::vector<size_t>> active_sites) {
    if(not active_sites) active_sites = state.active_sites;
    if(active_sites.value().empty()) return Eigen::DSizes<long, 3>{0, 0, 0};
    Eigen::DSizes<long, 3> dimensions;
    std::sort(active_sites.value().begin(), active_sites.value().end());
    if(active_sites.value().front() > active_sites.value().back())
        throw std::runtime_error(fmt::format("Active site list is not increasing: {}", active_sites.value()));

    dimensions[1] = state.get_mps_site(active_sites.value().front()).get_M().dimension(1);
    dimensions[2] = state.get_mps_site(active_sites.value().back()).get_M().dimension(2);
    dimensions[0] = 1;
    for(auto &site : active_sites.value()) { dimensions[0] *= state.get_mps_site(site).get_M().dimension(0); }
    return dimensions;
}

Eigen::DSizes<long, 4> tools::finite::multisite::get_dimensions(const class_model_finite &model, std::optional<std::vector<size_t>> active_sites) {
    if(not active_sites) active_sites = model.active_sites;
    if(active_sites.value().empty()) return Eigen::DSizes<long, 4>{0, 0, 0, 0};
    if(active_sites.value().front() > active_sites.value().back())
        throw std::runtime_error(fmt::format("Active site list is not increasing: {}", active_sites.value()));
    Eigen::DSizes<long, 4> dimensions;
    dimensions[0] = model.get_mpo(active_sites.value().front()).MPO().dimension(0);
    dimensions[1] = model.get_mpo(active_sites.value().back()).MPO().dimension(1);
    dimensions[2] = 1;
    dimensions[3] = 1;
    for(auto &site : active_sites.value()) {
        dimensions[2] *= model.get_mpo(site).MPO().dimension(2);
        dimensions[3] *= model.get_mpo(site).MPO().dimension(3);
    }
    return dimensions;
}

long tools::finite::multisite::get_problem_size(const class_state_finite &state, std::optional<std::vector<size_t>> active_sites) {
    auto dims = get_dimensions(state, std::move(active_sites));
    return (dims[0] * dims[1] * dims[2]);
}

std::vector<size_t> tools::finite::multisite::generate_site_list(class_state_finite &state, long threshold, size_t max_sites, const size_t min_sites) {
    if(max_sites < min_sites) throw std::runtime_error("generate site list: asked for max sites < min sites");
    tools::log->trace("Multisite activation: site {} | direction {} | sites min {} max {} | max problem size {}", state.get_position(),
                      state.get_direction(),min_sites, max_sites, threshold);
    using namespace Textra;
    int                                 direction = state.get_direction();
    int                                 position  = static_cast<int>(state.get_position());
    int                                 length    = static_cast<int>(state.get_length());
    std::vector<long>                   sizes;
    std::vector<size_t>                 sites;
    std::vector<Eigen::DSizes<long, 3>> shape;
    if(direction == -1) position = std::min(position + 1, length - 1); // If going to the left, take position to be the site on the right of the center bond.
    while(true) {
        sites.emplace_back(position);
        sizes.emplace_back(get_problem_size(state, sites));
        shape.emplace_back(get_dimensions(state,sites));
        position += direction;
        if(sites.size() >= max_sites) break;
        if(position == -1 or position == length) break;
    }
    tools::log->trace("Candidate sites {}", sites);
    tools::log->trace("Candidate sizes {}", sizes);
    // Evaluate best cost. Threshold depends on optSpace
    // Case 1: All costs are equal              -> take all sites
    // Case 2: Costs increase indefinitely      -> take until threshold
    // Case 3: Costs increase and saturate      -> take until threshold

    std::string reason;
    while(true) {
        bool allequal = std::all_of(sizes.begin(), sizes.end(), [sizes](long c) { return c == sizes.front(); });
        auto size     = sizes.back();
        if(size < threshold and sites.size() == max_sites) {
            reason = "reached max sites";
            break;
        }
        if(size <= threshold and sites.size() <= max_sites) {
            reason = fmt::format("good problem shape found");
            break;
        } else if(sites.size() <= min_sites) {
            reason = fmt::format("at least {} sites had to be kept", min_sites);
            break;
        } else if(allequal and sites.size() <= max_sites) {
            reason = fmt::format("problem sizes are equal: {}", size);
            break;
        } else if(sites.size() == 1) {
            throw std::logic_error("At least two sites required!");
        } else if(sites.empty()) {
            throw std::logic_error("No sites for a jump");
        } else {
            sites.pop_back();
            sizes.pop_back();
            shape.pop_back();
        }
    }
    std::sort(sites.begin(), sites.end());
    tools::log->debug(
        "Multisite activation: site {} | direction {} | sites min {} max {} | max problem size {} | chosen sites {} | shape {} = {} | reason {}", state.get_position(),
                      state.get_direction(),min_sites, max_sites, threshold, sites, shape.back(), sizes.back(),reason);
    if(sites.size() < 2) throw std::runtime_error("Less than 2 active sites");
    return sites;
}

std::vector<size_t> tools::finite::multisite::generate_truncated_site_list(class_state_finite &state, long threshold, long chi_lim, const size_t max_sites,
                                                                           const size_t min_sites) {
    auto active_sites = generate_site_list(state, threshold, max_sites, min_sites);
    if(active_sites.size() == max_sites) return active_sites; // Fastest outcome
    state.active_sites.clear();

    // Consider that the active site list may be limited by the edge
    if(state.get_direction() > 0 and active_sites.back() == state.get_length() - 1) return active_sites;
    if(state.get_direction() < 0 and active_sites.front() == 0) return active_sites;

    // Try activating more sites by truncating some sites ahead.
    size_t best_num_sites = min_sites;
    for(size_t num_sites = min_sites; num_sites <= max_sites; num_sites++) {
        auto state_copy = state;
        tools::finite::mps::truncate_next_sites(state_copy, chi_lim, num_sites);
        auto new_active_sites = generate_site_list(state_copy, threshold, num_sites, min_sites);
        if(new_active_sites.size() > active_sites.size()) {
            best_num_sites = num_sites;
            active_sites   = new_active_sites;
        }
    }
    tools::finite::mps::truncate_next_sites(state, chi_lim, best_num_sites);
    auto resulting_sites = generate_site_list(state, threshold, max_sites, min_sites);
    tools::log->info("Problem size {} | active sites: {}", get_problem_size(state, resulting_sites), resulting_sites);
    tools::log->info("Bond dimensions {}", tools::finite::measure::bond_dimensions(state));
    return resulting_sites;
}

//
// using namespace Textra;
// using Scalar = class_state_finite::Scalar;

//
// double tools::finite::measure::multisite::internal::significant_digits(double H2, double E2){
//    double max_digits    = std::numeric_limits<double>::max_digits10;
//    double lost_bits     = -std::log2(1.0 - std::abs(std::min(H2,E2)/std::max(H2,E2)));
//    double lost_digits   = std::log10(std::pow(2.0,lost_bits));
////    tools::log->trace("Significant digits: {}",std::floor(max_digits - lost_digits));
//    return digits = std::floor(max_digits - lost_digits);
//}
