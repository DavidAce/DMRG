//
// Created by david on 2019-06-24.
//

#include <tensors/model/class_model_finite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/multisite.h>

std::array<long, 3> tools::finite::multisite::get_dimensions(const class_state_finite &state, std::optional<std::vector<size_t>> sites) {
    if(not sites) sites = state.active_sites;
    if(sites.value().empty()) return std::array<long, 3>{0, 0, 0};
    std::array<long, 3> dimensions{};
    std::sort(sites.value().begin(), sites.value().end());
    if(sites.value().front() > sites.value().back())
        throw std::runtime_error(fmt::format("Given site list is not increasing: {}", sites.value()));

    dimensions[1] = state.get_mps_site(sites.value().front()).get_M().dimension(1);
    dimensions[2] = state.get_mps_site(sites.value().back()).get_M().dimension(2);
    dimensions[0] = 1;
    for(auto &site : sites.value()) { dimensions[0] *= state.get_mps_site(site).get_M().dimension(0); }
    return dimensions;
}

std::array<long, 4> tools::finite::multisite::get_dimensions(const class_model_finite &model, std::optional<std::vector<size_t>> sites) {
    if(not sites) sites = model.active_sites;
    if(sites.value().empty()) return std::array<long, 4>{0, 0, 0, 0};
    if(sites.value().front() > sites.value().back())
        throw std::runtime_error(fmt::format("Given site list is not increasing: {}", sites.value()));
    std::array<long, 4> dimensions{};
    dimensions[0] = model.get_mpo(sites.value().front()).MPO().dimension(0);
    dimensions[1] = model.get_mpo(sites.value().back()).MPO().dimension(1);
    dimensions[2] = 1;
    dimensions[3] = 1;
    for(auto &site : sites.value()) {
        dimensions[2] *= model.get_mpo(site).MPO().dimension(2);
        dimensions[3] *= model.get_mpo(site).MPO().dimension(3);
    }
    return dimensions;
}

std::array<long, 4> tools::finite::multisite::get_dimensions_squared(const class_model_finite &model, std::optional<std::vector<size_t>> active_sites) {
    if(not active_sites) active_sites = model.active_sites;
    if(active_sites.value().empty()) return std::array<long, 4>{0, 0, 0, 0};
    if(active_sites.value().front() > active_sites.value().back())
        throw std::runtime_error(fmt::format("Active site list is not increasing: {}", active_sites.value()));
    std::array<long, 4> dimensions{};
    dimensions[0] = model.get_mpo(active_sites.value().front()).MPO2().dimension(0);
    dimensions[1] = model.get_mpo(active_sites.value().back()).MPO2().dimension(1);
    dimensions[2] = 1;
    dimensions[3] = 1;
    for(auto &site : active_sites.value()) {
        dimensions[2] *= model.get_mpo(site).MPO2().dimension(2);
        dimensions[3] *= model.get_mpo(site).MPO2().dimension(3);
    }
    return dimensions;
}

long tools::finite::multisite::get_problem_size(const class_state_finite &state, std::optional<std::vector<size_t>> active_sites) {
    auto dims = get_dimensions(state, std::move(active_sites));
    return (dims[0] * dims[1] * dims[2]);
}

std::vector<size_t> tools::finite::multisite::generate_site_list(class_state_finite &state, long threshold, size_t max_sites, const size_t min_sites) {
    if(max_sites < min_sites) throw std::runtime_error("generate site list: asked for max sites < min sites");
    tools::log->trace("Multisite activation: site {} | direction {} | sites min {} max {} | max problem size {}", state.get_position<long>(),
                      state.get_direction(), min_sites, max_sites, threshold);

    const auto                          initial_position = state.get_position<long>();
    auto                                direction        = state.get_direction();
    long                                position         = initial_position;
    long                                length           = state.get_length<long>();
    bool                                at_edge          = position <= -1 or position >= length;
    std::vector<long>                   sizes;
    std::vector<size_t>                 sites;
    std::vector<std::array<long, 3>> shape;
    while(true) {
        if(position >= 0) {
            sites.emplace_back(position);
            sizes.emplace_back(get_problem_size(state, sites));
            shape.emplace_back(get_dimensions(state, sites));
            if(sites.size() >= max_sites) break;
        }
        position += direction;
        if(position <= -1 or position >= length) break;
    }
    tools::log->trace("Candidate sites {}", sites);
    tools::log->trace("Candidate sizes {}", sizes);
    // Evaluate best cost. Threshold depends on optSpace
    // Case 1: All costs are equal              -> take all sites
    // Case 2: Costs increase indefinitely      -> take until threshold
    // Case 3: Costs increase and saturate      -> take until threshold

    std::string reason;
    while(true) {
        if(sites.empty() and at_edge) {
            reason = "Can't select sites beyond the edge";
            break;
        }
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
    if(at_edge or shape.empty() or sizes.empty())
        tools::log->debug("Multisite activation: site {} | direction {} | sites min {} max {} | max problem size {} | chosen sites {} | reason {}",
                          initial_position, direction, min_sites, max_sites, threshold, sites, reason);
    else
        tools::log->debug(
            "Multisite activation: site {} | direction {} | sites min {} max {} | max problem size {} | chosen sites {} | shape {} = {} | reason {}",
            initial_position, direction, min_sites, max_sites, threshold, sites, shape.back(), sizes.back(), reason);

    if(not at_edge and sites.size() < min_sites) throw std::runtime_error(fmt::format("Activated sites ({}) < min_sites ({})", sites.size(), min_sites));
    if(not at_edge and sites.size() > max_sites) throw std::runtime_error(fmt::format("Activated sites ({}) > max_sites ({})", sites.size(), max_sites));
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
