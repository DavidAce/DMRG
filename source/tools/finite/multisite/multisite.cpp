#include "../multisite.h"
#include "debug/exceptions.h"
#include <math/num.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/site/mpo/MpoSite.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateFinite.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>

std::array<long, 3> tools::finite::multisite::get_dimensions(const StateFinite &state, std::optional<std::vector<size_t>> sites) {
    if(not sites) sites = state.active_sites;
    if(sites.value().empty()) return std::array<long, 3>{0, 0, 0};
    std::array<long, 3> dimensions{};
    std::sort(sites.value().begin(), sites.value().end());
    if(sites.value().front() > sites.value().back()) throw except::logic_error("Given site list is not increasing: {}", sites.value());

    dimensions[1] = state.get_mps_site(sites.value().front()).get_M().dimension(1);
    dimensions[2] = state.get_mps_site(sites.value().back()).get_M().dimension(2);
    dimensions[0] = 1;
    for(auto &site : sites.value()) { dimensions[0] *= state.get_mps_site(site).get_M().dimension(0); }
    return dimensions;
}

std::array<long, 4> tools::finite::multisite::get_dimensions(const ModelFinite &model, std::optional<std::vector<size_t>> sites) {
    if(not sites) sites = model.active_sites;
    if(sites.value().empty()) return std::array<long, 4>{0, 0, 0, 0};
    if(sites.value().front() > sites.value().back()) throw except::logic_error("Given site list is not increasing: {}", sites.value());
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

std::array<long, 4> tools::finite::multisite::get_dimensions_squared(const ModelFinite &model, std::optional<std::vector<size_t>> active_sites) {
    if(not active_sites) active_sites = model.active_sites;
    if(active_sites.value().empty()) return std::array<long, 4>{0, 0, 0, 0};
    if(active_sites.value().front() > active_sites.value().back()) throw except::logic_error("Active site list is not increasing: {}", active_sites.value());
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

long tools::finite::multisite::get_problem_size(const StateFinite &state, std::optional<std::vector<size_t>> active_sites) {
    auto dims = get_dimensions(state, std::move(active_sites));
    return (dims[0] * dims[1] * dims[2]);
}

std::vector<size_t> tools::finite::multisite::generate_site_list(StateFinite &state, long threshold, size_t max_sites, const size_t min_sites,
                                                                 const std::string &tag) {
    if(max_sites < min_sites) throw std::runtime_error("generate site list: asked for max sites < min sites");
    tools::log->trace("Multisite activation {} | site {} | direction {} | sites min {} max {} | max problem size {}", tag, state.get_position<long>(),
                      state.get_direction(), min_sites, max_sites, threshold);

    const auto initial_position = state.get_position<long>();
    auto       direction        = state.get_direction();
    //    long                                position         = initial_position;
    long                             length  = state.get_length<long>();
    bool                             at_edge = initial_position <= -1 or initial_position >= length;
    std::vector<size_t>              sites;
    std::vector<long>                sizes;
    std::vector<std::array<long, 3>> shape;

    //    long max_pos = std::clamp(initial_position , initial_position, std::min<long>(length-1, static_cast<const long>(max_sites - 1)));
    //    long min_pos = std::clamp(initial_position, initial_position, std::min<long>(length-1, static_cast<const long>(max_sites - 1)));

    if(not at_edge) {
        long max_pos = initial_position;
        long min_pos = initial_position;
        if(direction > 0) {
            max_pos = std::clamp<long>(initial_position + static_cast<long>(max_sites), initial_position, length - 1);
            min_pos = initial_position;
        } else {
            max_pos = std::clamp<long>(initial_position + 1, initial_position, length - 1);
            min_pos = std::clamp<long>(max_pos - static_cast<long>(max_sites) + 1, 0, initial_position);
        }

        auto range = num::range<size_t>(min_pos, max_pos + 1); // +1 to include last position
        if(direction < 0) std::reverse(range.begin(), range.end());

        sites.emplace_back(initial_position); // Current position is always included
        sizes.emplace_back(get_problem_size(state, sites));
        shape.emplace_back(get_dimensions(state, sites));
        for(auto &pos : range) {
            if(std::find(sites.begin(), sites.end(), pos) != sites.end()) continue; // Skip the first site
            sites.emplace_back(pos);
            sizes.emplace_back(get_problem_size(state, sites));
            shape.emplace_back(get_dimensions(state, sites));
        }
    }
    tools::log->trace("Candidate sites {}", sites);
    tools::log->trace("Candidate sizes {}", sizes);
    // Evaluate best cost. Threshold depends on optSolver
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
            throw except::logic_error("At least two sites required!");
        } else if(sites.empty()) {
            throw except::logic_error("No sites for a jump");
        } else {
            sites.pop_back();
            sizes.pop_back();
            shape.pop_back();
        }
    }

    std::sort(sites.begin(), sites.end());
    if(at_edge or shape.empty() or sizes.empty())
        tools::log->trace("Multisite activation: site {} | direction {} | sites min {} max {} | max problem size {} | chosen sites {} | reason {}",
                          initial_position, direction, min_sites, max_sites, threshold, sites, reason);
    else
        tools::log->trace("Multisite activation: site {} | direction {} | sites min {} max {} | max problem size {} | chosen sites {} | "
                          "shape {} = {} | reason {}",
                          initial_position, direction, min_sites, max_sites, threshold, sites, shape.back(), sizes.back(), reason);

    if(not at_edge and sites.size() < min_sites) throw except::logic_error("Activated sites ({}) < min_sites ({})", sites.size(), min_sites);
    if(not at_edge and sites.size() > max_sites) throw except::logic_error("Activated sites ({}) > max_sites ({})", sites.size(), max_sites);
    return sites;
}
