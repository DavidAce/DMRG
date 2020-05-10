//
// Created by david on 2020-03-28.
//

#include <h5pp/h5pp.h>
#include <tools/common/io.h>

std::string tools::common::io::h5resume::find_resumable_state(const h5pp::File &h5ppFile, const std::string &sim_name, const std::string &search) {
    std::vector<std::string> prefix_candidates;
    for(const auto &candidate : h5ppFile.getAttributeNames("common/storage_level"))
        if(candidate.find(sim_name) != std::string::npos and h5ppFile.readAttribute<std::string>(candidate, "common/storage_level") == "FULL")
            prefix_candidates.push_back(candidate);

    // Apply the search filter
    auto search_filter = [search](std::string_view x) { return x.find(search) == std::string::npos; };
    prefix_candidates.erase(std::remove_if(prefix_candidates.begin(), prefix_candidates.end(), search_filter), prefix_candidates.end());

    // Return the results if done
    if(prefix_candidates.empty()) return "";
    if(prefix_candidates.size() == 1) return prefix_candidates[0];

    // We have a number of possible candidates coming from results, journal or projection. The most suitable one follows from a list of priorities,
    //  - Latest state: I.e. compare number in state_emin, state_emax, state_0, state_1....etc and take the latest
    //  - Latest step: I.e. compare the step-number in all the candidates (step number may be reset between states)

    // Check if there are multiple states. In that case, collect them and sort them
    std::vector<std::string> state_nums;
    for(const auto &candidate : prefix_candidates) {
        std::string::size_type pos = candidate.find("state_");
        if(pos != std::string::npos) state_nums.push_back(candidate.substr(pos, candidate.find('/', pos))); // E.g. "/xDMRG/state_0/..." would match [0]
    }

    // Erase all items that aren't numbers (i.e. filter away emin and emax)
    if(not state_nums.empty()) {
        auto has_letter_filter = [](std::string_view s) { return std::find_if(s.begin(), s.end(), [](unsigned char c) { return !std::isdigit(c); }) != s.end(); };
        state_nums.erase(std::remove_if(state_nums.begin(), state_nums.end(), has_letter_filter), state_nums.end());
    }

    // Erase all but the latest state from the candidate list
    if(not state_nums.empty()) {
        std::sort(state_nums.begin(), state_nums.end());
        auto state_filter = [state_nums](std::string_view x) { return x.find(state_nums.back()) == std::string::npos; };
        prefix_candidates.erase(std::remove_if(prefix_candidates.begin(), prefix_candidates.end(), state_filter), prefix_candidates.end());
    }

    // Return the results if done
    if(prefix_candidates.empty()) return "";
    if(prefix_candidates.size() == 1) return prefix_candidates[0];

    // By now we have narrowed it down such that we need to compare the step numbers of each candidate
    std::vector<std::pair<size_t, std::string>> step_sorted_candidates;
    for(const auto &candidate : prefix_candidates) {
        auto steps = h5ppFile.readAttribute<size_t>(candidate, "common/step");
        step_sorted_candidates.emplace_back(std::make_pair(steps, candidate));
    }

    // Sort according to step number in decreasing order
    std::sort(step_sorted_candidates.begin(), step_sorted_candidates.end(), [](auto &left, auto &right) { return left.first > right.first; });

    // Return the first element, which should contain the latest possible state
    return step_sorted_candidates.front().second;
}
