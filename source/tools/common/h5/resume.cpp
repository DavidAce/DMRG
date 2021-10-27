
#include <config/enums.h>
#include <h5pp/h5pp.h>
#include <tools/common/h5.h>
#include <tools/common/log.h>

std::string tools::common::h5::resume::extract_state_name(std::string_view state_prefix) {
    constexpr std::string_view state_pattern = "state_";
    std::string                state_name;
    std::string::size_type     start_pos = state_prefix.find(state_pattern);
    if(start_pos != std::string::npos) {
        std::string::size_type len = state_prefix.find('/', start_pos) - start_pos;
        state_name                 = state_prefix.substr(start_pos, len); // E.g. "/xDMRG/state_0/results/..." would match state_0
    }
    return state_name;
}

std::optional<size_t> tools::common::h5::resume::extract_state_number(std::string_view state_prefix) {
    std::string state_name = extract_state_name(state_prefix);
    std::string state_number;
    for(const auto &c : state_name) {
        if(std::isdigit(c)) state_number.push_back(c);
    }
    try {
        size_t number = std::stoul(state_number);
        return number;
    } catch(const std::exception &err) {
        tools::log->info("Could not convert {} to a number: {}", state_number, err.what());
        return std::nullopt;
    }
}

std::vector<std::string> tools::common::h5::resume::find_resumable_states(const h5pp::File &h5file, AlgorithmType algo_type, std::string_view search) {
    std::string_view         algo_name = enum2sv(algo_type);
    std::vector<std::string> state_prefix_candidates;

    if(search.empty())
        tools::log->info("Searching for resumable states from algorithm [{}] in file [{}]", algo_name, h5file.getFilePath());
    else
        tools::log->info("Searching for resumable states with keyword [{}] from algorithm [{}] in file [{}]", search, algo_name, h5file.getFilePath());

    for(const auto &candidate : h5file.getAttributeNames("common/storage_level")) {
        if(candidate.find(algo_name) != std::string::npos and h5file.readAttribute<std::string>(candidate, "common/storage_level") == "FULL")
            state_prefix_candidates.push_back(candidate);
    }

    tools::log->info("Found state candidates: {}", state_prefix_candidates);

    // Apply the search filter
    if(not search.empty()) {
        auto search_filter = [search](std::string_view x) {
            return x.find(search) == std::string::npos;
        };
        state_prefix_candidates.erase(std::remove_if(state_prefix_candidates.begin(), state_prefix_candidates.end(), search_filter),
                                      state_prefix_candidates.end());
        tools::log->info("States matching keyword [{}]:  {}", search, state_prefix_candidates);
    }

    // Return the results if done
    if(state_prefix_candidates.empty()) return {};
    if(state_prefix_candidates.size() == 1) return {state_prefix_candidates[0]};

    // Here we collect unique state names
    std::vector<std::string> state_names;
    for(const auto &pfx : state_prefix_candidates) {
        auto name = h5file.readAttribute<std::string>(pfx, "common/state_name");
        if(not name.empty()) {
            if(std::find(state_names.begin(), state_names.end(), name) == state_names.end()) { state_names.emplace_back(name); }
        }
    }

    // Sort state names in ascending order
    std::sort(state_names.begin(), state_names.end(), std::less<>());

    // For each state name, we find the one with highest step number
    std::vector<std::string> state_candidates_latest;
    for(const auto &name : state_names) {
        // Collect the candidates that have the current state name and their step
        std::vector<std::pair<size_t, std::string>> matching_candidates;
        for(const auto &candidate : state_prefix_candidates) {
            if(candidate.find(name) == std::string::npos) continue;
            auto steps = h5file.readAttribute<size_t>(candidate, "common/step");
            matching_candidates.emplace_back(std::make_pair(steps, candidate));
        }

        // Sort according to step number in decreasing order
        std::sort(matching_candidates.begin(), matching_candidates.end(), [](auto &lhs, auto &rhs) {
            if(lhs.first != rhs.first) return lhs.first > rhs.first;
            // Take care of cases with repeated step numbers
            else if(lhs.second.find("finished") != std::string::npos and rhs.second.find("finished") == std::string::npos)
                return true;
            else if(lhs.second.find("iter") != std::string::npos and rhs.second.find("iter") == std::string::npos)
                return true;
            return false;
        });
        for(const auto &candidate : matching_candidates) tools::log->info("Candidate step {} : [{}]", candidate.first, candidate.second);

        // Add the front candidate
        state_candidates_latest.emplace_back(matching_candidates.front().second);
    }
    return state_candidates_latest;
}
