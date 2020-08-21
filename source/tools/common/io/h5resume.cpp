//
// Created by david on 2020-03-28.
//

#include <config/enums.h>
#include <h5pp/h5pp.h>
#include <tools/common/io.h>
#include <tools/common/log.h>

std::string tools::common::io::h5resume::extract_state_name(const std::string &state_prefix) {
    std::string_view       state_pattern = "state_";
    std::string            state_name;
    std::string::size_type start_pos = state_prefix.find(state_pattern);
    if(start_pos != std::string::npos) {
        std::string::size_type len = state_prefix.find('/', start_pos) - start_pos;
        state_name                 = state_prefix.substr(start_pos, len); // E.g. "/xDMRG/state_0/results/..." would match state_0
    }
    return state_name;
}

std::optional<size_t> tools::common::io::h5resume::extract_state_number(const std::string &state_prefix) {
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

std::string tools::common::io::h5resume::find_resumable_state(const h5pp::File &h5ppFile, AlgorithmType algo_type, const std::string &search) {
    std::string_view         algo_name = enum2str(algo_type);
    std::vector<std::string> state_prefix_candidates;

    if(search.empty()) tools::log->info("Searching for resumable states from algorithm [{}] in file [{}]", algo_name, h5ppFile.getFilePath());
    else
        tools::log->info("Searching for resumable states with keyword [{}] from algorithm [{}] in file [{}]", search, algo_name, h5ppFile.getFilePath());

    for(const auto &candidate : h5ppFile.getAttributeNames("common/storage_level"))
        if(candidate.find(algo_name) != std::string::npos and h5ppFile.readAttribute<std::string>(candidate, "common/storage_level") == "FULL")
            state_prefix_candidates.push_back(candidate);

    tools::log->info("Found state candidates: {}", state_prefix_candidates);


    // Apply the search filter
    if(not search.empty()){
        auto search_filter = [search](std::string_view x) { return x.find(search) == std::string::npos; };
        state_prefix_candidates.erase(std::remove_if(state_prefix_candidates.begin(), state_prefix_candidates.end(), search_filter), state_prefix_candidates.end());
        tools::log->info("States matching keyword [{}]:  {}", search, state_prefix_candidates);
    }

    // Return the results if done
    if(state_prefix_candidates.empty()) return "";
    if(state_prefix_candidates.size() == 1) return state_prefix_candidates[0];

    // We have a number of possible candidates coming from results, journal or projection. The most suitable one follows from a list of priorities,
    //  - Latest state: I.e. compare number in state_0, state_1....etc and take the latest
    //  - Latest step: I.e. compare the step-number in all the candidates (step number may be reset between states)

    // Here we collect the state numbers
    std::vector<std::string> state_nums;
    for(const auto &candidate : state_prefix_candidates) {
        // E.g. "/xDMRG/state_0/..." would extract "0"
        auto num = extract_state_number(candidate);
        if(num) state_nums.emplace_back(std::to_string(num.value()));
    }

    // Here we sort the state numbers, then erase all but the latest state from the candidate list
    if(not state_nums.empty()) {
        std::sort(state_nums.begin(), state_nums.end());
        auto state_filter = [state_nums](std::string_view x) { return x.find(state_nums.back()) == std::string::npos; };
        state_prefix_candidates.erase(std::remove_if(state_prefix_candidates.begin(), state_prefix_candidates.end(), state_filter),
                                      state_prefix_candidates.end());
    }

    // Return the results if done
    if(state_prefix_candidates.empty()) return "";
    if(state_prefix_candidates.size() == 1) return state_prefix_candidates[0];

    // By now we have narrowed it down such that we need to compare the step numbers of each candidate
    std::vector<std::pair<size_t, std::string>> step_sorted_candidates;
    for(const auto &candidate : state_prefix_candidates) {
        auto steps = h5ppFile.readAttribute<size_t>(candidate, "common/step");
        step_sorted_candidates.emplace_back(std::make_pair(steps, candidate));
    }

    // Sort according to step number in decreasing order
    std::sort(step_sorted_candidates.begin(), step_sorted_candidates.end(), [](auto &left, auto &right) { return left.first > right.first; });

    // Return the first element, which should contain the latest possible state
    return step_sorted_candidates.front().second;
}
