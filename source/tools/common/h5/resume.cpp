
#include "config/enums.h"
#include "debug/exceptions.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include <h5pp/h5pp.h>

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

std::vector<std::string> tools::common::h5::resume::find_state_prefixes(const h5pp::File &h5file, AlgorithmType algo_type, std::string_view name) {
    std::string_view algo_name = enum2sv(algo_type);
    if(name.empty())
        tools::log->info("Searching for states for algorithm [{}] in file [{}]", algo_name, h5file.getFilePath());
    else
        tools::log->info("Searching for states matching [{}] for algorithm [{}] in file [{}]", name, algo_name, h5file.getFilePath());
    auto                     algo_prefix = fmt::format("/{}", algo_name);
    auto                     matches     = h5file.findGroups(name, algo_prefix, -1, 1);
    std::vector<std::string> state_prefixes;
    for(const auto &match : matches) {
        state_prefixes.emplace_back(fmt::format("{}/{}", algo_prefix, match));
        tools::log->debug(" -- found state prefix: {}", state_prefixes.back());
    }
    return state_prefixes;
}

std::vector<tools::common::h5::MpsInfo> tools::common::h5::resume::find_fully_stored_mps(const h5pp::File &h5file, std::string_view state_prefix) {
    std::vector<MpsInfo> mps_info;
    for(const auto &pfx : h5file.findGroups("mps", state_prefix)) {
        auto path  = fmt::format("{}/{}", state_prefix, pfx);
        auto iter  = h5file.readAttribute<size_t>(path, "iter");
        auto step  = h5file.readAttribute<size_t>(path, "step");
        auto event = h5file.readAttribute<StorageEvent>(path, "storage_event");
        auto level = h5file.readAttribute<StorageLevel>(path, "storage_level");
        if(level == StorageLevel::FULL) mps_info.emplace_back(MpsInfo{path, iter, step, event});
    }
    // Sort according to step
    std::sort(begin(mps_info), end(mps_info), [](const MpsInfo &i1, const MpsInfo &i2) { return i1.step < i2.step; });
    return mps_info;
}