#pragma once
#include <optional>
#include <array>
#include <vector>

struct state_measure_infinite {
    std::optional<double> norm                         = std::nullopt;
    std::optional<long>   bond_dimension               = std::nullopt;
    std::optional<double> entanglement_entropy         = std::nullopt;
    std::optional<double> truncation_error             = std::nullopt;
};