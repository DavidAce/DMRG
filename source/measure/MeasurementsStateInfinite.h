#pragma once
#include <array>
#include <optional>
#include <vector>

struct MeasurementsStateInfinite {
    std::optional<double> norm                 = std::nullopt;
    std::optional<long>   bond_dimension       = std::nullopt;
    std::optional<double> entanglement_entropy = std::nullopt;
    std::optional<double> truncation_error     = std::nullopt;
};