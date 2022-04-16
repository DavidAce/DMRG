#pragma once
#include <array>
#include <optional>
#include <vector>

struct MeasurementsTensorsFinite {
    std::optional<size_t> length                    = std::nullopt;
    std::optional<double> energy                    = std::nullopt;
    std::optional<double> energy_variance           = std::nullopt;
    std::optional<double> energy_shift              = std::nullopt;
    std::optional<double> energy_minus_energy_shift = std::nullopt;
};
