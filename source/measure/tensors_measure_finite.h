#pragma once
#include <optional>
#include <array>
#include <vector>

struct tensors_measure_finite {
    std::optional<size_t>                length                        = std::nullopt;
    std::optional<double>                energy                        = std::nullopt;
    std::optional<double>                energy_per_site               = std::nullopt;
    std::optional<double>                energy_variance               = std::nullopt;
    std::optional<double>                energy_variance_per_site      = std::nullopt;
};
