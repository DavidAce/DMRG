#pragma once
#include <optional>
#include <array>
#include <vector>

struct tensors_measure_infinite {
    std::optional<size_t> length                            = std::nullopt;
    std::optional<double> energy_mpo                        = std::nullopt;
    std::optional<double> energy_per_site_mpo               = std::nullopt;
    std::optional<double> energy_variance_mpo               = std::nullopt;
    std::optional<double> energy_per_site_ham               = std::nullopt;
    std::optional<double> energy_per_site_mom               = std::nullopt;
    std::optional<double> energy_variance_per_site_mpo      = std::nullopt;
    std::optional<double> energy_variance_per_site_ham      = std::nullopt;
    std::optional<double> energy_variance_per_site_mom      = std::nullopt;
};
