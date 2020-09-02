#pragma once
#include <array>
#include <optional>
#include <vector>

struct state_measure_finite {
    std::optional<size_t>                length                        = std::nullopt;
    std::optional<long>                  bond_dimension_midchain       = std::nullopt;
    std::optional<long>                  bond_dimension_current        = std::nullopt;
    std::optional<std::vector<long>>     bond_dimensions               = std::nullopt;
    std::optional<double>                norm                          = std::nullopt;
    std::optional<std::array<double, 3>> spin_components               = std::nullopt;
    std::optional<std::vector<double>>   truncation_errors             = std::nullopt;
    std::optional<double>                entanglement_entropy_midchain = std::nullopt;
    std::optional<double>                entanglement_entropy_current  = std::nullopt;
    std::optional<std::vector<double>>   entanglement_entropies        = std::nullopt;
    std::optional<std::vector<double>>   renyi_2                       = std::nullopt;
    std::optional<std::vector<double>>   renyi_3                       = std::nullopt;
    std::optional<std::vector<double>>   renyi_4                       = std::nullopt;
    std::optional<std::vector<double>>   renyi_100                     = std::nullopt;
};
