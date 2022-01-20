#pragma once
#include <array>
#include <optional>
#include <vector>

struct MeasurementsStateFinite {
    std::optional<size_t>                   length                        = std::nullopt;
    std::optional<long>                     bond_dimension_midchain       = std::nullopt;
    std::optional<long>                     bond_dimension_current        = std::nullopt;
    std::optional<std::vector<long>>        bond_dimensions               = std::nullopt;
    std::optional<double>                   norm                          = std::nullopt;
    std::optional<std::array<double, 3>>    spin_components               = std::nullopt;
    std::optional<std::vector<double>>      truncation_errors             = std::nullopt;
    std::optional<double>                   entanglement_entropy_midchain = std::nullopt;
    std::optional<double>                   entanglement_entropy_current  = std::nullopt;
    std::optional<std::vector<double>>      entanglement_entropies        = std::nullopt;
    std::optional<double>                   number_entropy_midchain       = std::nullopt;
    std::optional<double>                   number_entropy_current        = std::nullopt;
    std::optional<std::vector<double>>      number_entropies              = std::nullopt;
    std::optional<std::vector<double>>      renyi_2                       = std::nullopt;
    std::optional<std::vector<double>>      renyi_3                       = std::nullopt;
    std::optional<std::vector<double>>      renyi_4                       = std::nullopt;
    std::optional<std::vector<double>>      renyi_inf                     = std::nullopt;
    std::optional<Eigen::Tensor<double, 1>> expectation_values_sx         = std::nullopt;
    std::optional<Eigen::Tensor<double, 1>> expectation_values_sy         = std::nullopt;
    std::optional<Eigen::Tensor<double, 1>> expectation_values_sz         = std::nullopt;
    std::optional<Eigen::Tensor<double, 2>> correlation_matrix_sx         = std::nullopt;
    std::optional<Eigen::Tensor<double, 2>> correlation_matrix_sy         = std::nullopt;
    std::optional<Eigen::Tensor<double, 2>> correlation_matrix_sz         = std::nullopt;
    std::optional<double>                   structure_factor_x            = std::nullopt;
    std::optional<double>                   structure_factor_y            = std::nullopt;
    std::optional<double>                   structure_factor_z            = std::nullopt;
};
