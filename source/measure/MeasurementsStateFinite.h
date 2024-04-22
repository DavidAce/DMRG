#pragma once
#include <array>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>
struct MeasurementsStateFinite {
    std::optional<size_t>                   length                           = std::nullopt;
    std::optional<long>                     bond_mid                         = std::nullopt;
    std::optional<long>                     bond_dim                         = std::nullopt;
    std::optional<std::vector<long>>        bond_dimensions                  = std::nullopt;
    std::optional<real>                     norm                             = std::nullopt;
    std::optional<std::array<real, 3>>      spin_components                  = std::nullopt;
    std::optional<std::vector<real>>        truncation_errors                = std::nullopt;
    std::optional<real>                     entanglement_entropy_midchain    = std::nullopt;
    std::optional<real>                     entanglement_entropy_current     = std::nullopt;
    std::optional<std::vector<real>>        entanglement_entropies           = std::nullopt;
    std::optional<real>                     number_entropy_midchain          = std::nullopt;
    std::optional<real>                     number_entropy_current           = std::nullopt;
    std::optional<std::vector<real>>        number_entropies                 = std::nullopt;
    std::optional<Eigen::Tensor<real, 2>>   number_probabilities             = std::nullopt;
    std::optional<std::vector<real>>        renyi_2                          = std::nullopt;
    std::optional<std::vector<real>>        renyi_3                          = std::nullopt;
    std::optional<std::vector<real>>        renyi_4                          = std::nullopt;
    std::optional<std::vector<real>>        renyi_inf                        = std::nullopt;
    std::optional<Eigen::Tensor<real, 1>>   expectation_values_sx            = std::nullopt;
    std::optional<Eigen::Tensor<real, 1>>   expectation_values_sy            = std::nullopt;
    std::optional<Eigen::Tensor<real, 1>>   expectation_values_sz            = std::nullopt;
    std::optional<Eigen::Tensor<real, 2>>   correlation_matrix_sx            = std::nullopt;
    std::optional<Eigen::Tensor<real, 2>>   correlation_matrix_sy            = std::nullopt;
    std::optional<Eigen::Tensor<real, 2>>   correlation_matrix_sz            = std::nullopt;
    std::optional<double>                   structure_factor_x               = std::nullopt;
    std::optional<double>                   structure_factor_y               = std::nullopt;
    std::optional<double>                   structure_factor_z               = std::nullopt;
    std::optional<Eigen::Tensor<real, 1>>   opdm_spectrum                 = std::nullopt;
    std::optional<Eigen::Tensor<cplx, 2>>   opdm                             = std::nullopt;
    std::optional<Eigen::ArrayXXd>          subsystem_entanglement_entropies = std::nullopt;
};
