#pragma once
#include <array>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>
struct MeasurementsStateFinite {
    std::optional<size_t>                 length                           = std::nullopt;
    std::optional<long>                   bond_mid                         = std::nullopt;
    std::optional<long>                   bond_dim                         = std::nullopt;
    std::optional<std::vector<long>>      bond_dimensions                  = std::nullopt;
    std::optional<fp64>                   norm                             = std::nullopt;
    std::optional<std::array<fp64, 3>>    spin_components                  = std::nullopt;
    std::optional<std::vector<fp64>>      truncation_errors                = std::nullopt;
    std::optional<fp64>                   entanglement_entropy_midchain    = std::nullopt;
    std::optional<fp64>                   entanglement_entropy_current     = std::nullopt;
    std::optional<std::vector<fp64>>      entanglement_entropies           = std::nullopt;
    std::optional<fp64>                   number_entropy_midchain          = std::nullopt;
    std::optional<fp64>                   number_entropy_current           = std::nullopt;
    std::optional<std::vector<fp64>>      number_entropies                 = std::nullopt;
    std::optional<Eigen::Tensor<fp64, 2>> number_probabilities             = std::nullopt;
    std::optional<std::vector<fp64>>      renyi_2                          = std::nullopt;
    std::optional<std::vector<fp64>>      renyi_3                          = std::nullopt;
    std::optional<std::vector<fp64>>      renyi_4                          = std::nullopt;
    std::optional<std::vector<fp64>>      renyi_inf                        = std::nullopt;
    std::optional<Eigen::Tensor<fp64, 1>> expectation_values_sx            = std::nullopt;
    std::optional<Eigen::Tensor<fp64, 1>> expectation_values_sy            = std::nullopt;
    std::optional<Eigen::Tensor<fp64, 1>> expectation_values_sz            = std::nullopt;
    std::optional<Eigen::Tensor<fp64, 2>> correlation_matrix_sx            = std::nullopt;
    std::optional<Eigen::Tensor<fp64, 2>> correlation_matrix_sy            = std::nullopt;
    std::optional<Eigen::Tensor<fp64, 2>> correlation_matrix_sz            = std::nullopt;
    std::optional<double>                 structure_factor_x               = std::nullopt;
    std::optional<double>                 structure_factor_y               = std::nullopt;
    std::optional<double>                 structure_factor_z               = std::nullopt;
    std::optional<Eigen::Tensor<fp64, 1>> opdm_spectrum                    = std::nullopt;
    std::optional<Eigen::Tensor<cx64, 2>> opdm                             = std::nullopt;
    std::optional<Eigen::ArrayXXd>        subsystem_entanglement_entropies = std::nullopt;
    std::optional<Eigen::ArrayXXd>        information_lattice              = std::nullopt;
    std::optional<Eigen::ArrayXd>         information_per_scale            = std::nullopt;
    std::optional<double>                 information_center_of_mass       = std::nullopt;
    std::optional<double>                 see_time                         = std::nullopt; /*! The time it took to calculate the last subsystem_entanglement_entropies */
};
