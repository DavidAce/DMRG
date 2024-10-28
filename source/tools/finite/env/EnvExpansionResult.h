#pragma once
#include <array>
#include <functional>
#include <optional>
#include <string>
#include <unsupported/Eigen/CXX11/Tensor>
struct EnvExpansionResult {
    bool                                   ok       = false; /*!< True if the expansion took place (false on queries, e.g., envExpandMode == NONE) */
    long                                   posL     = -1;    /*!< Position left of the expanded bond */
    long                                   posR     = -1;    /*!< Position right of the expanded bond */
    std::array<long, 3>                    dimL_old = {};    /*!< Dimensions of the left site before the expansion */
    std::array<long, 3>                    dimL_new = {};    /*!< Dimensions of the left site after expanding and truncating */
    std::array<long, 3>                    dimR_old = {};    /*!< Dimensions of the right site before the expansion */
    std::array<long, 3>                    dimR_new = {};    /*!< Dimensions of the right site after expanding and truncating */
    std::vector<size_t>                    sites;
    std::vector<std::array<long, 3>>       dims_old, dims_new;
    std::vector<long>                      bond_old, bond_new;
    double                                 ene_old   = std::numeric_limits<double>::quiet_NaN(); /*!< The old expectation value  of energy */
    double                                 ene_new   = std::numeric_limits<double>::quiet_NaN(); /*!< The old expectation value  of energy */
    double                                 var_old   = std::numeric_limits<double>::quiet_NaN(); /*!< The new expectation value  of energy variance*/
    double                                 var_new   = std::numeric_limits<double>::quiet_NaN(); /*!< The new expectation value  of energy variance*/
    double                                 alpha_mps = std::numeric_limits<double>::quiet_NaN(); /*!< The mixing factor used in the expansion */
    double                                 alpha_h1v = std::numeric_limits<double>::quiet_NaN(); /*!< The mixing factor for the H-term in the expansion */
    double                                 alpha_h2v = std::numeric_limits<double>::quiet_NaN(); /*!< The mixing factor for the H²-term in the expansion */
    Eigen::Tensor<std::complex<double>, 3> mixed_blk;                                            /*!< The mixed mps block */
    std::optional<std::reference_wrapper<EnvEne>> env_ene = std::nullopt;                        /*!< Reference to the environment expanded (either L or R) */
    std::optional<std::reference_wrapper<EnvVar>> env_var = std::nullopt;                        /*!< Reference to the environment expanded (either L or R) */
    std::string                                   msg;
};