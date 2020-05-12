#pragma once
#include <string>
#include <complex>
#include <optional>
#include <tools/finite/opt-internals/enum_classes.h>
enum class ModelType;

/* clang-format off */

class class_state_finite;

namespace tools::finite::mps {
    using Scalar = std::complex<double>;

    extern void initialize                          (class_state_finite & state, ModelType model_type, size_t num_sites, size_t position);
    extern void normalize                           (class_state_finite & state, std::optional<size_t> chi_lim = std::nullopt,std::optional<double> svd_threshold = std::nullopt);
    extern void random_product_state                (class_state_finite & state, const std::string & parity_sector, const long state_number, const bool use_pauli_eigenstates = false);
    extern void random_current_state                (class_state_finite & state, const std::string & parity_sector1, const std::string & parity_sector2);
    extern void rebuild_environments                (class_state_finite & state);
    extern void move_center_point                   (class_state_finite & state, std::optional<size_t> chi_lim = std::nullopt, std::optional<double> svd_threshold = std::nullopt); /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */
    extern void truncate_all_sites                  (class_state_finite & state, const size_t & chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_active_sites               (class_state_finite & state, const size_t & chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_next_sites                 (class_state_finite & state, const size_t & chi_lim, size_t num_sites = 4, std::optional<double> svd_threshold = std::nullopt);
    extern void project_to_closest_parity_sector    (class_state_finite & state, std::string paulistring);

    namespace internals{
        extern void set_product_state_in_parity_sector_from_bitset(class_state_finite & state, const std::string &parity_sector, const long state_number);
        extern void set_product_state_in_parity_sector_randomly(class_state_finite & state, const std::string &parity_sector);
        extern void set_product_state_randomly(class_state_finite & state, const std::string &parity_sector, bool use_pauli_eigenstates);
    }
}

/* clang-format on */
