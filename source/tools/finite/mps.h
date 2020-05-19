#pragma once
#include <complex>
#include <list>
#include <optional>
#include <string>
#include <tools/finite/opt-internals/enum_classes.h>
#include <unsupported/Eigen/CXX11/Tensor>
enum class ModelType;

/* clang-format off */
class class_state_finite;
namespace tools::finite::mps {
    using Scalar = std::complex<double>;
    extern void move_center_point               (class_state_finite & state, std::optional<long> chi_lim = std::nullopt, std::optional<double> svd_threshold = std::nullopt); /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */
    extern void merge_multisite_tensor          (class_state_finite & state, const Eigen::Tensor<Scalar,3> & multisite_tensor, const std::list<size_t> & sites, size_t center_position, std::optional<long> chi_lim = std::nullopt, std::optional<double> svd_threshold = std::nullopt);
    extern void normalize_state                 (class_state_finite & state, std::optional<long> chi_lim = std::nullopt,std::optional<double> svd_threshold = std::nullopt);

    extern void random_product_state            (class_state_finite & state, const std::string & sector, long state_number, bool use_pauli_eigenspinors = false);
    extern void random_current_state            (class_state_finite & state, const std::string & axis1, const std::string & axis2);
    extern void truncate_all_sites              (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_active_sites           (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_next_sites             (class_state_finite & state, long chi_lim, size_t num_sites = 4, std::optional<double> svd_threshold = std::nullopt);
    extern void project_to_closest_sector       (class_state_finite & state, const std::string & sector);
    namespace internals{
        extern void set_random_product_state_in_sector_using_bitset(class_state_finite & state, const std::string &sector, long state_number);
        extern void set_random_product_state_in_sector_using_eigenspinors(class_state_finite & state, const std::string &sector);
        extern void set_random_product_state(class_state_finite & state, const std::string &sector, bool use_pauli_eigenspinors);
    }
}

/* clang-format on */
