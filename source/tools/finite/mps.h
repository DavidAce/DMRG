#pragma once
#include <Eigen/Core>
#include <complex>
#include <config/enums.h>
#include <general/eigen_tensor_fwd_decl.h>
#include <list>
#include <optional>
#include <set>
#include <string>

/* clang-format off */
class class_state_finite;
namespace tools::finite::mps {
    using Scalar = std::complex<double>;
    extern void move_center_point               (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt); /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */
    extern void move_center_point_to_edge       (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt); /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */
    extern void merge_multisite_tensor          (class_state_finite & state, const Eigen::Tensor<Scalar,3> & multisite_tensor, const std::vector<size_t> & sites, size_t center_position, long chi_lim, std::optional<double> svd_threshold = std::nullopt, std::optional<LogPolicy> log_policy = std::nullopt);
    extern bool normalize_state                 (class_state_finite & state,long chi_lim, std::optional<double> svd_threshold = std::nullopt, NormPolicy norm_policy = NormPolicy::IFNEEDED);
    extern void randomize_state                 (class_state_finite & state, StateInit state_type, StateInitType type,const std::string & sector, long chi_lim, bool use_eigenspinors, std::optional<long> bitfield = std::nullopt);
    extern void apply_random_paulis             (class_state_finite & state, const std::vector<Eigen::Matrix2cd> & paulimatrices);
    extern void apply_random_paulis             (class_state_finite & state, const std::vector<std::string> & paulistrings);
    extern void truncate_all_sites              (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_active_sites           (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_next_sites             (class_state_finite & state, long chi_lim, size_t num_sites = 4, std::optional<double> svd_threshold = std::nullopt);
    extern void apply_twosite_operators         (class_state_finite & state, const std::vector<Eigen::Tensor<Scalar,2>> & twosite_operators, long chi_lim, std::optional<double> svd_threshold = std::nullopt);

    namespace internal{
        inline std::set<long> used_bitfields;
        extern bool bitfield_is_valid (std::optional<long> bitfield);
        extern int get_sign(const std::string &sector);
        extern std::string get_axis(const std::string &sector);
        extern Eigen::Vector2cd get_spinor(const std::string &axis, int sign);
        extern Eigen::Vector2cd get_spinor(const std::string &sector);
        extern Eigen::Matrix2cd get_pauli(const std::string &axis);
        extern std::vector<long> get_valid_bond_dimensions(size_t sizeplusone, long spin_dim,long chi_lim);

        extern void random_product_state (class_state_finite & state, StateInitType type, const std::string & sector, bool use_eigenspinors = false, std::optional<long> bitfield = std::nullopt);
        extern void random_entangled_state (class_state_finite & state, StateInitType type, const std::string & sector, long chi_lim, bool use_eigenspinors = false);

        // Product states
        extern void set_random_product_state_with_random_spinors(class_state_finite & state, StateInitType type);
        extern void set_random_product_state_on_axis_using_bitfield(class_state_finite & state, StateInitType type, const std::string &sector, long bitfield);
        extern void set_random_product_state_in_sector_using_eigenspinors(class_state_finite & state, StateInitType type, const std::string &sector);
        extern void set_random_product_state_on_axis(class_state_finite & state, StateInitType type, const std::string &sector);
        extern void set_product_state(class_state_finite & state, StateInitType type, const std::string &sector);

        // Entangled states
        extern void randomize_given_state (class_state_finite & state, StateInitType type, double factor = 1.0);
        extern void set_random_entangled_state_in_sector_using_eigenspinors(class_state_finite & state, StateInitType type, const std::string &sector, long chi_lim);
        extern void set_random_entangled_state_with_random_spinors(class_state_finite & state, StateInitType type, long chi_lim);
    }
}

/* clang-format on */
