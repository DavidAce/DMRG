#pragma once
#include <Eigen/Core>
#include <complex>
#include <config/enums.h>
#include <general/eigen_tensor_fwd_decl.h>
#include <list>
#include <optional>
#include <set>
#include <string>


namespace qm{
    class Gate;
}

/* clang-format off */
class class_state_finite;
namespace tools::finite::mps {
    using Scalar = std::complex<double>;
//    extern void move_center_point_single_site_fast   (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void move_center_point_single_site        (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void move_center_point                    (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void move_center_point_to_edge            (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void move_center_point_to_middle          (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void merge_multisite_tensor               (class_state_finite & state, const Eigen::Tensor<Scalar,3> & multisite_tensor, const std::vector<size_t> & sites, long center_position, long chi_lim, std::optional<double> svd_threshold = std::nullopt, std::optional<LogPolicy> log_policy = std::nullopt);
    extern bool normalize_state                      (class_state_finite & state,long chi_lim, std::optional<double> svd_threshold = std::nullopt, NormPolicy norm_policy = NormPolicy::IFNEEDED);
    extern void randomize_state                      (class_state_finite & state, StateInit state_type, StateInitType type,const std::string & sector, long chi_lim, bool use_eigenspinors, std::optional<long> bitfield = std::nullopt);
    extern void apply_random_paulis                  (class_state_finite & state, const std::vector<Eigen::Matrix2cd> & paulimatrices);
    extern void apply_random_paulis                  (class_state_finite & state, const std::vector<std::string> & paulistrings);
    extern void truncate_all_sites                   (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_active_sites                (class_state_finite & state, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void truncate_next_sites                  (class_state_finite & state, long chi_lim, size_t num_sites = 4, std::optional<double> svd_threshold = std::nullopt);
//    extern void apply_twosite_gates                  (class_state_finite & state, const std::vector<qm::Gate> & twosite_gates, bool inverse, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
//    extern void apply_twosite_gates                  (class_state_finite & state, const std::vector<Eigen::Tensor<Scalar,2>> & twosite_gates, bool inverse, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void apply_gates                          (class_state_finite & state, const std::vector<Eigen::Tensor<Scalar,2>> & nsite_tensors, size_t gate_size, bool reverse, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    extern void apply_gates                          (class_state_finite & state, const std::vector<qm::Gate> & gates, bool reverse, long chi_lim, std::optional<double> svd_threshold = std::nullopt);

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
        extern void set_product_state_aligned(class_state_finite & state, StateInitType type, const std::string &sector);
        extern void set_product_state_neel(class_state_finite & state, StateInitType type, const std::string &sector);

        // Entangled states
        extern void randomize_given_state (class_state_finite & state, StateInitType type);
        extern void set_random_entangled_state_in_sector_using_eigenspinors(class_state_finite & state, StateInitType type, const std::string &sector, long chi_lim);
        extern void set_random_entangled_state_with_random_spinors(class_state_finite & state, StateInitType type, long chi_lim);
    }
}

/* clang-format on */
