#pragma once

#include <complex>
#include <config/enums.h>
#include <general/eigen_tensor_fwd_decl.h>
#include <math/svd/settings.h>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace qm {
    class Gate;
}

/* clang-format off */
class StateFinite;
namespace tools::finite::mps {
    using Scalar = std::complex<double>;
//    extern void move_center_point_single_site_fast   (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_single_site      (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point                  (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_to_pos           (StateFinite & state, long pos, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_to_pos_dir       (StateFinite & state, long pos, int dir, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_to_edge          (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_to_middle        (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t merge_multisite_tensor             (StateFinite & state, const Eigen::Tensor<Scalar,3> & multisite_tensor, const std::vector<size_t> & positions, long center_position, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt, std::optional<LogPolicy> logPolicy = std::nullopt);
    extern bool normalize_state                      (StateFinite & state, std::optional<long> chi_lim = std::nullopt, std::optional<svd::settings> svd_settings = std::nullopt, NormPolicy norm_policy = NormPolicy::IFNEEDED);
    extern void randomize_state                      (StateFinite & state, StateInit state_type, StateInitType type,const std::string & sector, long chi_lim, bool use_eigenspinors, std::optional<long> bitfield = std::nullopt);
    extern void apply_random_paulis                  (StateFinite & state, const std::vector<Eigen::Matrix2cd> & paulimatrices);
    extern void apply_random_paulis                  (StateFinite & state, const std::vector<std::string> & paulistrings);
    extern void truncate_all_sites                   (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void truncate_active_sites                (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void truncate_next_sites                  (StateFinite & state, long chi_lim, size_t num_sites = 4, std::optional<svd::settings> svd_settings = std::nullopt);
//    extern void apply_twosite_gates                  (StateFinite & state, const std::vector<qm::Gate> & twosite_gates, bool inverse, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
//    extern void apply_twosite_gates                  (StateFinite & state, const std::vector<Eigen::Tensor<Scalar,2>> & twosite_gates, bool inverse, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_gates                          (StateFinite & state, const std::vector<Eigen::Tensor<Scalar,2>> & nsite_tensors, size_t gate_size, bool reverse, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_gates                          (StateFinite & state, const std::vector<qm::Gate> & gates, bool reverse, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);

    namespace init{
        inline std::set<long> used_bitfields;
        extern bool bitfield_is_valid (std::optional<long> bitfield);
        extern int get_sign(const std::string &sector);
        extern std::string get_axis(const std::string &sector);
        extern Eigen::Vector2cd get_spinor(const std::string &axis, int sign);
        extern Eigen::Vector2cd get_spinor(const std::string &sector);
        extern Eigen::Matrix2cd get_pauli(const std::string &axis);
        extern std::vector<long> get_valid_bond_dimensions(size_t sizeplusone, long spin_dim,long chi_lim);

        extern void random_product_state (StateFinite & state, StateInitType type, const std::string & sector, bool use_eigenspinors = false, std::optional<long> bitfield = std::nullopt);
        extern void random_entangled_state (StateFinite & state, StateInitType type, const std::string & sector, long chi_lim, bool use_eigenspinors = false);

        // Product states
        extern void set_random_product_state_with_random_spinors(StateFinite & state, StateInitType type);
        extern void set_random_product_state_on_axis_using_bitfield(StateFinite & state, StateInitType type, const std::string &sector, long bitfield);
        extern void set_random_product_state_in_sector_using_eigenspinors(StateFinite & state, StateInitType type, const std::string &sector);
        extern void set_random_product_state_on_axis(StateFinite & state, StateInitType type, const std::string &sector);
        extern void set_product_state_aligned(StateFinite & state, StateInitType type, const std::string &sector);
        extern void set_product_state_neel(StateFinite & state, StateInitType type, const std::string &sector);

        // Entangled states
        extern void randomize_given_state (StateFinite & state, StateInitType type);
        extern void set_random_entangled_state_in_sector_using_eigenspinors(StateFinite & state, StateInitType type, const std::string &sector, long chi_lim);
        extern void set_random_entangled_state_with_random_spinors(StateFinite & state, StateInitType type, long chi_lim);
    }
}

/* clang-format on */
