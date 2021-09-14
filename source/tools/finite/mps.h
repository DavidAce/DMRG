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
    class SwapGate;
}

/* clang-format off */
class StateFinite;
namespace tools::finite::mps {
    using cplx = std::complex<double>;
//    extern void move_center_point_single_site_fast   (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_single_site      (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point                  (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_to_pos           (StateFinite & state, long pos, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_to_pos_dir       (StateFinite & state, long pos, int dir, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_to_edge          (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_to_middle        (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t merge_multisite_tensor             (StateFinite & state, const Eigen::Tensor<cplx,3> & multisite_tensor, const std::vector<size_t> & positions, long center_position, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt, std::optional<LogPolicy> logPolicy = std::nullopt);
    extern bool normalize_state                      (StateFinite & state, std::optional<long> chi_lim = std::nullopt, std::optional<svd::settings> svd_settings = std::nullopt, NormPolicy norm_policy = NormPolicy::IFNEEDED);
    extern void randomize_state                      (StateFinite & state, StateInit state_type, StateInitType type, std::string_view sector, long chi_lim, bool use_eigenspinors, std::optional<long> bitfield = std::nullopt);
    extern void apply_random_paulis                  (StateFinite & state, const std::vector<Eigen::Matrix2cd> & paulimatrices);
    extern void apply_random_paulis                  (StateFinite & state, const std::vector<std::string> & paulistrings);
    extern void truncate_all_sites                   (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void truncate_active_sites                (StateFinite & state, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void truncate_next_sites                  (StateFinite & state, long chi_lim, size_t num_sites = 4, std::optional<svd::settings> svd_settings = std::nullopt);
//    extern void apply_twosite_gates                  (StateFinite & state, const std::vector<qm::Gate> & twosite_gates, bool inverse, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
//    extern void apply_twosite_gates                  (StateFinite & state, const std::vector<Eigen::Tensor<cplx,2>> & twosite_gates, bool inverse, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_gate                           (StateFinite & state, const qm::Gate & gate, Eigen::Tensor<cplx, 3> & temp, bool reverse, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_swap_gate                      (StateFinite & state, qm::SwapGate & gate, Eigen::Tensor<cplx, 3> & temp, bool reverse, long chi_lim, std::vector<size_t> & order, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_gates                          (StateFinite & state, const std::vector<Eigen::Tensor<cplx,2>> & nsite_tensors, size_t gate_size, bool reverse, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_gates                          (StateFinite & state, const std::vector<qm::Gate> & gates, bool reverse, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void swap_sites                           (StateFinite & state, size_t posL, size_t posR, std::vector<size_t> & order, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_swap_gates                     (StateFinite & state, std::vector<qm::SwapGate> & gates, bool reverse, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    namespace init{
        inline std::set<long> used_bitfields;
        extern bool bitfield_is_valid (std::optional<long> bitfield);
        extern int get_sign(std::string_view sector);
        extern std::string_view get_axis(std::string_view sector);
        extern Eigen::Vector2cd get_spinor(std::string_view axis, int sign);
        extern Eigen::Vector2cd get_spinor(std::string_view sector);
        extern Eigen::Matrix2cd get_pauli(std::string_view axis);
        extern std::vector<long> get_valid_bond_dimensions(size_t sizeplusone, long spin_dim,long chi_lim);

        extern void random_product_state (StateFinite & state, StateInitType type, std::string_view sector, bool use_eigenspinors = false, std::optional<long> bitfield = std::nullopt);
        extern void random_entangled_state (StateFinite & state, StateInitType type, std::string_view sector, long chi_lim, bool use_eigenspinors = false);

        // Product states
        extern void set_random_product_state_with_random_spinors(StateFinite & state, StateInitType type);
        extern void set_random_product_state_on_axis_using_bitfield(StateFinite & state, StateInitType type, std::string_view sector, long bitfield);
        extern void set_random_product_state_in_sector_using_eigenspinors(StateFinite & state, StateInitType type, std::string_view sector);
        extern void set_random_product_state_on_axis(StateFinite & state, StateInitType type, std::string_view sector);
        extern void set_product_state_aligned(StateFinite & state, StateInitType type, std::string_view sector);
        extern void set_product_state_neel(StateFinite & state, StateInitType type, std::string_view sector);

        // Entangled states
        extern void randomize_given_state (StateFinite & state, StateInitType type);
        extern void set_random_entangled_state_in_sector_using_eigenspinors(StateFinite & state, StateInitType type, std::string_view sector, long chi_lim);
        extern void set_random_entangled_state_with_random_spinors(StateFinite & state, StateInitType type, long chi_lim);
    }
}

/* clang-format on */
