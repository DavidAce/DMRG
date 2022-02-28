#pragma once

#include <complex>
#include <config/enums.h>
#include <math/svd/settings.h>
#include <math/tenx/fwd_decl.h>
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
//    extern void move_center_point_single_site_fast   (StateFinite & state, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_single_site      (StateFinite & state, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point                  (StateFinite & state, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_to_pos           (StateFinite & state, long pos, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_to_pos_dir       (StateFinite & state, long pos, int dir, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_to_edge          (StateFinite & state, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t move_center_point_to_middle        (StateFinite & state, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt);
    extern size_t merge_multisite_mps                (StateFinite & state, const Eigen::Tensor<cplx,3> & multisite_mps, const std::vector<size_t> & sites, long center_position, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt, std::optional<LogPolicy> logPolicy = std::nullopt);
    extern bool normalize_state                      (StateFinite & state, std::optional<long> bond_limit = std::nullopt, std::optional<svd::settings> svd_settings = std::nullopt, NormPolicy norm_policy = NormPolicy::IFNEEDED);
    extern void randomize_state                      (StateFinite & state, StateInit state_type, StateInitType type, std::string_view sector, long bond_limit, bool use_eigenspinors, std::optional<long> bitfield = std::nullopt);
    extern void apply_random_paulis                  (StateFinite & state, const std::vector<Eigen::Matrix2cd> & paulimatrices);
    extern void apply_random_paulis                  (StateFinite & state, const std::vector<std::string> & paulistrings);
    extern void truncate_all_sites                   (StateFinite & state, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void truncate_active_sites                (StateFinite & state, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void truncate_next_sites                  (StateFinite & state, long bond_limit, size_t num_sites = 4, std::optional<svd::settings> svd_settings = std::nullopt);
//    extern void apply_twosite_gates                  (StateFinite & state, const std::vector<qm::Gate> & twosite_gates, bool inverse, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt);
//    extern void apply_twosite_gates                  (StateFinite & state, const std::vector<Eigen::Tensor<cplx,2>> & twosite_gates, bool inverse, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_gate                           (StateFinite & state, const qm::Gate & gate, Eigen::Tensor<cplx, 3> & temp, bool reverse, long bond_limit, GateMove gm, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_gates                          (StateFinite & state, const std::vector<Eigen::Tensor<cplx,2>> & nsite_tensors, size_t gate_size, bool reverse, long bond_limit, GateMove gm = GateMove::AUTO, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_gates                          (StateFinite & state, const std::vector<qm::Gate> & gates, bool reverse, long bond_limit, GateMove gm = GateMove::AUTO, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_gates_old                      (StateFinite &state, const std::vector<qm::Gate> &gates, bool reverse, long bond_limit, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void swap_sites                           (StateFinite & state, size_t posL, size_t posR, std::vector<size_t> & order, GateMove gm, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_swap_gate                      (StateFinite & state, qm::SwapGate & gate, Eigen::Tensor<cplx, 3> & temp, bool reverse, long bond_limit, std::vector<size_t> & order, GateMove gm, std::optional<svd::settings> svd_settings = std::nullopt);
    extern void apply_swap_gates                     (StateFinite & state, std::vector<qm::SwapGate> & gates, bool reverse, long bond_limit, GateMove gm = GateMove::AUTO, std::optional<svd::settings> svd_settings = std::nullopt);
    namespace init{
        inline std::set<long> used_bitfields;
        extern bool bitfield_is_valid (std::optional<long> bitfield);
        extern std::vector<long> get_valid_bond_dimensions(size_t sizeplusone, long spin_dim,long bond_limit);

        extern void random_product_state (StateFinite & state, StateInitType type, std::string_view sector, bool use_eigenspinors = false, std::optional<long> bitfield = std::nullopt);
        extern void random_entangled_state (StateFinite & state, StateInitType type, std::string_view sector, long bond_limit, bool use_eigenspinors = false);

        // Product states
        extern void set_random_product_state_with_random_spinors(StateFinite & state, StateInitType type);
        extern void set_random_product_state_on_axis_using_bitfield(StateFinite & state, StateInitType type, std::string_view sector, long bitfield);
        extern void set_random_product_state_in_sector_using_eigenspinors(StateFinite & state, StateInitType type, std::string_view sector);
        extern void set_random_product_state_on_axis(StateFinite & state, StateInitType type, std::string_view sector);
        extern void set_product_state_aligned(StateFinite & state, StateInitType type, std::string_view sector);
        extern void set_product_state_neel(StateFinite & state, StateInitType type, std::string_view sector);

        // Entangled states
        extern void randomize_given_state (StateFinite & state, StateInitType type);
        extern void set_random_entangled_state_in_sector_using_eigenspinors(StateFinite & state, StateInitType type, std::string_view sector, long bond_limit);
        extern void set_random_entangled_state_with_random_spinors(StateFinite & state, StateInitType type, long bond_limit);
    }
}

/* clang-format on */
