#pragma once

#include "config/enums.h"
#include "math/svd/config.h"
#include "math/tenx/fwd_decl.h"
#include <complex>
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
enum class BitOrder {AsIs, Reverse};
namespace tools{
    extern std::string get_bitfield(size_t nbits, const std::string &pattern, BitOrder bitOrder = BitOrder::AsIs);
    extern std::string get_bitfield(size_t nbits, size_t pattern, BitOrder bitOrder = BitOrder::AsIs);
}


namespace tools::finite::mps {
    using cplx = std::complex<double>;
    extern size_t move_center_point_single_site      (StateFinite & state, std::optional<svd::config> svd_cfg = std::nullopt);
    extern size_t move_center_point                  (StateFinite & state, std::optional<svd::config> svd_cfg = std::nullopt);
    extern size_t move_center_point_to_pos           (StateFinite & state, long pos, std::optional<svd::config> svd_cfg = std::nullopt);
    extern size_t move_center_point_to_pos_dir       (StateFinite & state, long pos, int dir, std::optional<svd::config> svd_cfg = std::nullopt);
    extern size_t move_center_point_to_inward_edge   (StateFinite & state, std::optional<svd::config> svd_cfg = std::nullopt);
    extern size_t move_center_point_to_middle        (StateFinite & state, std::optional<svd::config> svd_cfg = std::nullopt);
    extern size_t merge_multisite_mps                (StateFinite & state, const Eigen::Tensor<cplx,3> & multisite_mps, const std::vector<size_t> & sites, long center_position, std::optional<svd::config> svd_cfg = std::nullopt, std::optional<LogPolicy> logPolicy = std::nullopt);
    extern bool normalize_state                      (StateFinite & state, std::optional<svd::config> svd_cfg = std::nullopt, NormPolicy norm_policy = NormPolicy::IFNEEDED);
    extern void initialize_state                     (StateFinite & state, StateInit state_type, StateInitType type, std::string_view sector, bool use_eigenspinors, long bond_lim, std::string & pattern);
    extern void apply_random_paulis                  (StateFinite & state, const std::vector<Eigen::Matrix2cd> & paulimatrices);
    extern void apply_random_paulis                  (StateFinite & state, const std::vector<std::string> & paulistrings);
    extern void truncate_all_sites                   (StateFinite & state, std::optional<svd::config> svd_cfg = std::nullopt);
    extern void truncate_active_sites                (StateFinite & state, std::optional<svd::config> svd_cfg = std::nullopt);
    extern void truncate_next_sites                  (StateFinite & state, size_t num_sites = 4, std::optional<svd::config> svd_cfg = std::nullopt);
    template<typename GateType>
    extern std::vector<size_t>
                generate_gate_sequence               (const StateFinite &state, const std::vector<GateType> &gates, CircuitOp cop, bool range_long_to_short = false);
    extern void apply_gate                           (StateFinite & state, const qm::Gate & gate, Eigen::Tensor<cplx, 3> & temp, GateOp gop, GateMove gmov, std::optional<svd::config> svd_cfg = std::nullopt);
    extern void apply_gates                          (StateFinite & state, const std::vector<Eigen::Tensor<cplx,2>> & nsite_tensors, size_t gate_size, CircuitOp cop, bool moveback = true, GateMove gm = GateMove::AUTO, std::optional<svd::config> svd_cfg = std::nullopt);
    extern void apply_gates                          (StateFinite & state, const std::vector<qm::Gate> & gates,  CircuitOp cop, bool moveback = true, GateMove gm = GateMove::AUTO, std::optional<svd::config> svd_cfg = std::nullopt);
    extern void apply_gates_old                      (StateFinite & state, const std::vector<qm::Gate> & gates,  CircuitOp cop, bool moveback = true,  std::optional<svd::config> svd_cfg = std::nullopt);
    extern void apply_circuit                        (StateFinite & state, const std::vector<std::vector<qm::Gate>> & gates, CircuitOp gop, bool mark_as_used = false, bool moveback = true, GateMove gm = GateMove::AUTO, std::optional<svd::config> svd_cfg = std::nullopt);
    extern void swap_sites                           (StateFinite & state, size_t posL, size_t posR, std::vector<size_t> & sites, GateMove gm);
    extern void apply_swap_gate                      (StateFinite & state, qm::SwapGate & gate, Eigen::Tensor<cplx, 3> & temp, GateOp gop, std::vector<size_t> & sites, GateMove gm, std::optional<svd::config> svd_cfg = std::nullopt);
    extern void apply_swap_gates                     (StateFinite & state, std::vector<qm::SwapGate> & gates, CircuitOp cop, GateMove gm = GateMove::AUTO, std::optional<svd::config> svd_cfg = std::nullopt);
    namespace init{
        inline std::set<size_t> used_bitfields;
        extern bool bitfield_is_valid (size_t bitfield);
        extern std::vector<long> get_valid_bond_dimensions(size_t sizeplusone, long spin_dim, long bond_lim);

        extern void random_product_state        (StateFinite & state, StateInitType type, std::string_view axis, bool use_eigenspinors, std::string &pattern);
        extern void random_entangled_state      (StateFinite & state, StateInitType type, std::string_view axis, bool use_eigenspinors, long bond_lim);

        // Product states
        extern void set_random_product_state_with_random_spinors(StateFinite & state, StateInitType type, std::string &pattern);
        extern void set_random_product_state_on_axis_using_eigenspinors(StateFinite & state, StateInitType type, std::string_view axis, std::string &pattern);
        extern void set_random_product_state_on_axis(StateFinite & state, StateInitType type, std::string_view axis, std::string &pattern);
        extern void set_product_state_domain_wall(StateFinite & state, StateInitType type, std::string_view axis, std::string  & pattern);
        extern void set_product_state_aligned(StateFinite & state, StateInitType type, std::string_view axis, std::string &pattern);
        extern void set_product_state_neel(StateFinite & state, StateInitType type, std::string_view axis, std::string  &pattern);
        extern void set_product_state_neel_shuffled (StateFinite & state, StateInitType type, std::string_view axis, std::string & pattern);
        extern void set_product_state_neel_dislocated (StateFinite & state, StateInitType type, std::string_view axis, std::string & pattern);
        extern void set_product_state_on_axis_using_pattern(StateFinite & state, StateInitType type, std::string_view axis, std::string & pattern);

        // Entangled states
        extern void randomize_given_state (StateFinite & state, StateInitType type);
        extern void set_midchain_singlet_neel_state(StateFinite & state, StateInitType type, std::string_view axis, std::string &pattern);
        extern void set_random_entangled_state_on_axes_using_eigenspinors(StateFinite & state, StateInitType type, const std::vector<std::string> & axes, long bond_lim);
        extern void set_random_entangled_state_on_axis_using_eigenspinors(StateFinite & state, StateInitType type, std::string_view axis, long bond_lim);
        extern void set_random_entangled_state_with_random_spinors(StateFinite & state, StateInitType type, long bond_lim);
        extern void set_random_entangled_state_haar(StateFinite & state, StateInitType type, long bond_lim);
    }
}

/* clang-format on */
