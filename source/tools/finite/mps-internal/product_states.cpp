//
// Created by david on 2019-08-12.
//

#include "tools/finite/mps.h"
#include <bitset>
#include <config/nmspc_settings.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <math/rnd.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>

using Scalar = tools::finite::mps::Scalar;

int tools::finite::mps::internal::get_sign(const std::string &sector) {
    if(sector.at(0) == '+') return 1;
    else if(sector.at(0) == '-')
        return -1;
    else
        return 0;
}

std::string tools::finite::mps::internal::get_axis(const std::string &sector) {
    std::vector<std::string> valid_axis_str = {"x", "+x", "-x", "y", "+y", "-y", "z", "+z", "-z"};
    bool                     axis_is_valid  = std::find(valid_axis_str.begin(), valid_axis_str.end(), sector) != valid_axis_str.end();
    if(not axis_is_valid) throw std::runtime_error(fmt::format("Could not extract valid axis from sector string [{}]. Choose one of (+-) x,y or z.", sector));
    int sign = get_sign(sector);
    if(sign == 0) {
        return sector.substr(0, 1);
    } else {
        return sector.substr(1, 1);
    }
}

Eigen::Vector2cd tools::finite::mps::internal::get_spinor(const std::string &axis, int sign) {
    if(axis == "x" and sign > 0) return qm::spinOneHalf::sx_spinors[0];
    if(axis == "x" and sign <= 0) return qm::spinOneHalf::sx_spinors[1];
    if(axis == "y" and sign > 0) return qm::spinOneHalf::sy_spinors[0];
    if(axis == "y" and sign <= 0) return qm::spinOneHalf::sy_spinors[1];
    if(axis == "z" and sign > 0) return qm::spinOneHalf::sz_spinors[0];
    if(axis == "z" and sign <= 0) return qm::spinOneHalf::sz_spinors[1];
    throw std::runtime_error(fmt::format("get_spinor given invalid axis: {}", axis));
}

Eigen::Vector2cd tools::finite::mps::internal::get_spinor(const std::string &sector) { return get_spinor(get_axis(sector), get_sign(sector)); }

Eigen::Matrix2cd tools::finite::mps::internal::get_pauli(const std::string &axis) {
    if(axis.find('x') != std::string::npos) return qm::spinOneHalf::sx;
    if(axis.find('y') != std::string::npos) return qm::spinOneHalf::sy;
    if(axis.find('z') != std::string::npos) return qm::spinOneHalf::sz;
    if(axis.find('I') != std::string::npos) return qm::spinOneHalf::Id;
    throw std::runtime_error(fmt::format("get_pauli given invalid axis: {}", axis));
}

void tools::finite::mps::internal::random_product_state(class_state_finite &state, const std::string &sector, bool use_eigenspinors,
                                                        std::optional<long> bitfield)
/*!
 * There are many ways to generate an initial product state based on the
 * arguments (sector,use_eigenspinors, bitfield) = (string,long, optional<bool> ).
 * Let
 *      * sector="+-axis" mean one of {"x","+x","-x","y","+y","-y","z","+z","-z"},
 *                        where the sign on an axis implies the global spin parity.
 *                        In addition, sector can be defined as {"none", "random"}.
 *      * bitfield is only enabled if it is a non-negative number. The bits in this number are interpreted as [up/down].
 * Then
 *
 *       a) ("none", ignored, ignored)        Does not randomize.
 *       b) ("random" , ignored , ignored)    Set each spinor randomly on C2
 *       c) ("axis", ignored, enabled)        Interpret seed_state as bitfield "01100010110..." and interpret these as
 *                                            up(0)/down(1) eigenspinors of the pauli matrix corresponding to "axis".
 *                                            Note that the axis sign is ignored.
 *                                            Note that the bitfield is used only once. Subsequent calls with the same bitfield
 *                                            will verert to case
 *       d) ("+-axis", true, disabled)        Set each |spinor> = a|up> + (1-a)|down>, where a is either 0 or 1 with 50% probability,
 *                                            and |up/down> are eigenspinors of the pauli matrix specified in "axis". If the global
 *                                            sign (+-) is omitted, a random sign is chosen with equal probabilities. In the x and
 *                                            z cases the global state becomes entirely real, which improves performance.
 *       e) ("axis", false, disabled)         Set each |spinor> = a|up> + b|down>, where a and b are random,
 *                                            real and normalized weights, and |up/down> are eigenspinors of the pauli matrix
 *                                            specified in "axis". Note that the axis sign is ignored.

 * Note: we "use" the bitfield only once. Subsequent calls do not keep resetting the seed.
*/
{
    if(sector == "none") return; // a)
    tools::log->debug("Setting random product state in sector {}", sector);
    state.clear_measurements();
    state.clear_cache();
    if(sector == "random") {
        internal::set_random_product_state_with_spinors_in_c2(state); // b)
    } else if(bitfield_is_valid(bitfield)) {
        internal::set_random_product_state_on_axis_using_bitfield(state, sector, bitfield.value()); // c)
        internal::used_bitfields.insert(bitfield.value());
    } else if(use_eigenspinors) {
        internal::set_random_product_state_in_sector_using_eigenspinors(state, sector); // d)
    } else {
        internal::set_random_product_state_on_axis(state, sector); // e)
    }
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_have_been_updated(true); // This operation changes all sites
}

void tools::finite::mps::internal::set_product_state(class_state_finite &state, const std::string &sector) {
    Eigen::Tensor<Scalar, 1> L(1);
    L.setConstant(1.0);
    std::string axis      = get_axis(sector);
    int         sign      = get_sign(sector);
    int         last_sign = 1;
    auto        spinor    = Textra::MatrixTensorMap(get_spinor(axis, last_sign).normalized(), 2, 1, 1);
    tools::log->info("Setting product state using the |{}> eigenspinor of the pauli matrix σ{} on all sites...", sign, axis);
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps = *mps_ptr;
        mps.set_mps(spinor, L);
        if(mps.isCenter()) mps.set_LC(L);
    }
    state.clear_measurements();
    tools::log->info("Setting product state using the |{}> eigenspinor of the pauli matrix σ{} on all sites... OK", sign, axis);
}

void tools::finite::mps::internal::set_random_product_state_with_spinors_in_c2(class_state_finite &state) {
    Eigen::Tensor<Scalar, 1> L(1);
    L.setConstant(1.0);
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps = *mps_ptr;
        mps.set_mps(Textra::MatrixTensorMap(Eigen::VectorXcd::Random(2).normalized(), 2, 1, 1), L);
        if(mps.isCenter()) mps.set_LC(L);
    }
}

void tools::finite::mps::internal::set_random_product_state_on_axis_using_bitfield(class_state_finite &state, const std::string &sector, long bitfield) {
    auto axis = get_axis(sector);
    tools::log->trace("Setting product state from bitset of number {} in along axis {}", bitfield, axis);
    if(bitfield < 0) throw std::runtime_error(fmt::format("Can't set product state from bitfield of negative number: {}", bitfield));

    constexpr long maxbits = 64;
    if(maxbits < state.get_length()) throw std::range_error("Max supported state length for bitset is 64");
    std::bitset<maxbits>     bs(static_cast<size_t>(bitfield));
    std::vector<int>         bs_vec;
    std::vector<std::string> ud_vec;
    for(size_t i = 0; i < state.get_length(); i++) bs_vec.emplace_back(bs[i]);

    Eigen::Tensor<Scalar, 1> L(1);
    L.setConstant(1.0);
    tools::log->info("Initializing product state from bitfield of number {} with eigenspinors of σ{}: {}", bitfield, axis, bs_vec);
    int carry_sign = 1;
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps  = *mps_ptr;
        int   sign = 2 * bs[mps.get_position()] - 1;
        carry_sign *= sign;
        mps.set_mps(Textra::MatrixTensorMap(get_spinor(sector).normalized(), 2, 1, 1), L);
        std::string arrow = sign < 0 ? "↓" : "↑";
        ud_vec.emplace_back(arrow);
        if(mps.isCenter()) mps.set_LC(L);
    }
    tools::log->info("Initialized state from bitfield of number {} with eigenspinors of σ{} in sector {}: {}", bitfield, sector, carry_sign, ud_vec);
}

void tools::finite::mps::internal::set_random_product_state_in_sector_using_eigenspinors(class_state_finite &state, const std::string &sector) {
    Eigen::Tensor<Scalar, 1> L(1);
    L.setConstant(1.0);
    std::string axis      = get_axis(sector);
    int         sign      = get_sign(sector);
    int         last_sign = 1;
    tools::log->info("Setting random product state in sector {} using eigenspinors of the pauli matrix σ{}...", sector, axis);
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps = *mps_ptr;
        last_sign = 2 * rnd::uniform_integer_01() - 1;
        mps.set_mps(Textra::MatrixTensorMap(get_spinor(axis, last_sign).normalized(), 2, 1, 1), L);
        if(mps.isCenter()) mps.set_LC(L);
    }
    auto spin_component = tools::finite::measure::spin_component(state, axis);
    if(spin_component * sign < 0) {
        auto &mps = *state.mps_sites.back();
        mps.set_mps(Textra::MatrixTensorMap(get_spinor(axis, -last_sign).normalized(), 2, 1, 1), L);
        spin_component = tools::finite::measure::spin_component(state, axis);
    }
    state.clear_measurements();
    tools::log->info("Setting random product state in sector {} using eigenspinors of the pauli matrix σ{}... global spin component {}: OK", sector, axis,
                     spin_component);
    if(spin_component * sign < 0) throw std::logic_error("Could not initialize_state in the correct sector");
}

void tools::finite::mps::internal::set_random_product_state_on_axis(class_state_finite &state, const std::string &sector) {
    Eigen::Tensor<Scalar, 1> L(1);
    L.setConstant(1.0);
    std::string axis = get_axis(sector);
    tools::log->info("Setting random product state on axis {} using linear combinations of eigenspinors a|+> + b|-> of the pauli matrix σ{}...", sector, axis);
    auto spinor_up = internal::get_spinor(axis, 1);
    auto spinor_dn = internal::get_spinor(axis, -1);
    for(auto &mps_ptr : state.mps_sites) {
        auto &           mps    = *mps_ptr;
        Eigen::VectorXd  ab     = Eigen::VectorXd::Random(2).normalized();
        Eigen::VectorXcd spinor = ab(0) * spinor_up + ab(1) * spinor_dn;
        mps.set_mps(Textra::MatrixTensorMap(spinor.normalized(), 2, 1, 1), L);
        if(mps.isCenter()) mps.set_LC(L);
    }
}
