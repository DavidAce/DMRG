#include "../mps.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "math/num.h"
#include "math/rnd.h"
#include "math/tenx.h"
#include "qm/spin.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include <bitset>

void tools::finite::mps::init::random_product_state(StateFinite &state, StateInitType type, std::string_view sector, bool use_eigenspinors, long bitfield)
/*!
 * There are many ways to generate an initial product state based on the
 * arguments (sector,use_eigenspinors, bitfield) = (string,long, optional<bool> ).
 * Let
 *      * sector="axis"   means one of {"x","+x","-x","y","+y","-y","z","+z","-z"},
 *                        where the sign on an axis implies the global spin parity
 *                        sector. When the sign is present, the last spin is chosen
 *                        to end up in the correct sector.
 *                        In addition, sector can be defined as {"none", "random"}.
 *      * use_eigenspinors = true selects either |up> or |down> randomly on each site (rather than a mixture).
 *      * bitfield is only enabled if it is a non-negative number. The bits in this number are interpreted as [up/down].
 * Then
 *
 *       a) ("random" , ignored , ignored)    Set each spinor randomly on C2 (or R2 if type == REAL)
 *       b) ("axis", ignored, enabled)        Interpret seed_state as bitfield "01100010110..." and interpret these as
 *                                            up(0)/down(1) eigenspinors of the pauli matrix corresponding to "axis".
 *                                            Note that the axis sign is ignored.
 *                                            Note that the bitfield is used only once. Subsequent calls with the same bitfield
 *                                            will verert to case
 *       c) ("axis", true, disabled)          Set each |spinor> = |up> or |down> with 50% probability,
 *                                            and |up/down> are eigenspinors of the pauli matrix specified in "axis". If the global
 *                                            sign (+-) is omitted, a random sign is chosen with equal probabilities. In the x and
 *                                            z cases the global state becomes entirely real, which improves performance.
 *       d) ("axis", false, disabled)         Set each |spinor> = a|up> + b|down>, where a and b are random,
 *                                            real and normalized weights, and |up/down> are eigenspinors of the pauli matrix
 *                                            specified in "axis". Note that the axis sign is ignored.


 * Note: we "use" the bitfield only once. Subsequent calls do not keep resetting the seed.
*/
{
    tools::log->info("Setting random product state of type {} in sector {}", enum2sv(type), sector);
    state.clear_measurements();
    state.clear_cache();
    auto axis_valid = qm::spin::half::is_valid_axis(sector);
    if(sector == "random") {
        init::set_random_product_state_with_random_spinors(state, type); // a)
    } else if(init::bitfield_is_valid(bitfield) and axis_valid) {
        init::set_random_product_state_on_axis_using_bitfield(state, type, sector, bitfield); // b)
        init::used_bitfields.insert(bitfield);
    } else if(use_eigenspinors and axis_valid) {
        init::set_random_product_state_in_sector_using_eigenspinors(state, type, sector); // c)
    } else if(axis_valid) {
        init::set_random_product_state_on_axis(state, type, sector); // d)
    } else {
        throw except::runtime_error("Expected initial sector string: \"random\"|{}. Got \"{}\"", qm::spin::half::valid_axis_str, sector);
    }
}

void tools::finite::mps::init::set_product_state_aligned(StateFinite &state, StateInitType type, std::string_view sector) {
    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    auto axis = qm::spin::half::get_axis(sector);
    int  sign = qm::spin::half::get_sign(sector);
    if(type == StateInitType::REAL and axis == "y") throw std::runtime_error("StateInitType REAL incompatible with state in sector [y] which impliex CPLX");
    Eigen::Tensor<cplx, 3> spinor = tenx::TensorCast(qm::spin::half::get_spinor(axis, sign).normalized(), 2, 1, 1);
    tools::log->debug("Setting product state aligned using the |{}> eigenspinor of the pauli matrix σ{} on all sites", sign, axis);
    std::string label = "A";
    for(const auto &mps_ptr : state.mps_sites) {
        auto &mps = *mps_ptr;
        mps.set_mps(spinor, L, 0, label);
        if(mps.isCenter()) {
            mps.set_LC(L);
            label = "B";
        }
    }
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
}

void tools::finite::mps::init::set_product_state_neel(StateFinite &state, StateInitType type, std::string_view sector) {
    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    auto axis = qm::spin::half::get_axis(sector);
    if(type == StateInitType::REAL and axis == "y") throw std::runtime_error("StateInitType REAL incompatible with state in sector [y] which impliex CPLX");
    std::array<Eigen::Tensor<cplx, 3>, 2> spinors = {tenx::TensorCast(qm::spin::half::get_spinor(axis, +1).normalized(), 2, 1, 1),
                                                     tenx::TensorCast(qm::spin::half::get_spinor(axis, -1).normalized(), 2, 1, 1)};
    tools::log->debug("Setting product state neel using the |+-{}> eigenspinors of the pauli matrix σ{} on all sites", axis, axis);
    std::string label = "A";
    for(const auto &mps_ptr : state.mps_sites) {
        auto &&mps = *mps_ptr;
        auto   idx = num::mod<size_t>(mps.get_position(), 2);
        mps.set_mps(spinors.at(idx), L, 0, label);
        if(mps.isCenter()) {
            mps.set_LC(L);
            label = "B";
        }
    }
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
}

void tools::finite::mps::init::set_random_product_state_with_random_spinors(StateFinite &state, StateInitType type) {
    tools::log->info("Setting random product state with spinors in C²");
    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    std::string label = "A";
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps = *mps_ptr;
        if(type == StateInitType::CPLX)
            mps.set_mps(tenx::TensorCast(Eigen::VectorXcd::Random(2).normalized(), 2, 1, 1), L, 0, label);
        else if(type == StateInitType::REAL)
            mps.set_mps(tenx::TensorCast(Eigen::VectorXd::Random(2).normalized().cast<cplx>(), 2, 1, 1), L, 0, label);
        if(mps.isCenter()) {
            mps.set_LC(L);
            label = "B";
        }
    }
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
}

void tools::finite::mps::init::set_random_product_state_on_axis_using_bitfield(StateFinite &state, StateInitType type, std::string_view sector, long bitfield) {
    auto axis = qm::spin::half::get_axis(sector);
    tools::log->info("Setting random product state using the bitset of number {} to select eigenspinors of σ{}", bitfield, axis);

    if(bitfield < 0) throw except::runtime_error("Can't set product state from bitfield of negative number: {}", bitfield);
    if(type == StateInitType::REAL and axis == "y") throw except::runtime_error("StateInitType REAL incompatible with state in sector [y] which impliex CPLX");

    constexpr long maxbits = 64;
    if(maxbits < state.get_length()) throw except::range_error("Max supported state length for bitset is 64");
    std::bitset<maxbits>     bs(static_cast<size_t>(bitfield));
    std::vector<int>         bs_vec;
    std::vector<std::string> ud_vec;
    for(size_t i = 0; i < state.get_length(); i++) bs_vec.emplace_back(bs[i]);

    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    int         carry_sign = 1;
    std::string label      = "A";
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps  = *mps_ptr;
        int   sign = 2 * bs[mps.get_position()] - 1;
        carry_sign *= sign;
        mps.set_mps(tenx::TensorCast(qm::spin::half::get_spinor(sector).normalized(), 2, 1, 1), L, 0, label);
        std::string arrow = sign < 0 ? "↓" : "↑";
        ud_vec.emplace_back(arrow);
        if(mps.isCenter()) {
            mps.set_LC(L);
            label = "B";
        }
    }
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
}

void tools::finite::mps::init::set_random_product_state_in_sector_using_eigenspinors(StateFinite &state, StateInitType type, std::string_view sector) {
    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    auto axis      = qm::spin::half::get_axis(sector);
    int  sign      = qm::spin::half::get_sign(sector);
    int  last_sign = 1;
    if(type == StateInitType::REAL and axis == "y") throw std::runtime_error("StateInitType REAL incompatible with state in sector [y] which impliex CPLX");
    tools::log->info("Setting random product state in sector {} using eigenspinors of the pauli matrix σ{}", sector, axis);
    std::string label = "A";
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps = *mps_ptr;
        last_sign = 2 * rnd::uniform_integer_01() - 1;
        mps.set_mps(tenx::TensorCast(qm::spin::half::get_spinor(axis, last_sign).normalized(), 2, 1, 1), L, 0, label);
        if(mps.isCenter()) {
            mps.set_LC(L);
            label = "B";
        }
    }
    auto spin_component = tools::finite::measure::spin_component(state, axis);
    if(spin_component * sign < 0) {
        auto &mps = *state.mps_sites.back();
        mps.set_mps(tenx::TensorCast(qm::spin::half::get_spinor(axis, -last_sign).normalized(), 2, 1, 1), L, 0, label);
        spin_component = tools::finite::measure::spin_component(state, axis);
    }
    state.clear_measurements();
    if(spin_component * sign < 0) throw except::logic_error("Could not initialize_state in the correct sector");
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
}

void tools::finite::mps::init::set_random_product_state_on_axis(StateFinite &state, StateInitType type, std::string_view sector) {
    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    auto axis = qm::spin::half::get_axis(sector);
    tools::log->info("Setting random product state on axis {} using linear combinations of eigenspinors a|+> + b|-> of the pauli matrix σ{}", axis, axis);
    if(type == StateInitType::REAL and axis == "y") throw std::runtime_error("StateInitType REAL incompatible with state in sector [y] which impliex CPLX");
    auto        spinor_up = qm::spin::half::get_spinor(axis, 1);
    auto        spinor_dn = qm::spin::half::get_spinor(axis, -1);
    std::string label     = "A";
    for(auto &mps_ptr : state.mps_sites) {
        auto            &mps = *mps_ptr;
        Eigen::Vector2cd ab;
        if(type == StateInitType::REAL)
            ab = Eigen::Vector2d::Random().normalized();
        else
            ab = Eigen::Vector2cd::Random().normalized();

        Eigen::VectorXcd spinor = ab(0) * spinor_up + ab(1) * spinor_dn;
        mps.set_mps(tenx::TensorCast(spinor.normalized(), 2, 1, 1), L, 0, label);
        if(mps.isCenter()) {
            mps.set_LC(L);
            label = "B";
        }
    }
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
}
