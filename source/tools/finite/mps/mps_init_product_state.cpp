#include "../mps.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "math/num.h"
#include "math/rnd.h"
#include "math/tenx.h"
#include "qm/spin.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include <bitset>

void tools::finite::mps::init::random_product_state(StateFinite &state, StateInitType type, std::string_view axis, bool use_eigenspinors, size_t bitfield)
/*!
 * There are many ways to generate an initial product state based on the
 * arguments (axis ,use_eigenspinors, bitfield) = (string,bool,size_t).
 * Let
 *      * axis="str"         where str is valid if it is {"x","+x","-x","y","+y","-y","z","+z","-z"},
 *                           and the sign corresponds to the sign of the global spin compoonent (parity sector).
 *                           When the sign is present, the spin on the last site is set manually to end up on the correct axis and parity sector.
 *      * use_eigenspinors = true selects on of {|up>, |down>} along axis randomly on each site,
 *                           whereas false selects a|up> + b|down> with random a,b on each site.
 *      * bitfield is only enabled if it is a non-negative number. The bits in this number are interpreted as [up/down].
 * Then
 *
 *       a) ("random" , ignored , ignored)    Set each spinor randomly on C2 (or R2 if type == REAL)
 *       b) (valid axis, ignored, enabled)    Interpret seed_state as bitfield "01100010110..." and interpret these as
 *                                            up(0)/down(1) eigenspinors of the pauli matrix corresponding to "axis".
 *                                            Note that the axis sign is ignored.
 *                                            Note that the bitfield is used only once. Subsequent calls with the same bitfield
 *                                            will verert to case
 *       c) (valid axis, true, disabled)      Set each |spinor> = |up> or |down> with 50% probability,
 *                                            and |up/down> are eigenspinors of the pauli matrix specified in "axis". If the global
 *                                            sign (+-) is omitted, a random sign is chosen with equal probabilities. In the x and
 *                                            z cases the global state becomes entirely real, which improves performance.
 *       d) (valid axis, false, disabled)     Set each |spinor> = a|up> + b|down>, where a and b are random,
 *                                            real and normalized weights, and |up/down> are eigenspinors of the pauli matrix
 *                                            specified in "axis". Note that the axis sign is ignored.


 * Note: we "use" the bitfield only once. Subsequent calls do not keep resetting the seed.
*/
{
    tools::log->info("Setting random product state of type {} on axis {}", enum2sv(type), axis);
    state.clear_measurements();
    state.clear_cache();
    auto axis_valid = qm::spin::half::is_valid_axis(axis);
    if(axis == "random") {
        init::set_random_product_state_with_random_spinors(state, type);                    // a)
    } else if(init::bitfield_is_valid(bitfield) and axis_valid) {
        init::set_random_product_state_on_axis_using_bitfield(state, type, axis, bitfield); // b)
        init::used_bitfields.insert(bitfield);
    } else if(use_eigenspinors and axis_valid) {
        init::set_random_product_state_on_axis_using_eigenspinors(state, type, axis); // c)
    } else if(axis_valid) {
        init::set_random_product_state_on_axis(state, type, axis);                    // d)
    } else {
        throw except::runtime_error("Expected initial axis string: \"random\"|{}. Got \"{}\"", qm::spin::half::valid_axis_str, axis);
    }
}

void tools::finite::mps::init::product_state_neel_shuffled(StateFinite &state, StateInitType type, std::string_view axis, std::vector<size_t> &pattern) {
    tools::log->info("Setting randomly shuffled Néel state of type {} on axis {} {}", enum2sv(type), axis,
                     pattern.empty() ? "" : fmt::format(" | from pattern: {}", pattern));
    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    auto axus = qm::spin::half::get_axis_unsigned(axis);
    if(type == StateInitType::REAL and axus == "y") throw std::runtime_error("StateInitType REAL incompatible with state on axis [y] which impliex CPLX");
    std::array<Eigen::Tensor<cplx, 3>, 2> spinors = {tenx::TensorCast(qm::spin::half::get_spinor(axus, +1).normalized(), 2, 1, 1),
                                                     tenx::TensorCast(qm::spin::half::get_spinor(axus, -1).normalized(), 2, 1, 1)};
    std::array<std::string_view, 2>       arrows  = {"↓", "↑"};
    std::string                           str;
    if(pattern.size() != state.get_length()) {
        pattern.resize(state.get_length(), 0);
        for(auto &&[i, s] : iter::enumerate(pattern)) s = num::mod<size_t>(i, 2); // Set Neel pattern 010101010101...
        std::shuffle(pattern.begin(), pattern.end(), rnd::internal::rng);         // Shuffle the sequence randomly like 101011110100...
    }

    std::string label = "A";
    for(const auto &mps_ptr : state.mps_sites) {
        auto &&mps = *mps_ptr;
        auto   idx = pattern.at(mps.get_position<size_t>());
        mps.set_mps(spinors.at(idx), L, 0, label);
        str.append(arrows.at(idx));
        if(mps.isCenter()) {
            mps.set_LC(L);
            label = "B";
        }
    }
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
    tools::log->info("Initial state: {}", str);
}

void tools::finite::mps::init::set_product_state_aligned(StateFinite &state, StateInitType type, std::string_view axis) {
    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    auto axus = qm::spin::half::get_axis_unsigned(axis);
    int  sign = qm::spin::half::get_sign(axis);
    if(type == StateInitType::REAL and axis == "y") throw std::runtime_error("StateInitType REAL incompatible with state in axis [y] which impliex CPLX");
    Eigen::Tensor<cplx, 3> spinor = tenx::TensorCast(qm::spin::half::get_spinor(axus, sign).normalized(), 2, 1, 1);
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

void tools::finite::mps::init::set_product_state_neel(StateFinite &state, StateInitType type, std::string_view axis) {
    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    auto axus = qm::spin::half::get_axis_unsigned(axis);
    if(type == StateInitType::REAL and axus == "y") throw std::runtime_error("StateInitType REAL incompatible with state on axis [y] which impliex CPLX");
    std::array<Eigen::Tensor<cplx, 3>, 2> spinors = {tenx::TensorCast(qm::spin::half::get_spinor(axus, +1).normalized(), 2, 1, 1),
                                                     tenx::TensorCast(qm::spin::half::get_spinor(axus, -1).normalized(), 2, 1, 1)};
    tools::log->debug("Setting product state neel using the |+-{}> eigenspinors of the pauli matrix σ{} on all sites", axus, axus);
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

void tools::finite::mps::init::set_random_product_state_with_gaussian_spinors(StateFinite &state, StateInitType type) {
    tools::log->info("Setting random product state with gaussian spinors");
    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    std::string label    = "A";
    auto        gaussian = [&type]() -> std::complex<double> {
        switch(type) {
            case StateInitType::REAL: return std::complex<double>(rnd::normal(0, 1), 0);
            case StateInitType::CPLX: {
                auto re = rnd::normal(0, 1);
                auto im = rnd::normal(0, 1);
                return std::complex<double>(re, im);
            }
            default: return std::complex<double>(rnd::normal(0, 1), 0);
        }
    };
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps = *mps_ptr;
        mps.set_mps(tenx::TensorCast(Eigen::VectorXcd::NullaryExpr(2, gaussian).normalized(), 2, 1, 1), L, 0, label);
        if(mps.isCenter()) {
            mps.set_LC(L);
            label = "B";
        }
    }
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
}

void tools::finite::mps::init::set_random_product_state_on_axis_using_bitfield(StateFinite &state, StateInitType type, std::string_view axis, size_t bitfield,
                                                                               LogPolicy logPolicy) {
    auto axus = qm::spin::half::get_axis_unsigned(axis);
    if(bitfield == -1ul) throw except::logic_error("Can't set product state from bitfield == -1ul");
    if(type == StateInitType::REAL and axus == "y") throw except::runtime_error("StateInitType REAL incompatible with state on axis [y] which impliex CPLX");

    constexpr size_t maxbits = 64;
    if(maxbits < state.get_length()) throw except::range_error("Max supported state length for bitset is 64");
    std::bitset<maxbits>     bs(bitfield);
    std::vector<int>         bs_vec;
    std::vector<std::string> ud_vec;
    for(size_t i = 0; i < state.get_length(); i++) bs_vec.emplace_back(bs[i]);

    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    std::string label = "A";
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps  = *mps_ptr;
        int   sign = 2 * bs[mps.get_position()] - 1;
        mps.set_mps(tenx::TensorCast(qm::spin::half::get_spinor(axus, sign).normalized(), 2, 1, 1), L, 0, label);
        ud_vec.emplace_back((sign < 0 ? "↓" : "↑"));
        if(mps.isCenter()) {
            mps.set_LC(L);
            label = "B";
        }
    }
    if(logPolicy == LogPolicy::NORMAL)
        tools::log->info("Set random product state using the bitset of number {} to select eigenspinors of σ{}: {}", bitfield, axus, fmt::join(ud_vec, ""));
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
}

void tools::finite::mps::init::set_random_product_state_on_axis_using_eigenspinors(StateFinite &state, StateInitType type, std::string_view axis,
                                                                                   LogPolicy logPolicy) {
    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    auto axus     = qm::spin::half::get_axis_unsigned(axis);
    int  sign_tgt = qm::spin::half::get_sign(axis);
    int  sign_now = 1;
    if(type == StateInitType::REAL and axus == "y") throw std::runtime_error("StateInitType REAL incompatible with state on axis [y] which impliex CPLX");
    std::string              label = "A";
    std::vector<std::string> ud_vec;
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps  = *mps_ptr;
        auto  sign = 2 * rnd::uniform_integer_01() - 1;
        mps.set_mps(tenx::TensorCast(qm::spin::half::get_spinor(axus, sign_now).normalized(), 2, 1, 1), L, 0, label);
        if(mps.isCenter()) {
            mps.set_LC(L);
            label = "B";
        }
        ud_vec.emplace_back((sign < 0 ? "↓" : "↑"));
        sign_now *= sign;
    }
    if(sign_now * sign_tgt < 0) {
        // Flip the last spin to get the correct overall sign
        auto &mps     = *state.mps_sites.back();
        auto  sign    = -sign_now;
        ud_vec.back() = (sign < 0 ? "↓" : "↑");
        mps.set_mps(tenx::TensorCast(qm::spin::half::get_spinor(axus, sign).normalized(), 2, 1, 1), L, 0, label);
    }
    if(logPolicy == LogPolicy::NORMAL)
        tools::log->info("Set random product state on axis {} using eigenspinors of the pauli matrix σ{}: {}", axis, axus, fmt::join(ud_vec, ""));
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
}

void tools::finite::mps::init::set_random_product_state_on_axis(StateFinite &state, StateInitType type, std::string_view axis) {
    Eigen::Tensor<cplx, 1> L(1);
    L.setConstant(1.0);
    auto axus = qm::spin::half::get_axis_unsigned(axis);
    tools::log->info("Setting random product state on axis {} using linear combinations of eigenspinors a|+> + b|-> of the pauli matrix σ{}", axus, axus);
    if(type == StateInitType::REAL and axus == "y") throw std::runtime_error("StateInitType REAL incompatible with state on axis [y] which impliex CPLX");
    auto        spinor_up = qm::spin::half::get_spinor(axus, 1);
    auto        spinor_dn = qm::spin::half::get_spinor(axus, -1);
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
