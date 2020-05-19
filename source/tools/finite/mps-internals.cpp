//
// Created by david on 2019-08-12.
//
#include "mps.h"
#include <bitset>
#include <config/nmspc_settings.h>
#include <general/nmspc_quantum_mechanics.h>
#include <math/nmspc_random.h>
#include <tensors/state/class_environment.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>

using Scalar = tools::finite::mps::Scalar;

int get_sign(const std::string &sector) {
    if(sector.at(0) == '+')
        return 1;
    else if(sector.at(0) == '-')
        return -1;
    else
        return 0;
}

std::string get_axis(const std::string &sector) {
    int sign = get_sign(sector);
    if(sign == 0) {
        return sector.substr(0, 1);
    } else {
        return sector.substr(1, 1);
    }
}

Eigen::Vector2cd get_spinor(const std::string &axis, int sign) {
    if(axis == "x" and sign > 0) return qm::spinOneHalf::sx_spinors[0];
    if(axis == "x" and sign <= 0) return qm::spinOneHalf::sx_spinors[1];
    if(axis == "y" and sign > 0) return qm::spinOneHalf::sy_spinors[0];
    if(axis == "y" and sign <= 0) return qm::spinOneHalf::sy_spinors[1];
    if(axis == "z" and sign > 0) return qm::spinOneHalf::sz_spinors[0];
    if(axis == "z" and sign <= 0) return qm::spinOneHalf::sz_spinors[1];
    throw std::runtime_error(fmt::format("get_spinor given invalid axis: {}", axis));
}

Eigen::Vector2cd get_spinor(const std::string &sector) { return get_spinor(get_axis(sector), get_sign(sector)); }

void tools::finite::mps::internals::set_random_product_state_in_sector_using_bitset(class_state_finite &state, const std::string &sector, long state_number) {
    tools::log->trace("Setting product state from bitset of number {} in sector {}", state_number, sector);
    if(state_number < 0) {
        throw std::runtime_error(fmt::format("Can't set sector from bitset with negative state number: {}", state_number));
    }
    std::vector<std::string> ok_sectors        = {"x", "+x", "-x", "y", "+y", "-y", "z", "+z", "-z"};
    bool                     sector_is_defined = std::find(ok_sectors.begin(), ok_sectors.end(), sector) != ok_sectors.end();
    if(not sector_is_defined)
        throw std::logic_error(
            fmt::format("Can't use bitfield of state_number {} to set product state in sector [{}]. Choose one of (+-) x,y or z.", state_number, sector));

    constexpr long maxbits = 128;
    if(maxbits < state.get_length()) throw std::range_error("Max supported state length for bitset is 128");
    std::bitset<maxbits>     bs(static_cast<size_t>(state_number));
    std::vector<int>         bs_vec;
    std::vector<std::string> ud_vec;
    for(size_t i = 0; i < state.get_length(); i++) bs_vec.emplace_back(bs[i]);

    int                      sector_sign = get_sign(sector);
    Eigen::Tensor<Scalar, 1> L(1);
    L.setConstant(1.0);
    tools::log->info("Initializing state from bitfield of state number {} with eigvecs of σ{} in sector {}: {}", state_number, sector, sector, bs_vec);
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

    if(sector_sign * carry_sign == -1) {
        // Flip the last spin to get the correct total sign.
        auto &mps  = *state.mps_sites.back();
        int   sign = 2 * bs[mps.get_position()] - 1;
        sign *= -1;
        mps.set_mps(Textra::MatrixTensorMap(get_spinor(sector).normalized(), 2, 1, 1), L);
        std::string arrow = sign < 0 ? "↓" : "↑";
        ud_vec.back()     = arrow;
    }
    tools::log->info("Initialized  state from bitfield of state number {} with eigvecs of σ{} in sector {}: {}", state_number, sector, carry_sign, ud_vec);
}

void tools::finite::mps::internals::set_random_product_state_in_sector_using_eigenspinors(class_state_finite &state, const std::string &sector) {
    state.clear_measurements();
    state.clear_cache();
    Eigen::Tensor<Scalar, 1> L(1);
    std::string              axis      = get_axis(sector);
    int                      sign      = get_sign(sector);
    int                      last_sign = 1;
    tools::log->info("Setting random product state in sector {} using eigenspinors of the pauli matrix s{}...",sector,axis);
    L.setConstant(1.0);
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps = *mps_ptr;
        last_sign = 2 * rn::uniform_integer_01() - 1;
        mps.set_mps(Textra::MatrixTensorMap(get_spinor(axis).normalized(), 2, 1, 1), L);
        if(mps.isCenter()) mps.set_LC(L);
    }
    auto spin_component = tools::finite::measure::spin_component(state, axis);
    if(spin_component * sign < 0) {
        state.mps_sites.back()->set_mps(Textra::MatrixTensorMap(get_spinor(axis, -last_sign).normalized(), 2, 1, 1), L);
        spin_component = tools::finite::measure::spin_component(state, axis);
    }
    state.clear_measurements();
    tools::log->info("Setting random product state in sector {} using eigenspinors of the pauli matrix s{}... global spin component {}: OK",sector,axis,spin_component);
    if(spin_component * sign < 0) throw std::logic_error("Could not initialize_state in the correct sector");
}

void tools::finite::mps::internals::set_random_product_state(class_state_finite &state, const std::string &sector, bool use_pauli_eigenspinors) {
    std::vector<std::string> ok_sectors        = {"x", "+x", "-x", "y", "+y", "-y", "z", "+z", "-z"};
    bool                     sector_is_defined = std::find(ok_sectors.begin(), ok_sectors.end(), sector) != ok_sectors.end();
    if(sector_is_defined and use_pauli_eigenspinors) {
        // Case a)
        set_random_product_state_in_sector_using_eigenspinors(state, sector);
    } else if(sector_is_defined and not use_pauli_eigenspinors) {
        set_random_product_state(state, "random", false);
        state = tools::finite::ops::get_projection_to_closest_parity_sector(state, sector);
    } else if(sector == "randomAxis") {
        std::vector<std::string> possibilities = {"x", "y", "z"};
        std::string              chosen_axis   = possibilities[static_cast<unsigned long>(rn::uniform_integer_box(0, 2))];
        set_random_product_state_in_sector_using_eigenspinors(state, chosen_axis);
    } else if(sector == "random") {
        Eigen::Tensor<Scalar, 1> L(1);
        L.setConstant(1.0);
        for(auto &mps_ptr : state.mps_sites) {
            auto &mps = *mps_ptr;
            mps.set_mps(Textra::MatrixTensorMap(Eigen::VectorXcd::Random(2).normalized(), 2, 1, 1), L);
            if(mps.isCenter()) mps.set_LC(L);
        }
    } else if(sector == "none") {
        return;
    } else {
        throw std::runtime_error(fmt::format("Wrong sector string. Expected one of (+-) x, y, z, randomAxis, random or none. Got: {}", sector));
    }
}
