//
// Created by david on 2020-06-07.
//

#include "tools/finite/mps.h"
#include <bitset>
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <math/rnd.h>
#include <physics/nmspc_quantum_mechanics.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>

using Scalar = tools::finite::mps::Scalar;

std::vector<long> tools::finite::mps::internal::get_valid_bond_dimensions(size_t sizeplusone, long spin_dim, long chi_lim) {
    // Construct a valid set of bond dimensions
    std::vector<long> bond_dimensions(sizeplusone, chi_lim);
    bond_dimensions.front() = 1;
    bond_dimensions.back()  = 1;
    for(size_t i = 1; i < bond_dimensions.size() / 2; i++) {
        auto max_bond_dim  = spin_dim * std::min(bond_dimensions[i - 1], bond_dimensions[i + 1]);
        bond_dimensions[i] = std::min(bond_dimensions[i], max_bond_dim);
    }
    for(size_t i = bond_dimensions.size() - 1; i >= bond_dimensions.size() / 2; i--) {
        auto max_bond_dim  = spin_dim * std::min(bond_dimensions[i - 1], bond_dimensions[i + 1]);
        bond_dimensions[i] = std::min(bond_dimensions[i], max_bond_dim);
    }
    return bond_dimensions;
}

void tools::finite::mps::internal::random_entangled_state(class_state_finite &state, StateInitType type, [[maybe_unused]] const std::string &sector, long chi_lim,
                                                          bool use_eigenspinors) {
    if(use_eigenspinors)
        set_random_entangled_state_with_random_spinors(state, type, chi_lim);
    else
        set_random_entangled_state_with_random_spinors(state, type, chi_lim);
}

void tools::finite::mps::internal::set_random_entangled_state_with_random_spinors(class_state_finite &state, StateInitType type, long chi_lim) {
    tools::log->info("Setting random entangled state with random unit spinors...");
    const auto  spin_dim        = state.get_mps_site<size_t>(0).spin_dim();
    auto        bond_dimensions = internal::get_valid_bond_dimensions(state.get_length() + 1, spin_dim, chi_lim);
    bool        pastCenter      = false;
    std::string label           = "A";
    for(auto &mps_ptr : state.mps_sites) {
        auto &          mps  = *mps_ptr;
        auto            chiL = bond_dimensions[mps.get_position()];
        auto            chiR = bond_dimensions[mps.get_position() + 1];
        auto            size = spin_dim * chiL * chiR;
        auto            smdt = pastCenter ? chiR : chiL;
        Eigen::VectorXd Ltmp = Eigen::VectorXd(smdt).unaryExpr([]([[maybe_unused]] auto dummy) { return rnd::uniform_double_01(); });
        std::sort(Ltmp.data(), Ltmp.data() + Ltmp.size(), std::greater<double>());
        Eigen::Tensor<Scalar, 1> L = Textra::TensorCast(Ltmp.normalized()).cast<Scalar>();
        Eigen::Tensor<Scalar, 3> G(spin_dim, chiL, chiR);
        if(type == StateInitType::REAL) {
            Eigen::VectorXd Gtmp = Eigen::VectorXd(size).unaryExpr([]([[maybe_unused]] auto dummy) { return rnd::uniform_double_box(-1.0, 1.0); });
            G                    = Textra::TensorCast(Gtmp.normalized(), spin_dim, chiL, chiR).cast<Scalar>();
        } else if(type == StateInitType::CPLX) {
            Eigen::VectorXcd Gtmp = Eigen::VectorXcd(size).unaryExpr([]([[maybe_unused]] auto dummy) { return rnd::uniform_complex_in_unit_circle(); });
            G                     = Textra::TensorCast(Gtmp.normalized(), spin_dim, chiL, chiR).cast<Scalar>();
        }

        mps.set_mps(G, L, 0, label);
        if(mps.isCenter()) {
            Ltmp = Eigen::VectorXd(chiR).unaryExpr([]([[maybe_unused]] auto dummy) { return rnd::uniform_double_01(); });
            std::sort(Ltmp.data(), Ltmp.data() + Ltmp.size(), std::greater<double>());
            Eigen::Tensor<Scalar, 1> LC = Textra::TensorCast(Ltmp.normalized()).cast<Scalar>();
            mps.set_LC(LC);
            pastCenter = true;
            label      = "B";
        }
    }
    tools::log->info("Setting random entangled state with random unit spinors... OK");
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
}

void tools::finite::mps::internal::set_random_entangled_state_in_sector_using_eigenspinors(class_state_finite &state, StateInitType type,
                                                                                           const std::string &sector, long chi_lim) {
    const auto spin_dim        = state.get_mps_site<size_t>(0).spin_dim();
    auto       bond_dimensions = internal::get_valid_bond_dimensions(state.get_length() + 1, spin_dim, chi_lim);
    auto       axis            = internal::get_axis(sector);
    auto       sign            = internal::get_sign(sector);
    tools::log->info("Setting random entangled state in sector {} using eigenspinors of the pauli matrix σ{}...", sector, axis);
    tools::log->info("Target bond dimensions: {}", bond_dimensions);
    if(type == StateInitType::REAL and axis == "y") throw std::runtime_error("StateInitType REAL incompatible with state in sector [y] which impliex CPLX");
    bool        past_center = false;
    std::string label       = "A";
    for(auto &mps_ptr : state.mps_sites) {
        auto &                                                            mps  = *mps_ptr;
        auto                                                              chiL = bond_dimensions[mps.get_position()];
        auto                                                              chiR = bond_dimensions[mps.get_position() + 1];
        auto                                                              rows = past_center ? chiL : chiL * spin_dim;
        auto                                                              cols = past_center ? chiR * spin_dim : chiR;
        Eigen::Tensor<Scalar, 1>                                          L = Textra::TensorCast(Eigen::VectorXd::Ones(chiL).normalized()).cast<Scalar>();
        Eigen::Tensor<Scalar, 3>                                          G(spin_dim, chiL, chiR);
        Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> G_mat(G.data(), rows, cols);
        // Here we construct the set of spinors for each site using randomly selected eigenspinors
        Eigen::array<long, 3> offset3;
        Eigen::array<long, 3> extent3{spin_dim, 1, 1};
        Eigen::array<long, 1> extent1{spin_dim};
        for(long col = 0; col < chiR; col++) {
            for(long row = 0; row < chiL; row++) {
                auto rnd_sign                              = 2 * rnd::uniform_integer_01() - 1;
                auto spinor                                = internal::get_spinor(axis, rnd_sign);
                offset3                                    = {0, row, col};
                G.slice(offset3, extent3).reshape(extent1) = Textra::TensorMap(spinor);
            }
        }
        G_mat.colwise().normalize();

        //        Textra::normalize(G);
        //        G = Textra::TensorMap(Textra::VectorCast(G).normalized(), spin_dim, chiL, chiR);
        mps.set_mps(G, L, 0, label);
        if(mps.isCenter()) {
            Eigen::Tensor<Scalar, 1> LC = Textra::TensorCast(Eigen::VectorXd::Ones(chiR).normalized()).cast<Scalar>();
            mps.set_LC(LC);
            past_center = true;
            label       = "B";
        }
    }

    // Here we can make sure the state ends up in the correct spin parity sector by flipping the last
    // spin if necessary
    auto spin_component = tools::finite::measure::spin_component(state, axis);
    if(spin_component * sign < 0) {
        auto &                                                            mps  = *state.mps_sites.back();
        auto                                                              chiL = bond_dimensions[mps.get_position()];
        auto                                                              chiR = bond_dimensions[mps.get_position() + 1];
        auto                                                              rows = past_center ? chiL : chiL * spin_dim;
        auto                                                              cols = past_center ? chiR * spin_dim : chiR;
        Eigen::array<long, 3>                                             offset3;
        Eigen::array<long, 3>                                             extent3{spin_dim, 1, 1};
        Eigen::array<long, 1>                                             extent1{spin_dim};
        Eigen::Tensor<Scalar, 3>                                          G(spin_dim, chiL, chiR);
        Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> G_mat(G.data(), rows, cols);
        for(long col = 0; col < chiR; col++) {
            for(long row = 0; row < chiL; row++) {
                offset3                                    = {0, row, col};
                Eigen::Tensor<Scalar, 1> old_spinor_tensor = mps.get_M().slice(offset3, extent3).reshape(extent1);
                Eigen::VectorXcd         old_spinor_vector = Textra::VectorCast(old_spinor_tensor);
                Scalar                   old_spinor_expval = old_spinor_vector.dot(internal::get_pauli(axis) * old_spinor_vector);
                int                      old_spinor_sign   = std::real(old_spinor_expval) > 0 ? 1 : -1;
                auto                     new_spinor_vector = get_spinor(axis, -old_spinor_sign);
                G.slice(offset3, extent3).reshape(extent1) = Textra::TensorMap(new_spinor_vector);
            }
        }
        G_mat.colwise().normalize();
        mps.set_M(G);
    }
    if(spin_component * sign < 0) throw std::logic_error("Could not initialize_state in the correct sector");
    state.clear_measurements();
    state.clear_cache();
    state.tag_all_sites_normalized(false); // This operation denormalizes all sites
}

void tools::finite::mps::internal::randomize_given_state(class_state_finite &state, StateInitType type) {
    using namespace qm::spinHalf;
    switch(type) {
        case StateInitType::REAL: tools::finite::mps::apply_random_paulis(state, std::vector<Eigen::Matrix2cd>{id, 1.0 / std::sqrt(2.0) * (sx + sz)}); break;
        case StateInitType::CPLX:
            tools::finite::mps::apply_random_paulis(state, std::vector<Eigen::Matrix2cd>{id, 1.0 / std::sqrt(3.0) * (sx + sy + sz)});
            break;
    }
}