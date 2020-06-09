//
// Created by david on 2020-06-07.
//

#include "tools/finite/mps.h"
#include <bitset>
#include <config/nmspc_settings.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <math/rnd.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
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

void tools::finite::mps::internal::random_entangled_state(class_state_finite &state, const std::string &sector, long chi_lim, bool use_eigenspinors,
                                                          bool real) {
    if(use_eigenspinors) set_random_entangled_state_with_spinors_in_c2(state, chi_lim, real);
    else
        set_random_entangled_state_with_spinors_in_c2(state, chi_lim, real);
}

void tools::finite::mps::internal::set_random_entangled_state_with_spinors_in_c2(class_state_finite &state, long chi_lim, bool real) {
    tools::log->info("Setting random entangled state with spinors randomly in C2 ...");
    const auto spin_dim        = state.get_mps_site(0).spin_dim();
    auto       bond_dimensions = internal::get_valid_bond_dimensions(state.get_length() + 1, spin_dim, chi_lim);
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps  = *mps_ptr;
        auto  chiL = bond_dimensions[mps.get_position()];
        auto  chiR = bond_dimensions[mps.get_position() + 1];
        auto  size = spin_dim * chiL * chiR;

        Eigen::Tensor<Scalar, 1> L = Textra::MatrixTensorMap(Eigen::VectorXd::Ones(chiL).normalized()).cast<Scalar>();
        Eigen::Tensor<Scalar, 3> G(spin_dim, chiL, chiR);
        if(real) G = Textra::MatrixTensorMap(Eigen::VectorXd::Random(size).normalized(), spin_dim, chiL, chiR).cast<Scalar>();
        else
            G = Textra::MatrixTensorMap(Eigen::VectorXcd::Random(size).normalized(), spin_dim, chiL, chiR);

        mps.set_mps(G, L);
        if(mps.isCenter()){
            Eigen::Tensor<Scalar, 1> LC = Textra::MatrixTensorMap(Eigen::VectorXd::Ones(chiR).normalized()).cast<Scalar>();
            mps.set_LC(LC);
        }
    }
    tools::log->info("Setting random entangled state with spinors randomly in C2 ... OK");
}

void tools::finite::mps::internal::set_random_entangled_state_in_sector_using_eigenspinors(class_state_finite &state, const std::string &sector, long chi_lim) {
    const auto spin_dim        = state.get_mps_site(0).spin_dim();
    auto       bond_dimensions = internal::get_valid_bond_dimensions(state.get_length() + 1, spin_dim, chi_lim);
    auto axis = internal::get_axis(sector);
    auto sign = internal::get_sign(sector);
    tools::log->info("Setting random entangled state in sector {} using eigenspinors of the pauli matrix σ{}...", sector, axis);
    tools::log->info("Target bond dimensions: {}", bond_dimensions);
    bool past_center = false;
    for(auto &mps_ptr : state.mps_sites) {
        auto &mps  = *mps_ptr;
        auto  chiL = bond_dimensions[mps.get_position()];
        auto  chiR = bond_dimensions[mps.get_position() + 1];
        auto  size = spin_dim * chiL * chiR;
        auto  rows = past_center ? chiL : chiL * spin_dim;
        auto  cols = past_center ? chiR * spin_dim : chiR;
        Eigen::Tensor<Scalar, 1> L = Textra::MatrixTensorMap(Eigen::VectorXd::Ones(chiL).normalized()).cast<Scalar>();
        Eigen::Tensor<Scalar, 3> G(spin_dim, chiL, chiR);
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> G_mat (G.data(),rows,cols);
        // Here we construct the set of spinors for each site using randomly selected eigenspinors
        Eigen::array<long, 3> offset3;
        Eigen::array<long, 3> extent3{spin_dim, 1, 1};
        Eigen::array<long, 1> extent1{spin_dim};
        for(long col = 0; col < chiR; col++) {
            for(long row = 0; row < chiL; row++){
                auto rnd_sign                              = 2 * rnd::uniform_integer_01() - 1;
                auto spinor                                = internal::get_spinor(axis, rnd_sign);
                offset3                                    = {0, row, col};
                G.slice(offset3, extent3).reshape(extent1) = Textra::MatrixTensorMap(spinor);
            }
        }
        G_mat.colwise().normalize();

//        Textra::normalize(G);
//        G = Textra::MatrixTensorMap(Textra::Tensor_to_Vector(G).normalized(), spin_dim, chiL, chiR);
        mps.set_mps(G, L);
        if(mps.isCenter()){
            Eigen::Tensor<Scalar, 1> LC = Textra::MatrixTensorMap(Eigen::VectorXd::Ones(chiR).normalized()).cast<Scalar>();
            mps.set_LC(LC);
            past_center = true;
        }
    }

    // Here we can make sure the state ends up in the correct spin parity sector by flipping the last
    // spin if necessary
    auto spin_component = tools::finite::measure::spin_component(state, axis);
    past_center = false;
    if(spin_component * sign < 0) {
        auto &                   mps  = *state.mps_sites.back();
        auto                     chiL = bond_dimensions[mps.get_position()];
        auto                     chiR = bond_dimensions[mps.get_position() + 1];
        auto  rows = past_center ? chiL : chiL * spin_dim;
        auto  cols = past_center ? chiR * spin_dim : chiR;
        Eigen::array<long, 3>    offset3;
        Eigen::array<long, 3>    extent3{spin_dim, 1, 1};
        Eigen::array<long, 1>    extent1{spin_dim};
        Eigen::Tensor<Scalar, 3> G(spin_dim, chiL, chiR);
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> G_mat (G.data(),rows,cols);
        for(long col = 0; col < chiR; col++) {
            for(long row = 0; row < chiL; row++) {
                offset3                                    = {0, row, col};
                Eigen::Tensor<Scalar, 1> old_spinor_tensor = mps.get_M().slice(offset3, extent3).reshape(extent1);
                Eigen::VectorXcd         old_spinor_vector = Textra::Tensor_to_Vector(old_spinor_tensor);
                Scalar                   old_spinor_expval = old_spinor_vector.dot(internal::get_pauli(axis) * old_spinor_vector);
                int                      old_spinor_sign   = std::real(old_spinor_expval) > 0 ? 1 : -1;
                auto                     new_spinor_vector = get_spinor(axis, -old_spinor_sign);
                G.slice(offset3, extent3).reshape(extent1) = Textra::MatrixTensorMap(new_spinor_vector);
            }
        }
        G_mat.colwise().normalize();
        mps.set_M(G);
        if(mps.isCenter()){
            past_center = true;
        }

    }
    if(spin_component * sign < 0) throw std::logic_error("Could not initialize_state in the correct sector");
}


void tools::finite::mps::internal::randomize_given_state(class_state_finite & state) {
    tools::finite::mps::apply_random_paulis(state, {"x","z"});
}