//
// Created by david on 2019-02-01.
//

#include "measure.h"
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>
#include <math/linalg/tensor.h>
#include <math/num.h>
#include <physics/nmspc_quantum_mechanics.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/contraction.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

using namespace std;
using namespace Textra;
using Scalar = class_state_finite::Scalar;

void tools::finite::measure::do_all_measurements(const class_tensors_finite &tensors) {
    tensors.measurements.length                   = measure::length(tensors);
    tensors.measurements.energy                   = measure::energy(tensors); // This number is needed for variance calculation!
    tensors.measurements.energy_per_site          = measure::energy_per_site(tensors);
    tensors.measurements.energy_variance          = measure::energy_variance(tensors);
    tensors.measurements.energy_variance_per_site = measure::energy_variance_per_site(tensors);
}

void tools::finite::measure::do_all_measurements(const class_state_finite &state) {
    state.measurements.length                        = measure::length(state);
    state.measurements.norm                          = measure::norm(state);
    state.measurements.bond_dimension_current        = measure::bond_dimension_current(state);
    state.measurements.bond_dimension_midchain       = measure::bond_dimension_midchain(state);
    state.measurements.bond_dimensions               = measure::bond_dimensions(state);
    state.measurements.truncation_errors             = measure::truncation_errors(state);
    state.measurements.entanglement_entropy_current  = measure::entanglement_entropy_current(state);
    state.measurements.entanglement_entropy_midchain = measure::entanglement_entropy_midchain(state);
    state.measurements.entanglement_entropies        = measure::entanglement_entropies(state);
    state.measurements.number_entropy_current        = measure::number_entropy_current(state);
    state.measurements.number_entropy_midchain       = measure::number_entropy_midchain(state);
    state.measurements.number_entropies              = measure::number_entropies(state);
    state.measurements.renyi_2                       = measure::renyi_entropies(state, 2);
    state.measurements.renyi_3                       = measure::renyi_entropies(state, 3);
    state.measurements.renyi_4                       = measure::renyi_entropies(state, 4);
    state.measurements.renyi_100                     = measure::renyi_entropies(state, 100);
    state.measurements.spin_components               = measure::spin_components(state);
}

size_t tools::finite::measure::length(const class_tensors_finite &tensors) { return tensors.get_length(); }
size_t tools::finite::measure::length(const class_state_finite &state) { return state.get_length(); }

double tools::finite::measure::norm(const class_state_finite &state) {
    if(state.measurements.norm) return state.measurements.norm.value();
    double norm;
    if(state.is_normalized_on_all_sites()) {
        // We know the all sites are normalized. We can check that the current position is normalized
        const auto  pos = std::clamp(state.get_position<long>(), 0l, state.get_length<long>());
        const auto &mps = state.get_mps_site(pos);
        tools::log->trace("Measuring norm using site {} with dimensions {}", pos, mps.dimensions());
        Eigen::Tensor<Scalar, 0> MM = mps.get_M().contract(mps.get_M().conjugate(), Textra::idx({0, 1, 2}, {0, 1, 2}));
        norm                        = std::real(MM(0));

    } else if(state.is_normalized_on_non_active_sites() and not state.active_sites.empty()) {
        tools::log->trace("Measuring norm using active sites {}", state.active_sites);
        Eigen::Tensor<Scalar, 2> chain;
        Eigen::Tensor<Scalar, 2> temp;
        bool                     first = true;
        for(const auto &pos : state.active_sites) {
            const Eigen::Tensor<Scalar, 3> &M = state.get_mps_site(pos).get_M();
            if(first) {
                chain = M.contract(M.conjugate(), idx({0, 1}, {0, 1}));
                first = false;
                continue;
            }
            temp.resize(Textra::array2{M.dimension(2), M.dimension(2)});
            temp.device(Textra::omp::getDevice()) = chain.contract(M, idx({0}, {1})).contract(M.conjugate(), idx({0, 1}, {1, 0}));

            chain = temp;
        }
        norm = std::abs(Textra::MatrixMap(chain).trace());
    } else {
        tools::log->trace("Measuring norm on full chain");
        Eigen::Tensor<Scalar, 2> chain;
        Eigen::Tensor<Scalar, 2> temp;
        bool                     first = true;
        for(const auto &mps : state.mps_sites) {
            const auto &M = mps->get_M();
            if(first) {
                chain = M.contract(M.conjugate(), idx({0, 1}, {0, 1}));
                first = false;
                continue;
            }
            temp.resize(Textra::array2{M.dimension(2), M.dimension(2)});
            temp.device(Textra::omp::getDevice()) = chain.contract(M, idx({0}, {1})).contract(M.conjugate(), idx({0, 1}, {1, 0}));

            chain = temp;
        }
        norm = std::abs(Textra::MatrixMap(chain).trace());
    }

    if(std::abs(norm - 1.0) > settings::precision::max_norm_error) tools::log->debug("Norm far from unity: {:.16f}", norm);
    state.measurements.norm = norm;
    return state.measurements.norm.value();
}

long tools::finite::measure::bond_dimension_current(const class_state_finite &state) {
    if(state.measurements.bond_dimension_current) return state.measurements.bond_dimension_current.value();
    if(state.has_center_point())
        state.measurements.bond_dimension_current = state.current_bond().dimension(0);
    else
        state.measurements.bond_dimension_current = 1;
    return state.measurements.bond_dimension_current.value();
}

long tools::finite::measure::bond_dimension_midchain(const class_state_finite &state) {
    if(state.measurements.bond_dimension_midchain) return state.measurements.bond_dimension_midchain.value();
    state.measurements.bond_dimension_midchain = state.midchain_bond().dimension(0);
    return state.measurements.bond_dimension_midchain.value();
}

std::vector<long> tools::finite::measure::bond_dimensions(const class_state_finite &state) {
    if(state.measurements.bond_dimensions) return state.measurements.bond_dimensions.value();
    std::vector<long> bond_dimensions;
    bond_dimensions.reserve(state.get_length() + 1);
    if(not state.has_center_point()) bond_dimensions.emplace_back(state.mps_sites.front()->get_chiL());
    for(const auto &mps : state.mps_sites) {
        bond_dimensions.emplace_back(mps->get_L().dimension(0));
        if(mps->isCenter()) { bond_dimensions.emplace_back(mps->get_LC().dimension(0)); }
    }
    if(bond_dimensions.size() != state.get_length() + 1) throw std::logic_error("bond_dimensions.size() should be length+1");
    state.measurements.bond_dimensions = bond_dimensions;
    return state.measurements.bond_dimensions.value();
}

std::vector<long> tools::finite::measure::bond_dimensions_merged(const class_state_finite &state) {
    std::vector<long> bond_dimensions;
    if(state.active_sites.empty()) return bond_dimensions;

    for(const auto &pos : state.active_sites) {
        bond_dimensions.emplace_back(state.get_mps_site(pos).get_L().dimension(0));
        if(state.get_mps_site(pos).isCenter()) { bond_dimensions.emplace_back(state.get_mps_site(pos).get_LC().dimension(0)); }
    }
    if(state.active_sites.size() > 1) {
        bond_dimensions.pop_back();
        bond_dimensions.erase(bond_dimensions.begin());
    } else if(state.get_direction() == 1)
        bond_dimensions.pop_back();
    return bond_dimensions;
}

double tools::finite::measure::entanglement_entropy_current(const class_state_finite &state) {
    if(state.measurements.entanglement_entropy_current) return state.measurements.entanglement_entropy_current.value();
    auto t_ent = tools::common::profile::get_default_prof()["t_ent"]->tic_token();
    if(state.has_center_point()) {
        auto &                   LC                     = state.current_bond();
        Eigen::Tensor<Scalar, 0> SE                     = -LC.square().contract(LC.square().log().eval(), idx({0}, {0}));
        state.measurements.entanglement_entropy_current = std::abs(SE(0));
    } else
        state.measurements.entanglement_entropy_current = 0;
    return state.measurements.entanglement_entropy_current.value();
}

double tools::finite::measure::entanglement_entropy_midchain(const class_state_finite &state) {
    if(state.measurements.entanglement_entropy_midchain) return state.measurements.entanglement_entropy_midchain.value();
    auto                     t_ent                   = tools::common::profile::get_default_prof()["t_ent"]->tic_token();
    auto &                   LC                      = state.midchain_bond();
    Eigen::Tensor<Scalar, 0> SE                      = -LC.square().contract(LC.square().log().eval(), idx({0}, {0}));
    state.measurements.entanglement_entropy_midchain = std::abs(SE(0));
    return state.measurements.entanglement_entropy_midchain.value();
}

std::vector<double> tools::finite::measure::entanglement_entropies(const class_state_finite &state) {
    if(state.measurements.entanglement_entropies) return state.measurements.entanglement_entropies.value();
    auto                t_ent = tools::common::profile::get_default_prof()["t_ent"]->tic_token();
    std::vector<double> entanglement_entropies;
    entanglement_entropies.reserve(state.get_length() + 1);
    if(not state.has_center_point()) entanglement_entropies.emplace_back(0);
    for(const auto &mps : state.mps_sites) {
        auto &                   L  = mps->get_L();
        Eigen::Tensor<Scalar, 0> SE = -L.square().contract(L.square().log().eval(), idx({0}, {0}));
        entanglement_entropies.emplace_back(std::abs(SE(0)));
        if(mps->isCenter()) {
            auto &LC = mps->get_LC();
            SE       = -LC.square().contract(LC.square().log().eval(), idx({0}, {0}));
            entanglement_entropies.emplace_back(std::abs(SE(0)));
            state.measurements.entanglement_entropy_current = std::abs(SE(0));
        }
    }
    if(entanglement_entropies.size() != state.get_length() + 1) throw std::logic_error("entanglement_entropies.size() should be length+1");
    if(entanglement_entropies.front() != 0.0) throw std::logic_error(fmt::format("First entropy should be 0. Got: {:.16f}", entanglement_entropies.front()));
    if(entanglement_entropies.back() != 0.0) throw std::logic_error(fmt::format("Last entropy should be 0. Got: {:.16f}", entanglement_entropies.back()));
    state.measurements.entanglement_entropies = entanglement_entropies;
    return state.measurements.entanglement_entropies.value();
}

std::vector<double> tools::finite::measure::renyi_entropies(const class_state_finite &state, double q) {
    if(q == 1.0) return entanglement_entropies(state);
    if(q == 2.0 and state.measurements.renyi_2) return state.measurements.renyi_2.value();
    if(q == 3.0 and state.measurements.renyi_3) return state.measurements.renyi_3.value();
    if(q == 4.0 and state.measurements.renyi_4) return state.measurements.renyi_4.value();
    if(q == 100.0 and state.measurements.renyi_100) return state.measurements.renyi_100.value();
    auto                t_ent = tools::common::profile::get_default_prof()["t_ent"]->tic_token();
    std::vector<double> renyi_q;
    renyi_q.reserve(state.get_length() + 1);
    if(not state.has_center_point()) renyi_q.emplace_back(0);
    for(const auto &mps : state.mps_sites) {
        const auto &             L = mps->get_L();
        Eigen::Tensor<Scalar, 0> RE;
        RE = (1.0 / 1.0 - q) * L.pow(2.0 * q).sum().log();
        renyi_q.emplace_back(std::abs(RE(0)));
        if(mps->isCenter()) {
            const auto &LC = mps->get_LC();
            RE             = (1.0 / 1.0 - q) * LC.pow(2.0 * q).sum().log();
            renyi_q.emplace_back(std::abs(RE(0)));
        }
    }
    if(renyi_q.size() != state.get_length() + 1) throw std::logic_error("renyi_q.size() should be length+1");
    if(q == 2.0) {
        state.measurements.renyi_2 = renyi_q;
        return state.measurements.renyi_2.value();
    }
    if(q == 3.0) {
        state.measurements.renyi_3 = renyi_q;
        return state.measurements.renyi_3.value();
    }
    if(q == 4.0) {
        state.measurements.renyi_4 = renyi_q;
        return state.measurements.renyi_4.value();
    }
    if(q == 100.0) {
        state.measurements.renyi_100 = renyi_q;
        return state.measurements.renyi_100.value();
    }
    return renyi_q;
}

std::array<double, 3> tools::finite::measure::spin_components(const class_state_finite &state) {
    if(state.measurements.spin_components) return state.measurements.spin_components.value();
    double spin_x                      = measure::spin_component(state, qm::spinHalf::sx);
    double spin_y                      = measure::spin_component(state, qm::spinHalf::sy);
    double spin_z                      = measure::spin_component(state, qm::spinHalf::sz);
    state.measurements.spin_components = {spin_x, spin_y, spin_z};
    return state.measurements.spin_components.value();
}

double tools::finite::measure::spin_component(const class_state_finite &state, const Eigen::Matrix2cd &paulimatrix) {
    auto t_spn       = tools::common::profile::get_default_prof()["t_spn"]->tic_token();
    auto [mpo, L, R] = qm::mpo::pauli_mpo(paulimatrix);
    Eigen::Tensor<Scalar, 3> temp;
    for(const auto &mps : state.mps_sites) {
        const auto &M = mps->get_M();
        temp.resize(M.dimension(2), M.dimension(2), mpo.dimension(1));
        temp.device(Textra::omp::getDevice()) = L.contract(M, idx({0}, {1})).contract(M.conjugate(), idx({0}, {1})).contract(mpo, idx({0, 1, 3}, {0, 2, 3}));
        L                                     = temp;
    }

    if(L.dimensions() != R.dimensions()) throw std::runtime_error("spin_component(): L and R dimension mismatch");
    Eigen::Tensor<Scalar, 0> spin_tmp = L.contract(R, idx({0, 1, 2}, {0, 1, 2}));
    double                   spin     = std::real(spin_tmp(0));
    return spin;
}

double tools::finite::measure::spin_component(const class_state_finite &state, const std::string &axis) {
    if(axis.find('x') != std::string::npos) return measure::spin_component(state, qm::spinHalf::sx);
    if(axis.find('y') != std::string::npos) return measure::spin_component(state, qm::spinHalf::sy);
    if(axis.find('z') != std::string::npos) return measure::spin_component(state, qm::spinHalf::sz);
    throw std::runtime_error("Unexpected axis: " + axis);
}

std::vector<double> tools::finite::measure::truncation_errors(const class_state_finite &state) {
    if(state.measurements.truncation_errors) return state.measurements.truncation_errors.value();
    std::vector<double> truncation_errors;
    if(not state.has_center_point()) truncation_errors.emplace_back(0);
    for(const auto &mps : state.mps_sites) {
        truncation_errors.emplace_back(mps->get_truncation_error());
        if(mps->isCenter()) truncation_errors.emplace_back(mps->get_truncation_error_LC());
    }
    if(truncation_errors.size() != state.get_length() + 1) throw std::logic_error("truncation_errors.size() should be length+1");
    state.measurements.truncation_errors = truncation_errors;
    return state.measurements.truncation_errors.value();
}

std::vector<double> tools::finite::measure::truncation_errors_active(const class_state_finite &state) {
    std::vector<double> truncation_errors;
    for(const auto &site : state.active_sites) {
        const auto &mps = state.get_mps_site(site);
        truncation_errors.emplace_back(mps.get_truncation_error());
        if(mps.isCenter()) truncation_errors.emplace_back(mps.get_truncation_error_LC());
    }
    return truncation_errors;
}

Eigen::Tensor<Scalar, 1> tools::finite::measure::mps_wavefn(const class_state_finite &state) {
    Eigen::Tensor<Scalar, 2> temp;
    Eigen::Tensor<Scalar, 2> chain(1, 1);
    chain.setConstant(1.0);
    // The "chain" is a matrix whose 0 index keeps growing.
    // For each site that passes, it grows by GA.dimension(0) = phys dim
    // Say the state is a 16x7 matrix (having contracted 4 particles, and the latest
    // chi was 7). Then contracting the next site, with dimensions 2x7x9 will get you a
    // 16x2x9 tensor. Now the reshaping convert it into a 32 x 9 matrix. Because
    // Eigen is column major, the doubling 16->32 will stack the third index twice.

    for(auto &mps : state.mps_sites) {
        long dim0 = mps->spin_dim();
        long dimR = mps->get_chiR();
        long dimL = chain.dimension(0);
        temp      = chain.contract(mps->get_M(), idx({1}, {1})).reshape(array2{dimL * dim0, dimR});
        chain     = temp;
    }

    Eigen::Tensor<Scalar, 1> mps_chain  = chain.reshape(array1{chain.dimension(0)});
    double                   norm_chain = Textra::VectorMap(chain).norm();
    if(std::abs(norm_chain - 1.0) > settings::precision::max_norm_error) {
        tools::log->warn("Norm far from unity: {}", norm_chain);
        throw std::runtime_error("Norm too far from unity: " + std::to_string(norm_chain));
    }
    return mps_chain;
}

template<typename state_or_mps_type>
double tools::finite::measure::energy_minus_energy_reduced(const state_or_mps_type &state, const class_model_finite &model, const class_edges_finite &edges,
                                                           tensors_measure_finite *measurements) {
    if(measurements != nullptr and measurements->energy_minus_energy_reduced) return measurements->energy_minus_energy_reduced.value();
    if constexpr(std::is_same_v<state_or_mps_type, class_state_finite>) {
        if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
            throw std::runtime_error(fmt::format("Could not compute energy: active sites are not equal: state {} | model {} | edges {}", state.active_sites,
                                                 model.active_sites, edges.active_sites));
        return tools::finite::measure::energy_minus_energy_reduced(state.get_multisite_mps(), model, edges, measurements);
    } else {
        const auto &mpo = model.get_multisite_mpo();
        const auto &env = edges.get_multisite_ene_blk();
        tools::log->trace("Measuring energy: state dims {} | model sites {} dims {} | edges sites {} dims [L{} R{}]", state.dimensions(), model.active_sites,
                          mpo.dimensions(), edges.active_sites, env.L.dimensions(), env.R.dimensions());
        auto   t_ene        = tools::common::profile::get_default_prof()["t_ene"]->tic_token();
        double e_minus_ered = tools::common::contraction::expectation_value(state, mpo, env.L, env.R);
        if(measurements != nullptr) measurements->energy_minus_energy_reduced = e_minus_ered;
        return e_minus_ered;
    }
}

template double tools::finite::measure::energy_minus_energy_reduced(const class_state_finite &, const class_model_finite &model,
                                                                    const class_edges_finite &edges, tensors_measure_finite *measurements);
template double tools::finite::measure::energy_minus_energy_reduced(const Eigen::Tensor<Scalar, 3> &, const class_model_finite &model,
                                                                    const class_edges_finite &edges, tensors_measure_finite *measurements);

template<typename state_or_mps_type>
double tools::finite::measure::energy(const state_or_mps_type &state, const class_model_finite &model, const class_edges_finite &edges,
                                      tensors_measure_finite *measurements) {
    if(measurements != nullptr and measurements->energy) return measurements->energy.value();
    // This measures the actual energy of the system regardless of the reduced/non-reduced state of the MPO's
    // If they are reduced, then
    //      "Actual energy" = (E - E_reduced) + E_reduced = (~0) + E_reduced = E
    // Else
    //      "Actual energy" = (E - E_reduced) + E_reduced = E  + 0 = E
    double energy;
    if constexpr(std::is_same_v<state_or_mps_type, class_state_finite>)
        energy = tools::finite::measure::energy_minus_energy_reduced(state.get_multisite_mps(), model, edges, measurements) + model.get_energy_reduced();
    else
        energy = tools::finite::measure::energy_minus_energy_reduced(state, model, edges, measurements) + model.get_energy_reduced();

    if(measurements != nullptr) measurements->energy = energy;
    return energy;
}

template double tools::finite::measure::energy(const class_state_finite &, const class_model_finite &model, const class_edges_finite &edges,
                                               tensors_measure_finite *measurements);
template double tools::finite::measure::energy(const Eigen::Tensor<Scalar, 3> &, const class_model_finite &model, const class_edges_finite &edges,
                                               tensors_measure_finite *measurements);

template<typename state_or_mps_type>
double tools::finite::measure::energy_per_site(const state_or_mps_type &state, const class_model_finite &model, const class_edges_finite &edges,
                                               tensors_measure_finite *measurements) {
    double energy_per_site = tools::finite::measure::energy(state, model, edges, measurements) / static_cast<double>(model.get_length());
    if(measurements != nullptr) measurements->energy_per_site = energy_per_site;
    return energy_per_site;
}

template double tools::finite::measure::energy_per_site(const class_state_finite &, const class_model_finite &model, const class_edges_finite &edges,
                                                        tensors_measure_finite *measurements);
template double tools::finite::measure::energy_per_site(const Eigen::Tensor<Scalar, 3> &, const class_model_finite &model, const class_edges_finite &edges,
                                                        tensors_measure_finite *measurements);

template<typename state_or_mps_type>
double tools::finite::measure::energy_variance(const state_or_mps_type &state, const class_model_finite &model, const class_edges_finite &edges,
                                               tensors_measure_finite *measurements) {
    // Here we show that the variance calculated with reduced-energy mpo's is equivalent to the usual way.
    // If mpo's are reduced:
    //      Var H = <(H-E_red)²> - <H-E_red>²     = <H²>  - 2<H>E_red + E_red² - (<H> - E_red)²
    //                                            = H2    - 2*E*E_red + E_red² - E² + 2*E*E_red - E_red²
    //                                            = H2    - E²
    //      Note that in the last line, H2-E² is a subtraction of two large numbers --> catastrophic cancellation --> loss of precision.
    //      On the other hand Var H = <(H-E_red)²> - energy_minus_energy_reduced² = <(H-E_red)²> - ~dE², where both terms are always  << 1.
    //      The first term computed from a double-layer of reduced mpo's.
    //      In the second term dE is usually small, in fact identically zero immediately after an energy-reduction operation,
    //      but may grow if the optimization steps make significant progress refining E.
    // Else, if E_red = 0 (i.e. not reduced) we get the usual formula:
    //      Var H = <(H - 0)²> - <H - 0>² = H2 - E²
    if(measurements != nullptr and measurements->energy_variance) return measurements->energy_variance.value();

    if constexpr(std::is_same_v<state_or_mps_type, class_state_finite>) {
        if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
            throw std::runtime_error(fmt::format("Could not compute energy variance: active sites are not equal: state {} | model {} | edges {}",
                                                 state.active_sites, model.active_sites, edges.active_sites));
        if(state.active_sites.empty()) throw std::runtime_error("Could not compute energy variance: active sites are empty");
        return tools::finite::measure::energy_variance(state.get_multisite_mps(), model, edges, measurements);
    } else {
        double energy = 0;
        if(model.is_reduced())
            energy = tools::finite::measure::energy_minus_energy_reduced(state, model, edges, measurements);
        else
            energy = tools::finite::measure::energy(state, model, edges, measurements);
        double E2 = energy * energy;

        if(not num::all_equal(model.active_sites, edges.active_sites))
            throw std::runtime_error(
                fmt::format("Could not compute energy variance: active sites are not equal: model {} | edges {}", model.active_sites, edges.active_sites));

        const auto &mpo = model.get_multisite_mpo_squared();
        const auto &env = edges.get_multisite_var_blk();
        tools::log->trace("Measuring energy variance: state dims {} | model sites {} dims {} | edges sites {} dims [L{} R{}]", state.dimensions(),
                          model.active_sites, mpo.dimensions(), edges.active_sites, env.L.dimensions(), env.R.dimensions());

        if(state.dimension(0) != mpo.dimension(2))
            throw std::runtime_error(
                fmt::format("State and model have incompatible physical dimension: state dim {} | model dim {}", state.dimension(0), mpo.dimension(2)));
        auto   t_var = tools::common::profile::get_default_prof()["t_var"]->tic_token();
        double H2    = tools::common::contraction::expectation_value(state, mpo, env.L, env.R);
        double var   = std::abs(H2 - E2);
        if(measurements != nullptr) measurements->energy_variance = var;
        return var;
    }
}

template double tools::finite::measure::energy_variance(const class_state_finite &, const class_model_finite &model, const class_edges_finite &edges,
                                                        tensors_measure_finite *measurements);
template double tools::finite::measure::energy_variance(const Eigen::Tensor<Scalar, 3> &, const class_model_finite &model, const class_edges_finite &edges,
                                                        tensors_measure_finite *measurements);

template<typename state_or_mps_type>
double tools::finite::measure::energy_variance_per_site(const state_or_mps_type &state, const class_model_finite &model, const class_edges_finite &edges,
                                                        tensors_measure_finite *measurements) {
    double energy_variance_per_site = tools::finite::measure::energy_variance(state, model, edges, measurements) / static_cast<double>(model.get_length());
    if(measurements != nullptr) measurements->energy_variance_per_site = energy_variance_per_site;
    return energy_variance_per_site;
}

template double tools::finite::measure::energy_variance_per_site(const class_state_finite &, const class_model_finite &model, const class_edges_finite &edges,
                                                                 tensors_measure_finite *measurements);
template double tools::finite::measure::energy_variance_per_site(const Eigen::Tensor<Scalar, 3> &, const class_model_finite &model,
                                                                 const class_edges_finite &edges, tensors_measure_finite *measurements);

template<typename state_or_mps_type>
double tools::finite::measure::energy_normalized(const state_or_mps_type &state, const class_model_finite &model, const class_edges_finite &edges,
                                                 double energy_min_per_site, double energy_max_per_site, tensors_measure_finite *measurements) {
    return (tools::finite::measure::energy_per_site(state, model, edges, measurements) - energy_min_per_site) / (energy_max_per_site - energy_min_per_site);
}

template double tools::finite::measure::energy_normalized(const class_state_finite &, const class_model_finite &model, const class_edges_finite &edges, double,
                                                          double, tensors_measure_finite *measurements);
template double tools::finite::measure::energy_normalized(const Eigen::Tensor<Scalar, 3> &, const class_model_finite &model, const class_edges_finite &edges,
                                                          double, double, tensors_measure_finite *measurements);

extern double tools::finite::measure::energy_reduced(const class_tensors_finite &tensors) { return tensors.model->get_energy_reduced(); }
extern double tools::finite::measure::energy_per_site_reduced(const class_tensors_finite &tensors) { return tensors.model->get_energy_per_site_reduced(); }

double tools::finite::measure::energy_minus_energy_reduced(const class_tensors_finite &tensors) {
    tensors.assert_edges_ene();
    return energy_minus_energy_reduced(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
}

double tools::finite::measure::energy(const class_tensors_finite &tensors) {
    if(tensors.measurements.energy)
        return tensors.measurements.energy.value();
    else {
        tensors.assert_edges_ene();
        tensors.measurements.energy = tools::finite::measure::energy(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
        return tensors.measurements.energy.value();
    }
}

double tools::finite::measure::energy_per_site(const class_tensors_finite &tensors) {
    if(tensors.measurements.energy_per_site)
        return tensors.measurements.energy_per_site.value();
    else {
        tensors.assert_edges_ene();
        tensors.measurements.energy_per_site = tools::finite::measure::energy_per_site(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
        return tensors.measurements.energy_per_site.value();
    }
}

double tools::finite::measure::energy_variance(const class_tensors_finite &tensors) {
    if(tensors.measurements.energy_variance)
        return tensors.measurements.energy_variance.value();
    else {
        tensors.assert_edges_var();
        tensors.measurements.energy_variance = tools::finite::measure::energy_variance(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
        return tensors.measurements.energy_variance.value();
    }
}

double tools::finite::measure::energy_variance_per_site(const class_tensors_finite &tensors) {
    if(tensors.measurements.energy_variance_per_site)
        return tensors.measurements.energy_variance_per_site.value();
    else {
        tensors.assert_edges_var();
        tensors.measurements.energy_variance_per_site =
            tools::finite::measure::energy_variance_per_site(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
        return tensors.measurements.energy_variance_per_site.value();
    }
}

double tools::finite::measure::energy_normalized(const class_tensors_finite &tensors, double emin, double emax) {
    tensors.assert_edges_ene();
    return tools::finite::measure::energy_normalized(*tensors.state, *tensors.model, *tensors.edges, emin, emax);
}


double tools::finite::measure::energy_minus_energy_reduced(const class_state_finite & state, const class_tensors_finite &tensors) {
    return tools::finite::measure::energy_minus_energy_reduced(state, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy(const class_state_finite & state, const class_tensors_finite &tensors) {
    return tools::finite::measure::energy(state, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_per_site(const class_state_finite & state, const class_tensors_finite &tensors) {
    return tools::finite::measure::energy_per_site(state, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_variance(const class_state_finite & state, const class_tensors_finite &tensors) {
    return tools::finite::measure::energy_variance(state, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_variance_per_site(const class_state_finite & state, const class_tensors_finite &tensors) {
    return tools::finite::measure::energy_variance_per_site(state, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_normalized(const class_state_finite & state, const class_tensors_finite &tensors, double emin, double emax) {
    return tools::finite::measure::energy_normalized(state, *tensors.model, *tensors.edges, emin, emax);
}



double tools::finite::measure::energy_minus_energy_reduced(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_finite &tensors) {
    return tools::finite::measure::energy_minus_energy_reduced(mps, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_finite &tensors) {
    return tools::finite::measure::energy(mps, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_per_site(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_finite &tensors) {
    return tools::finite::measure::energy_per_site(mps, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_variance(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_finite &tensors) {
    return tools::finite::measure::energy_variance(mps, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_variance_per_site(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_finite &tensors) {
    return tools::finite::measure::energy_variance_per_site(mps, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_normalized(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_finite &tensors, double emin, double emax) {
    return tools::finite::measure::energy_normalized(mps, *tensors.model, *tensors.edges, emin, emax);
}
