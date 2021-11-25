
#include "../measure.h"
#include <config/settings.h>
#include <io/fmt.h>
#include <math/num.h>
#include <math/tenx.h>
#include <qm/mpo.h>
#include <qm/spin.h>
#include <tensors/edges/EdgesFinite.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateFinite.h>
#include <tensors/TensorsFinite.h>
#include <tid/tid.h>
#include <tools/common/contraction.h>
#include <tools/common/log.h>

using cplx = tools::finite::measure::cplx;
using real = tools::finite::measure::real;

void tools::finite::measure::do_all_measurements(const TensorsFinite &tensors) {
    tensors.measurements.length = measure::length(tensors);

    // No need for energy or variance in fLBIT simulations. In fact, we don't update the ene and var environments,
    // so this operation would give an error
    if(tensors.state->get_algorithm() == AlgorithmType::fLBIT) return;
    tensors.measurements.energy                   = measure::energy(tensors); // This number is needed for variance calculation!
    tensors.measurements.energy_per_site          = measure::energy_per_site(tensors);
    tensors.measurements.energy_variance          = measure::energy_variance(tensors);
    tensors.measurements.energy_variance_per_site = measure::energy_variance_per_site(tensors);
}

void tools::finite::measure::do_all_measurements(const StateFinite &state) {
    state.measurements.length                        = measure::length(state);
    state.measurements.norm                          = measure::norm(state);
    state.measurements.bond_dimension_current        = measure::bond_dimension_current(state);
    state.measurements.bond_dimension_midchain       = measure::bond_dimension_midchain(state);
    state.measurements.bond_dimensions               = measure::bond_dimensions(state);
    state.measurements.truncation_errors             = measure::truncation_errors(state);
    state.measurements.entanglement_entropy_current  = measure::entanglement_entropy_current(state);
    state.measurements.entanglement_entropy_midchain = measure::entanglement_entropy_midchain(state);
    state.measurements.entanglement_entropies        = measure::entanglement_entropies(state);
    if(state.get_algorithm() == AlgorithmType::fLBIT) {
        state.measurements.number_entropy_current  = measure::number_entropy_current(state);
        state.measurements.number_entropy_midchain = measure::number_entropy_midchain(state);
        state.measurements.number_entropies        = measure::number_entropies(state);
    }

    state.measurements.renyi_2         = measure::renyi_entropies(state, 2);
    state.measurements.renyi_3         = measure::renyi_entropies(state, 3);
    state.measurements.renyi_4         = measure::renyi_entropies(state, 4);
    state.measurements.renyi_inf       = measure::renyi_entropies(state, std::numeric_limits<double>::infinity());
    state.measurements.spin_components = measure::spin_components(state);

    if(state.get_algorithm() == AlgorithmType::xDMRG) correlation_matrix_xyz(state);
    if(state.get_algorithm() == AlgorithmType::xDMRG) expectation_values_xyz(state);
}

size_t tools::finite::measure::length(const TensorsFinite &tensors) { return tensors.get_length(); }
size_t tools::finite::measure::length(const StateFinite &state) { return state.get_length(); }

double tools::finite::measure::norm(const StateFinite &state, bool full) {
    if(state.measurements.norm) return state.measurements.norm.value();
    double norm;
    auto   t_chi = tid::tic_scope("norm");
    if(not full and state.is_normalized_on_all_sites()) {
        // We know the all sites are normalized. We can check that the current position is normalized
        const auto  pos = std::clamp(state.get_position<long>(), 0l, state.get_length<long>());
        const auto &mps = state.get_mps_site(pos);
        tools::log->trace("Measuring norm using site {} with dimensions {}", pos, mps.dimensions());
        norm = tools::common::contraction::contract_mps_norm(mps.get_M());
    } else if(not full and state.is_normalized_on_non_active_sites() and not state.active_sites.empty()) {
        tools::log->trace("Measuring norm using active sites {}", state.active_sites);
        Eigen::Tensor<cplx, 2> chain;
        Eigen::Tensor<cplx, 2> temp;
        bool                   first = true;
        for(const auto &pos : state.active_sites) {
            const Eigen::Tensor<cplx, 3> &M = state.get_mps_site(pos).get_M();
            if(first) {
                chain = tools::common::contraction::contract_mps_partial(M, {0, 1});
                first = false;
                continue;
            }
            temp.resize(tenx::array2{M.dimension(2), M.dimension(2)});
            temp.device(tenx::omp::getDevice()) = chain.contract(M, tenx::idx({0}, {1})).contract(M.conjugate(), tenx::idx({0, 1}, {1, 0}));

            chain = temp;
        }
        norm = std::abs(tenx::MatrixMap(chain).trace());
    } else {
        tools::log->trace("Measuring norm on full chain");
        Eigen::Tensor<cplx, 2> chain;
        Eigen::Tensor<cplx, 2> temp;
        bool                   first = true;
        for(const auto &mps : state.mps_sites) {
            const auto &M = mps->get_M();
            if(first) {
                chain = tools::common::contraction::contract_mps_partial(M, {0, 1});
                first = false;
                continue;
            }
            temp.resize(tenx::array2{M.dimension(2), M.dimension(2)});
            temp.device(tenx::omp::getDevice()) = chain.contract(M, tenx::idx({0}, {1})).contract(M.conjugate(), tenx::idx({0, 1}, {1, 0}));

            chain = temp;
        }
        norm = std::abs(tenx::MatrixMap(chain).trace());
    }

    if(std::abs(norm - 1.0) > settings::precision::max_norm_error) tools::log->debug("Norm far from unity: {:.16f}", norm);
    state.measurements.norm = norm;
    return state.measurements.norm.value();
}

long tools::finite::measure::bond_dimension_current(const StateFinite &state) {
    if(state.measurements.bond_dimension_current) return state.measurements.bond_dimension_current.value();
    if(state.has_center_point())
        state.measurements.bond_dimension_current = state.current_bond().dimension(0);
    else
        state.measurements.bond_dimension_current = 1;
    return state.measurements.bond_dimension_current.value();
}

long tools::finite::measure::bond_dimension_midchain(const StateFinite &state) {
    if(state.measurements.bond_dimension_midchain) return state.measurements.bond_dimension_midchain.value();
    state.measurements.bond_dimension_midchain = state.midchain_bond().dimension(0);
    return state.measurements.bond_dimension_midchain.value();
}

std::vector<long> tools::finite::measure::bond_dimensions(const StateFinite &state) {
    if(state.measurements.bond_dimensions) return state.measurements.bond_dimensions.value();
    auto              t_chi = tid::tic_scope("chi");
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

std::vector<long> tools::finite::measure::bond_dimensions_merged(const StateFinite &state) {
    // Here we calculate the bond dimensions of the bonds that were merged into the full state in the last step
    // For instance, if the active sites are {2,3,4,5,6} this returns the 4 bonds connecting {2,3}, {3,4}, {4,5} and {5,6}
    auto t_chi = tid::tic_scope("chi_merged");
    if(state.active_sites.empty()) return {};
    if(state.active_sites.size() == 1) {
        // Because of subspace expansion, the only bond dimension that grows is the one directly behind
        // mps, relative to the current direction.
        if(state.get_direction() == 1) return {state.get_mps_site(state.active_sites[0]).get_chiL()};
        if(state.get_direction() != 1) return {state.get_mps_site(state.active_sites[0]).get_chiR()};
    }
    if(state.active_sites.size() == 2) return {state.get_mps_site(state.active_sites[0]).get_chiR()};
    std::vector<long> bond_dimensions;
    for(const auto &pos : state.active_sites) {
        if(&pos == &state.active_sites.front()) continue;
        const auto &mps = state.get_mps_site(pos);
        bond_dimensions.push_back(mps.get_chiL());
    }
    return bond_dimensions;
}

double tools::finite::measure::entanglement_entropy_current(const StateFinite &state) {
    if(state.measurements.entanglement_entropy_current) return state.measurements.entanglement_entropy_current.value();
    auto t_ent = tid::tic_scope("neumann_entropy");
    if(state.has_center_point()) {
        auto                  &LC                       = state.current_bond();
        Eigen::Tensor<cplx, 0> SE                       = -LC.square().contract(LC.square().log().eval(), tenx::idx({0}, {0}));
        state.measurements.entanglement_entropy_current = std::abs(SE(0));
    } else
        state.measurements.entanglement_entropy_current = 0;
    return state.measurements.entanglement_entropy_current.value();
}

double tools::finite::measure::entanglement_entropy_midchain(const StateFinite &state) {
    if(state.measurements.entanglement_entropy_midchain) return state.measurements.entanglement_entropy_midchain.value();
    auto                   t_ent                     = tid::tic_scope("neumann_entropy");
    auto                  &LC                        = state.midchain_bond();
    Eigen::Tensor<cplx, 0> SE                        = -LC.square().contract(LC.square().log().eval(), tenx::idx({0}, {0}));
    state.measurements.entanglement_entropy_midchain = std::abs(SE(0));
    return state.measurements.entanglement_entropy_midchain.value();
}

std::vector<double> tools::finite::measure::entanglement_entropies(const StateFinite &state) {
    if(state.measurements.entanglement_entropies) return state.measurements.entanglement_entropies.value();
    auto                t_ent = tid::tic_scope("neumann_entropy");
    std::vector<double> entanglement_entropies;
    entanglement_entropies.reserve(state.get_length() + 1);
    if(not state.has_center_point()) entanglement_entropies.emplace_back(0);
    for(const auto &mps : state.mps_sites) {
        auto                  &L  = mps->get_L();
        Eigen::Tensor<cplx, 0> SE = -L.square().contract(L.square().log().eval(), tenx::idx({0}, {0}));
        entanglement_entropies.emplace_back(std::abs(SE(0)));
        if(mps->isCenter()) {
            auto &LC = mps->get_LC();
            SE       = -LC.square().contract(LC.square().log().eval(), tenx::idx({0}, {0}));
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

std::vector<double> tools::finite::measure::renyi_entropies(const StateFinite &state, double q) {
    auto inf = std::numeric_limits<double>::infinity();
    if(q == 1.0) return entanglement_entropies(state);
    if(q == 2.0 and state.measurements.renyi_2) return state.measurements.renyi_2.value();
    if(q == 3.0 and state.measurements.renyi_3) return state.measurements.renyi_3.value();
    if(q == 4.0 and state.measurements.renyi_4) return state.measurements.renyi_4.value();
    if(q == inf and state.measurements.renyi_inf) return state.measurements.renyi_inf.value();
    auto                t_ren = tid::tic_scope("renyi_entropy");
    std::vector<double> renyi_q;
    renyi_q.reserve(state.get_length() + 1);
    if(not state.has_center_point()) renyi_q.emplace_back(0);
    for(const auto &mps : state.mps_sites) {
        const auto            &L = mps->get_L();
        Eigen::Tensor<cplx, 0> RE;
        if(q == inf)
            RE(0) = -2.0 * std::log(L(0));
        else
            RE = 1.0 / (1.0 - q) * L.pow(2.0 * q).sum().log();
        renyi_q.emplace_back(std::abs(RE(0)));
        if(mps->isCenter()) {
            const auto &LC = mps->get_LC();
            if(q == inf)
                RE(0) = -2.0 * std::log(LC(0));
            else
                RE = 1.0 / (1.0 - q) * LC.pow(2.0 * q).sum().log();
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
    if(q == inf) {
        state.measurements.renyi_inf = renyi_q;
        return state.measurements.renyi_inf.value();
    }
    return renyi_q;
}

std::array<double, 3> tools::finite::measure::spin_components(const StateFinite &state) {
    if(state.measurements.spin_components) return state.measurements.spin_components.value();
    double spin_x                      = measure::spin_component(state, qm::spin::half::sx);
    double spin_y                      = measure::spin_component(state, qm::spin::half::sy);
    double spin_z                      = measure::spin_component(state, qm::spin::half::sz);
    state.measurements.spin_components = {spin_x, spin_y, spin_z};
    return state.measurements.spin_components.value();
}

double tools::finite::measure::spin_component(const StateFinite &state, const Eigen::Matrix2cd &paulimatrix) {
    auto t_spn       = tid::tic_scope("spin");
    auto [mpo, L, R] = qm::mpo::pauli_mpo(paulimatrix);
    Eigen::Tensor<cplx, 3> temp;
    for(const auto &mps : state.mps_sites) {
        const auto &M = mps->get_M();
        temp.resize(M.dimension(2), M.dimension(2), mpo.dimension(1));
        temp.device(tenx::omp::getDevice()) =
            L.contract(M, tenx::idx({0}, {1})).contract(M.conjugate(), tenx::idx({0}, {1})).contract(mpo, tenx::idx({0, 1, 3}, {0, 2, 3}));
        L = temp;
    }

    if(L.dimensions() != R.dimensions()) throw std::runtime_error("spin_component(): L and R dimension mismatch");
    Eigen::Tensor<cplx, 0> spin_tmp = L.contract(R, tenx::idx({0, 1, 2}, {0, 1, 2}));
    double                 spin     = std::real(spin_tmp(0));
    return spin;
}

double tools::finite::measure::spin_component(const StateFinite &state, std::string_view axis) {
    if(axis.find('x') != std::string_view::npos) return measure::spin_component(state, qm::spin::half::sx);
    if(axis.find('y') != std::string_view::npos) return measure::spin_component(state, qm::spin::half::sy);
    if(axis.find('z') != std::string_view::npos) return measure::spin_component(state, qm::spin::half::sz);
    throw std::runtime_error("Unexpected axis: " + std::string(axis));
}

std::vector<double> tools::finite::measure::truncation_errors(const StateFinite &state) {
    if(state.measurements.truncation_errors) return state.measurements.truncation_errors.value();
    auto                t_chi = tid::tic_scope("trunc");
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

std::vector<double> tools::finite::measure::truncation_errors_active(const StateFinite &state) {
    if(state.active_sites.empty()) return {};
    if(state.active_sites.size() == 1) {
        // Because of subspace expansion, the only bond dimension that grows is the one directly behind
        // mps, relative to the current direction.
        if(state.get_direction() == 1) return {state.get_mps_site(state.active_sites[0]).get_truncation_error()};
        if(state.get_direction() != 1) return {state.get_mps_site(state.active_sites[0]).get_truncation_error_LC()};
    }
    if(state.active_sites.size() == 2) return {state.get_mps_site(state.active_sites[0]).get_truncation_error_LC()};
    std::vector<double> truncation_errors;
    for(const auto &pos : state.active_sites) {
        if(&pos == &state.active_sites.front()) continue;
        const auto &mps = state.get_mps_site(pos);
        truncation_errors.push_back(mps.get_truncation_error());
    }
    return truncation_errors;
}

Eigen::Tensor<cplx, 1> tools::finite::measure::mps_wavefn(const StateFinite &state) {
    Eigen::Tensor<cplx, 2> temp;
    Eigen::Tensor<cplx, 2> chain(1, 1);
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
        temp      = chain.contract(mps->get_M(), tenx::idx({1}, {1})).reshape(tenx::array2{dimL * dim0, dimR});
        chain     = temp;
    }

    Eigen::Tensor<cplx, 1> mps_chain  = chain.reshape(tenx::array1{chain.dimension(0)});
    double                 norm_chain = tenx::VectorMap(chain).norm();
    if(std::abs(norm_chain - 1.0) > settings::precision::max_norm_error) {
        tools::log->warn("Norm far from unity: {}", norm_chain);
        throw std::runtime_error("Norm too far from unity: " + std::to_string(norm_chain));
    }
    return mps_chain;
}

template<typename state_or_mps_type>
double tools::finite::measure::energy_minus_energy_reduced(const state_or_mps_type &state, const ModelFinite &model, const EdgesFinite &edges,
                                                           MeasurementsTensorsFinite *measurements) {
    if(measurements != nullptr and measurements->energy_minus_energy_reduced) return measurements->energy_minus_energy_reduced.value();
    if constexpr(std::is_same_v<state_or_mps_type, StateFinite>) {
        if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
            throw std::runtime_error(fmt::format("Could not compute energy: active sites are not equal: state {} | model {} | edges {}", state.active_sites,
                                                 model.active_sites, edges.active_sites));
        return tools::finite::measure::energy_minus_energy_reduced(state.get_multisite_mps(), model, edges, measurements);
    } else {
        auto        t_msr = tid::tic_scope("measure");
        const auto &mpo   = model.get_multisite_mpo();
        const auto &env   = edges.get_multisite_env_ene_blk();
        if constexpr(settings::debug)
            tools::log->trace("Measuring energy: state dims {} | model sites {} dims {} | edges sites {} dims [L{} R{}]", state.dimensions(),
                              model.active_sites, mpo.dimensions(), edges.active_sites, env.L.dimensions(), env.R.dimensions());
        auto   t_ene        = tid::tic_scope("ene");
        double e_minus_ered = tools::common::contraction::expectation_value(state, mpo, env.L, env.R);
        if(measurements != nullptr) measurements->energy_minus_energy_reduced = e_minus_ered;
        return e_minus_ered;
    }
}

template double tools::finite::measure::energy_minus_energy_reduced(const StateFinite &, const ModelFinite &model, const EdgesFinite &edges,
                                                                    MeasurementsTensorsFinite *measurements);
template double tools::finite::measure::energy_minus_energy_reduced(const Eigen::Tensor<cplx, 3> &, const ModelFinite &model, const EdgesFinite &edges,
                                                                    MeasurementsTensorsFinite *measurements);

template<typename state_or_mps_type>
double tools::finite::measure::energy(const state_or_mps_type &state, const ModelFinite &model, const EdgesFinite &edges,
                                      MeasurementsTensorsFinite *measurements) {
    if(measurements != nullptr and measurements->energy) return measurements->energy.value();
    // This measures the actual energy of the system regardless of the reduced/non-reduced state of the MPO's
    // If they are reduced, then
    //      "Actual energy" = (E - E_reduced) + E_reduced = (~0) + E_reduced = E
    // Else
    //      "Actual energy" = (E - E_reduced) + E_reduced = E  + 0 = E
    double energy;
    if constexpr(std::is_same_v<state_or_mps_type, StateFinite>)
        energy = tools::finite::measure::energy_minus_energy_reduced(state.get_multisite_mps(), model, edges, measurements) + model.get_energy_reduced();
    else
        energy = tools::finite::measure::energy_minus_energy_reduced(state, model, edges, measurements) + model.get_energy_reduced();

    if(measurements != nullptr) measurements->energy = energy;
    return energy;
}

template double tools::finite::measure::energy(const StateFinite &, const ModelFinite &model, const EdgesFinite &edges,
                                               MeasurementsTensorsFinite *measurements);
template double tools::finite::measure::energy(const Eigen::Tensor<cplx, 3> &, const ModelFinite &model, const EdgesFinite &edges,
                                               MeasurementsTensorsFinite *measurements);

template<typename state_or_mps_type>
double tools::finite::measure::energy_per_site(const state_or_mps_type &state, const ModelFinite &model, const EdgesFinite &edges,
                                               MeasurementsTensorsFinite *measurements) {
    double energy_per_site = tools::finite::measure::energy(state, model, edges, measurements) / static_cast<double>(model.get_length());
    if(measurements != nullptr) measurements->energy_per_site = energy_per_site;
    return energy_per_site;
}

template double tools::finite::measure::energy_per_site(const StateFinite &, const ModelFinite &model, const EdgesFinite &edges,
                                                        MeasurementsTensorsFinite *measurements);
template double tools::finite::measure::energy_per_site(const Eigen::Tensor<cplx, 3> &, const ModelFinite &model, const EdgesFinite &edges,
                                                        MeasurementsTensorsFinite *measurements);

template<typename state_or_mps_type>
double tools::finite::measure::energy_variance(const state_or_mps_type &state, const ModelFinite &model, const EdgesFinite &edges,
                                               MeasurementsTensorsFinite *measurements) {
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

    if constexpr(std::is_same_v<state_or_mps_type, StateFinite>) {
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

        auto   t_msr = tid::tic_scope("measure");
        double E2    = energy * energy;

        if(not num::all_equal(model.active_sites, edges.active_sites))
            throw std::runtime_error(
                fmt::format("Could not compute energy variance: active sites are not equal: model {} | edges {}", model.active_sites, edges.active_sites));

        const auto &mpo = model.get_multisite_mpo_squared();
        const auto &env = edges.get_multisite_env_var_blk();
        if constexpr(settings::debug)
            tools::log->trace("Measuring energy variance: state dims {} | model sites {} dims {} | edges sites {} dims [L{} R{}]", state.dimensions(),
                              model.active_sites, mpo.dimensions(), edges.active_sites, env.L.dimensions(), env.R.dimensions());

        if(state.dimension(0) != mpo.dimension(2))
            throw std::runtime_error(
                fmt::format("State and model have incompatible physical dimension: state dim {} | model dim {}", state.dimension(0), mpo.dimension(2)));
        auto   t_var = tid::tic_scope("var");
        double H2    = tools::common::contraction::expectation_value(state, mpo, env.L, env.R);
        double var   = std::abs(H2 - E2);
        if(measurements != nullptr) measurements->energy_variance = var;
        return var;
    }
}

template double tools::finite::measure::energy_variance(const StateFinite &, const ModelFinite &model, const EdgesFinite &edges,
                                                        MeasurementsTensorsFinite *measurements);
template double tools::finite::measure::energy_variance(const Eigen::Tensor<cplx, 3> &, const ModelFinite &model, const EdgesFinite &edges,
                                                        MeasurementsTensorsFinite *measurements);

template<typename state_or_mps_type>
double tools::finite::measure::energy_variance_per_site(const state_or_mps_type &state, const ModelFinite &model, const EdgesFinite &edges,
                                                        MeasurementsTensorsFinite *measurements) {
    double energy_variance_per_site = tools::finite::measure::energy_variance(state, model, edges, measurements) / static_cast<double>(model.get_length());
    if(measurements != nullptr) measurements->energy_variance_per_site = energy_variance_per_site;
    return energy_variance_per_site;
}

template double tools::finite::measure::energy_variance_per_site(const StateFinite &, const ModelFinite &model, const EdgesFinite &edges,
                                                                 MeasurementsTensorsFinite *measurements);
template double tools::finite::measure::energy_variance_per_site(const Eigen::Tensor<cplx, 3> &, const ModelFinite &model, const EdgesFinite &edges,
                                                                 MeasurementsTensorsFinite *measurements);

template<typename state_or_mps_type>
double tools::finite::measure::energy_normalized(const state_or_mps_type &state, const ModelFinite &model, const EdgesFinite &edges, double energy_min_per_site,
                                                 double energy_max_per_site, MeasurementsTensorsFinite *measurements) {
    return (tools::finite::measure::energy_per_site(state, model, edges, measurements) - energy_min_per_site) / (energy_max_per_site - energy_min_per_site);
}

template double tools::finite::measure::energy_normalized(const StateFinite &, const ModelFinite &model, const EdgesFinite &edges, double, double,
                                                          MeasurementsTensorsFinite *measurements);
template double tools::finite::measure::energy_normalized(const Eigen::Tensor<cplx, 3> &, const ModelFinite &model, const EdgesFinite &edges, double, double,
                                                          MeasurementsTensorsFinite *measurements);

extern double tools::finite::measure::energy_reduced(const TensorsFinite &tensors) { return tensors.model->get_energy_reduced(); }
extern double tools::finite::measure::energy_per_site_reduced(const TensorsFinite &tensors) { return tensors.model->get_energy_per_site_reduced(); }

double tools::finite::measure::energy_minus_energy_reduced(const TensorsFinite &tensors) {
    tensors.assert_edges_ene();
    return energy_minus_energy_reduced(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
}

double tools::finite::measure::energy(const TensorsFinite &tensors) {
    if(tensors.measurements.energy)
        return tensors.measurements.energy.value();
    else {
        tensors.assert_edges_ene();
        tensors.measurements.energy = tools::finite::measure::energy(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
        return tensors.measurements.energy.value();
    }
}

double tools::finite::measure::energy_per_site(const TensorsFinite &tensors) {
    if(tensors.measurements.energy_per_site)
        return tensors.measurements.energy_per_site.value();
    else {
        tensors.assert_edges_ene();
        tensors.measurements.energy_per_site = tools::finite::measure::energy_per_site(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
        return tensors.measurements.energy_per_site.value();
    }
}

double tools::finite::measure::energy_variance(const TensorsFinite &tensors) {
    if(tensors.measurements.energy_variance)
        return tensors.measurements.energy_variance.value();
    else {
        tensors.assert_edges_var();
        tensors.measurements.energy_variance = tools::finite::measure::energy_variance(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
        return tensors.measurements.energy_variance.value();
    }
}

double tools::finite::measure::energy_variance_per_site(const TensorsFinite &tensors) {
    if(tensors.measurements.energy_variance_per_site)
        return tensors.measurements.energy_variance_per_site.value();
    else {
        tensors.assert_edges_var();
        tensors.measurements.energy_variance_per_site =
            tools::finite::measure::energy_variance_per_site(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
        return tensors.measurements.energy_variance_per_site.value();
    }
}

double tools::finite::measure::energy_normalized(const TensorsFinite &tensors, double emin, double emax) {
    tensors.assert_edges_ene();
    return tools::finite::measure::energy_normalized(*tensors.state, *tensors.model, *tensors.edges, emin, emax);
}

double tools::finite::measure::energy_minus_energy_reduced(const StateFinite &state, const TensorsFinite &tensors) {
    return tools::finite::measure::energy_minus_energy_reduced(state, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy(const StateFinite &state, const TensorsFinite &tensors) {
    return tools::finite::measure::energy(state, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_per_site(const StateFinite &state, const TensorsFinite &tensors) {
    return tools::finite::measure::energy_per_site(state, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_variance(const StateFinite &state, const TensorsFinite &tensors) {
    return tools::finite::measure::energy_variance(state, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_variance_per_site(const StateFinite &state, const TensorsFinite &tensors) {
    return tools::finite::measure::energy_variance_per_site(state, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_normalized(const StateFinite &state, const TensorsFinite &tensors, double emin, double emax) {
    return tools::finite::measure::energy_normalized(state, *tensors.model, *tensors.edges, emin, emax);
}

double tools::finite::measure::energy_minus_energy_reduced(const Eigen::Tensor<cplx, 3> &mps, const TensorsFinite &tensors) {
    return tools::finite::measure::energy_minus_energy_reduced(mps, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy(const Eigen::Tensor<cplx, 3> &mps, const TensorsFinite &tensors) {
    return tools::finite::measure::energy(mps, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_per_site(const Eigen::Tensor<cplx, 3> &mps, const TensorsFinite &tensors) {
    return tools::finite::measure::energy_per_site(mps, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_variance(const Eigen::Tensor<cplx, 3> &mps, const TensorsFinite &tensors) {
    return tools::finite::measure::energy_variance(mps, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_variance_per_site(const Eigen::Tensor<cplx, 3> &mps, const TensorsFinite &tensors) {
    return tools::finite::measure::energy_variance_per_site(mps, *tensors.model, *tensors.edges);
}
double tools::finite::measure::energy_normalized(const Eigen::Tensor<cplx, 3> &mps, const TensorsFinite &tensors, double emin, double emax) {
    return tools::finite::measure::energy_normalized(mps, *tensors.model, *tensors.edges, emin, emax);
}

double tools::finite::measure::max_gradient(const Eigen::Tensor<cplx, 3> &mps, const TensorsFinite &tensors) {
    auto v = Eigen::Map<const Eigen::VectorXcd>(mps.data(), mps.size());

    auto en1 = tensors.get_multisite_env_ene_blk();
    auto H1t = tools::common::contraction::matrix_vector_product(mps, tensors.get_multisite_mpo(), en1.L, en1.R);
    auto Hv  = Eigen::Map<const Eigen::VectorXcd>(H1t.data(), H1t.size());
    auto vHv = v.dot(Hv);

    auto en2  = tensors.get_multisite_env_var_blk();
    auto H2t  = tools::common::contraction::matrix_vector_product(mps, tensors.get_multisite_mpo_squared(), en2.L, en2.R);
    auto H2v  = Eigen::Map<const Eigen::VectorXcd>(H2t.data(), H2t.size());
    auto vH2v = v.dot(H2v);

    double var    = std::abs(vH2v - vHv * vHv);
    auto   var_1  = 1.0 / var / std::log(10);
    auto   norm_1 = 1.0 / v.norm();
    auto   pref   = tensors.is_real() ? 1.0 : 2.0;
    auto   grad   = pref * var_1 * norm_1 * (H2v - 2.0 * vHv * Hv - (vH2v - 2.0 * vHv * vHv) * v); // Factor 2 for complex
    return grad.template lpNorm<Eigen::Infinity>();
}

template<typename Scalar>
double tools::finite::measure::max_gradient(const Eigen::Tensor<Scalar, 3> &mps, const Eigen::Tensor<Scalar, 4> &mpo1, const Eigen::Tensor<Scalar, 3> &en1L,
                                            const Eigen::Tensor<Scalar, 3> &en1R, const Eigen::Tensor<Scalar, 4> &mpo2, const Eigen::Tensor<Scalar, 3> &en2L,
                                            const Eigen::Tensor<Scalar, 3> &en2R) {
    using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    auto v           = Eigen::Map<const VectorType>(mps.data(), mps.size());
    auto H1t         = tools::common::contraction::matrix_vector_product(mps, mpo1, en1L, en1R);
    auto Hv          = Eigen::Map<const VectorType>(H1t.data(), H1t.size());
    auto vHv         = v.dot(Hv);

    auto H2t  = tools::common::contraction::matrix_vector_product(mps, mpo2, en2L, en2R);
    auto H2v  = Eigen::Map<const VectorType>(H2t.data(), H2t.size());
    auto vH2v = v.dot(H2v);

    auto var    = std::abs(vH2v - vHv * vHv);
    auto var_1  = 1.0 / var / std::log(10);
    auto norm_1 = 1.0 / v.norm();
    auto pref   = std::is_same_v<Scalar, cplx> ? 2.0 : 1.0;
    auto grad   = pref * var_1 * norm_1 * (H2v - 2.0 * vHv * Hv - (vH2v - 2.0 * vHv * vHv) * v); // Factor 2 for complex
    return grad.template lpNorm<Eigen::Infinity>();
}

template double tools::finite::measure::max_gradient(const Eigen::Tensor<real, 3> &mps, const Eigen::Tensor<real, 4> &mpo1, const Eigen::Tensor<real, 3> &en1L,
                                                     const Eigen::Tensor<real, 3> &en1R, const Eigen::Tensor<real, 4> &mpo2, const Eigen::Tensor<real, 3> &en2L,
                                                     const Eigen::Tensor<real, 3> &en2R);
template double tools::finite::measure::max_gradient(const Eigen::Tensor<cplx, 3> &mps, const Eigen::Tensor<cplx, 4> &mpo1, const Eigen::Tensor<cplx, 3> &en1L,
                                                     const Eigen::Tensor<cplx, 3> &en1R, const Eigen::Tensor<cplx, 4> &mpo2, const Eigen::Tensor<cplx, 3> &en2L,
                                                     const Eigen::Tensor<cplx, 3> &en2R);

template<typename Scalar>
double tools::finite::measure::max_gradient(const Eigen::Tensor<Scalar, 3> &mps, const Eigen::Tensor<Scalar, 4> &mpo2, const Eigen::Tensor<Scalar, 3> &en2L,
                                            const Eigen::Tensor<Scalar, 3> &en2R) {
    using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    auto v           = Eigen::Map<const VectorType>(mps.data(), mps.size());
    auto H2t         = tools::common::contraction::matrix_vector_product(mps, mpo2, en2L, en2R);
    auto H2v         = Eigen::Map<const VectorType>(H2t.data(), H2t.size());
    auto vH2v        = v.dot(H2v);

    double var    = std::abs(vH2v);
    auto   var_1  = 1.0 / var / std::log(10);
    auto   norm_1 = 1.0 / v.norm();
    auto   pref   = std::is_same_v<Scalar, double> ? 2.0 : 1.0; // Factor 2 for double, from the derivative of v^2
    auto   grad   = pref * var_1 * norm_1 * (H2v - vH2v * v);
    return grad.template lpNorm<Eigen::Infinity>();
}
template double tools::finite::measure::max_gradient(const Eigen::Tensor<real, 3> &mps, const Eigen::Tensor<real, 4> &mpo2, const Eigen::Tensor<real, 3> &en2L,
                                                     const Eigen::Tensor<real, 3> &en2R);
template double tools::finite::measure::max_gradient(const Eigen::Tensor<cplx, 3> &mps, const Eigen::Tensor<cplx, 4> &mpo2, const Eigen::Tensor<cplx, 3> &en2L,
                                                     const Eigen::Tensor<cplx, 3> &en2R);

double tools::finite::measure::expectation_value(const StateFinite &state, const std::vector<LocalObservableOp> &ops) {
    if(state.mps_sites.empty()) throw std::runtime_error("expectation_value: state.mps_sites is empty");
    if(ops.empty()) throw std::runtime_error("expectation_value: obs is empty");
    auto                   d0    = state.mps_sites.front()->get_chiL();
    Eigen::Tensor<cplx, 2> chain = tenx::TensorCast(Eigen::MatrixXcd::Identity(d0, d0));
    Eigen::Tensor<cplx, 2> temp;
    for(const auto &mps : state.mps_sites) {
        const auto &M     = mps->get_M();
        const auto  pos   = mps->get_position<long>();
        const auto  ob_it = std::find_if(ops.begin(), ops.end(), [&pos](const auto &ob) { return ob.pos == pos and not ob.used; });
        temp.resize(tenx::array2{M.dimension(2), M.dimension(2)});
        if(ob_it != ops.end()) {
            auto M_op                           = ob_it->op.contract(M, tenx::idx({1}, {0}));
            temp.device(tenx::omp::getDevice()) = chain.contract(M_op, tenx::idx({0}, {1})).contract(M.conjugate(), tenx::idx({0, 1}, {1, 0}));
            ob_it->used                         = true;
        } else {
            temp.device(tenx::omp::getDevice()) = chain.contract(M, tenx::idx({0}, {1})).contract(M.conjugate(), tenx::idx({0, 1}, {1, 0}));
        }
        chain = temp;
    }
    auto expval = tenx::MatrixMap(chain).trace();
    if(std::imag(expval) > 1e-14) tools::log->warn("expectation_value: result has imaginary part: {:+8.2e}{:+8.2e} i", std::real(expval), std::imag(expval));
    if(std::imag(expval) > 1e-8)
        throw except::runtime_error("expectation_value: result has imaginary part: {:+8.2e}{:+8.2e} i", std::real(expval), std::imag(expval));
    return std::real(expval);
}

double tools::finite::measure::expectation_value(const StateFinite &state, const std::vector<LocalObservableMpo> &mpos) {
    if(state.mps_sites.empty()) throw std::runtime_error("expectation_value: state.mps_sites is empty");
    if(mpos.empty()) throw std::runtime_error("expectation_value: obs is empty");

    // Generate a string of mpos for each site. If a site has no local observable given in obs, insert an identity MPO there.
    auto mpodims = mpos.front().mpo.dimensions();

    for(auto &ob : mpos) {
        if(ob.mpo.dimension(2) != ob.mpo.dimension(3)) throw except::runtime_error("expectation_value: given mpo's of unequal spin dimension up and down");
        if(ob.mpo.dimension(0) != ob.mpo.dimension(1)) throw except::runtime_error("expectation_value: given mpo's of unequal bond dimension left and right");
        if(ob.mpo.dimensions() != mpodims) throw except::runtime_error("expectation_value: given mpo's of unequal dimensions");
    }

    // Create compatible edges
    Eigen::Tensor<cplx, 1> Ledge(mpodims[0]); // The left  edge
    Eigen::Tensor<cplx, 1> Redge(mpodims[1]); // The right edge
    Ledge(mpodims[0] - 1) = 1;
    Redge(0)              = 1;
    Eigen::Tensor<cplx, 3> Ledge3, Redge3;
    {
        auto mpsdims = state.mps_sites.front()->dimensions();
        long mpsDim  = mpsdims[1];
        long mpoDim  = mpodims[0];
        Ledge3.resize(tenx::array3{mpsDim, mpsDim, mpoDim});
        Ledge3.setZero();
        for(long i = 0; i < mpsDim; i++) {
            std::array<long, 1> extent1                     = {mpoDim};
            std::array<long, 3> offset3                     = {i, i, 0};
            std::array<long, 3> extent3                     = {1, 1, mpoDim};
            Ledge3.slice(offset3, extent3).reshape(extent1) = Ledge;
        }
    }
    {
        auto mpsdims = state.mps_sites.back()->dimensions();
        long mpsDim  = mpsdims[2];
        long mpoDim  = mpodims[1];
        Redge3.resize(tenx::array3{mpsDim, mpsDim, mpoDim});
        Redge3.setZero();
        for(long i = 0; i < mpsDim; i++) {
            std::array<long, 1> extent1                     = {mpoDim};
            std::array<long, 3> offset3                     = {i, i, 0};
            std::array<long, 3> extent3                     = {1, 1, mpoDim};
            Redge3.slice(offset3, extent3).reshape(extent1) = Redge;
        }
    }

    // Generate an identity mpo with the same dimensions as the ones in obs
    Eigen::Tensor<cplx, 1> ones(mpodims[0] * mpodims[2]);
    ones.setConstant(1.0);
    Eigen::Tensor<cplx, 4> mpoI = tenx::asDiagonal(ones).reshape(mpodims);

    // Start applying the mpo or identity on each site starting from Ledge3
    Eigen::Tensor<cplx, 3> temp;
    for(const auto &mps : state.mps_sites) {
        const auto  pos   = mps->get_position<long>();
        const auto &M     = mps->get_M();
        const auto  ob_it = std::find_if(mpos.begin(), mpos.end(), [&pos](const auto &ob) { return ob.pos == pos and not ob.used; });
        const auto &mpo   = ob_it != mpos.end() ? ob_it->mpo : mpoI; // Choose the operator or an identity
        if(ob_it != mpos.end()) ob_it->used = true;
        temp.resize(M.dimension(2), M.dimension(2), mpo.dimension(1));
        temp   = M.contract(Ledge3, tenx::idx({0}, {1})).contract(mpo, tenx::idx({0}, {1})).contract(M.conjugate(), tenx::idx({0, 1, 3}, {0, 2, 3}));
        Ledge3 = temp;
    }

    if(Ledge3.dimensions() != Redge3.dimensions())
        throw except::runtime_error("expectation_value: Ledge3 and Redge3 dimension mismatch: {} != {}", Ledge3.dimensions(), Redge3.dimensions());

    // Finish by contracting Redge3
    Eigen::Tensor<cplx, 0> expval = Ledge3.contract(Redge3, tenx::idx({0, 1, 2}, {0, 1, 2}));

    if(std::imag(expval(0)) > 1e-14)
        tools::log->warn("expectation_value: result has imaginary part: {:+8.2e}{:+8.2e} i", std::real(expval(0)), std::imag(expval(0)));
    if(std::imag(expval(0)) > 1e-8)
        throw except::runtime_error("expectation_value: result has imaginary part: {:+8.2e}{:+8.2e} i", std::real(expval(0)), std::imag(expval(0)));
    return std::real(expval(0));
}

Eigen::Tensor<double, 1> tools::finite::measure::expectation_values(const StateFinite &state, const Eigen::Tensor<cplx, 2> &op) {
    tools::log->trace("Measuring local expectation values");
    long                     len = state.get_length<long>();
    Eigen::Tensor<double, 1> expvals(len);
    expvals.setZero();
    for(long pos = 0; pos < len; pos++) {
        LocalObservableOp ob1 = {op, pos};
        expvals(pos)          = expectation_value(state, {ob1});
    }
    return expvals;
}

double tools::finite::measure::correlation(const StateFinite &state, const Eigen::Tensor<cplx, 2> &op1, const Eigen::Tensor<cplx, 2> &op2, long pos1,
                                           long pos2) {
    if(pos1 == pos2) {
        // Stack the operators
        Eigen::Tensor<cplx, 2> op12 = op1.contract(op2, tenx::idx({0}, {1}));
        LocalObservableOp      ob12 = {op12, pos1};
        return expectation_value(state, {ob12});
    } else {
        // No need to stack
        LocalObservableOp ob1 = {op1, pos1};
        LocalObservableOp ob2 = {op2, pos2};
        return expectation_value(state, {ob1, ob2});
    }
}

Eigen::Tensor<double, 2> tools::finite::measure::correlation_matrix(const StateFinite &state, const Eigen::Tensor<cplx, 2> &op1,
                                                                    const Eigen::Tensor<cplx, 2> &op2) {
    tools::log->trace("Measuring correlation matrix");

    long                     len = state.get_length<long>();
    bool                     eq  = tenx::MatrixMap(op1) == tenx::MatrixMap(op2);
    Eigen::Tensor<double, 2> C(len, len);
    C.setZero();

    for(long pos_j = 0; pos_j < len; pos_j++) {
        for(long pos_i = pos_j; pos_i < len; pos_i++) {
            C(pos_i, pos_j) = correlation(state, op1, op2, pos_i, pos_j);
            if(eq)
                C(pos_j, pos_i) = C(pos_i, pos_j);
            else
                C(pos_j, pos_i) = correlation(state, op1, op2, pos_j, pos_i);
        }
    }
    return C;
}

void tools::finite::measure::expectation_values_xyz(const StateFinite &state) {
    Eigen::Tensor<cplx, 2> sx = tenx::TensorMap(qm::spin::half::sx);
    Eigen::Tensor<cplx, 2> sy = tenx::TensorMap(qm::spin::half::sy);
    Eigen::Tensor<cplx, 2> sz = tenx::TensorMap(qm::spin::half::sz);
    if(not state.measurements.expectation_values_sx) state.measurements.expectation_values_sx = measure::expectation_values(state, sx);
    if(not state.measurements.expectation_values_sy) state.measurements.expectation_values_sy = measure::expectation_values(state, sy);
    if(not state.measurements.expectation_values_sz) state.measurements.expectation_values_sz = measure::expectation_values(state, sz);
}

void tools::finite::measure::correlation_matrix_xyz(const StateFinite &state) {
    Eigen::Tensor<cplx, 2> sx = tenx::TensorMap(qm::spin::half::sx);
    Eigen::Tensor<cplx, 2> sy = tenx::TensorMap(qm::spin::half::sy);
    Eigen::Tensor<cplx, 2> sz = tenx::TensorMap(qm::spin::half::sz);
    if(not state.measurements.correlation_matrix_sx) state.measurements.correlation_matrix_sx = measure::correlation_matrix(state, sx, sx);
    if(not state.measurements.correlation_matrix_sy) state.measurements.correlation_matrix_sy = measure::correlation_matrix(state, sy, sy);
    if(not state.measurements.correlation_matrix_sz) state.measurements.correlation_matrix_sz = measure::correlation_matrix(state, sz, sz);
}
