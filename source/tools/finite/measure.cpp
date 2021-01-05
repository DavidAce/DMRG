//
// Created by david on 2019-02-01.
//

//
// Created by david on 2017-11-12.
//

#include "measure.h"
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>
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
        const auto  pos  = state.get_position<size_t>();
        if (pos == state.get_length()-1){
            // We know the current position is an "A" matrix, so the norm must be what's left in LC.
            Eigen::Tensor<Scalar,0>  LCsum = state.get_mps_site(pos).get_LC().square().sum(Textra::array1{0});
            return std::abs(LCsum(0));
        }
        const auto &mpsL = state.get_mps_site(pos);
        const auto &mpsR = state.get_mps_site(pos + 1);
        tools::log->trace("Measuring norm using center sites [{},{}] with dimensions [{}, {}]", pos, pos + 1, mpsL.dimensions(), mpsR.dimensions());
        Eigen::Tensor<Scalar, 0> norm_contraction;

        norm_contraction.device(Textra::omp::getDevice()) = mpsL.get_M()
                                                                .contract(mpsR.get_M(), Textra::idx({2}, {1}))
                                                                .contract(mpsL.get_M().conjugate(), Textra::idx({0, 1}, {0, 1}))
                                                                .contract(mpsR.get_M().conjugate(), Textra::idx({2, 1, 0}, {1, 2, 0}));
        norm = std::abs(norm_contraction(0));

    } else if(state.is_normalized_on_non_active_sites() and not state.active_sites.empty()) {
        tools::log->trace("Measuring norm using active sites {}", state.active_sites);
        Eigen::Tensor<Scalar, 2> chain;
        Eigen::Tensor<Scalar, 2> temp;
        bool                     first = true;
        for(auto &&pos : state.active_sites) {
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
        norm = std::abs(Textra::TensorMatrixMap(chain).trace());
    } else {
        tools::log->trace("Measuring norm on full chain");
        Eigen::Tensor<Scalar, 2> chain;
        Eigen::Tensor<Scalar, 2> temp;
        bool                     first = true;
        for(size_t pos = 0; pos < state.get_length(); pos++) {
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
        norm = std::abs(Textra::TensorMatrixMap(chain).trace());
    }

    if(std::abs(norm - 1.0) > settings::precision::max_norm_error) tools::log->debug("Norm far from unity: {:.16f}", norm);
    state.measurements.norm = norm;
    return state.measurements.norm.value();
}

long tools::finite::measure::bond_dimension_current(const class_state_finite &state) {
    if(state.measurements.bond_dimension_current) { return state.measurements.bond_dimension_current.value(); }
    state.measurements.bond_dimension_current = state.current_bond().dimension(0);
    return state.measurements.bond_dimension_current.value();
}

long tools::finite::measure::bond_dimension_midchain(const class_state_finite &state) {
    if(state.measurements.bond_dimension_midchain) { return state.measurements.bond_dimension_midchain.value(); }
    state.measurements.bond_dimension_midchain = state.midchain_bond().dimension(0);
    return state.measurements.bond_dimension_midchain.value();
}

std::vector<long> tools::finite::measure::bond_dimensions(const class_state_finite &state) {
    if(state.measurements.bond_dimensions) { return state.measurements.bond_dimensions.value(); }
    state.measurements.bond_dimensions = std::vector<long>{};
    for(size_t pos = 0; pos < state.get_length(); pos++) {
        state.measurements.bond_dimensions.value().emplace_back(state.get_mps_site(pos).get_L().dimension(0));
        if(state.get_mps_site(pos).isCenter()) { state.measurements.bond_dimensions.value().emplace_back(state.get_mps_site(pos).get_LC().dimension(0)); }
    }
    return state.measurements.bond_dimensions.value();
}

std::vector<long> tools::finite::measure::bond_dimensions_merged(const class_state_finite &state) {
    std::vector<long> bond_dimensions;
    for(auto &&pos : state.active_sites) {
        bond_dimensions.emplace_back(state.get_mps_site(pos).get_L().dimension(0));
        if(state.get_mps_site(pos).isCenter()) { bond_dimensions.emplace_back(state.get_mps_site(pos).get_LC().dimension(0)); }
    }
    if(bond_dimensions.size() > 2) {
        bond_dimensions.pop_back();
        bond_dimensions.erase(bond_dimensions.begin());
    }
    return bond_dimensions;
}

double tools::finite::measure::entanglement_entropy_current(const class_state_finite &state) {
    if(state.measurements.entanglement_entropy_current) { return state.measurements.entanglement_entropy_current.value(); }
    tools::common::profile::get_default_prof()["t_ent"]->tic();
    auto &                   LC                     = state.current_bond();
    Eigen::Tensor<Scalar, 0> SE                     = -LC.square().contract(LC.square().log().eval(), idx({0}, {0}));
    state.measurements.entanglement_entropy_current = std::real(SE(0));
    tools::common::profile::get_default_prof()["t_ent"]->toc();
    return state.measurements.entanglement_entropy_current.value();
}

double tools::finite::measure::entanglement_entropy_midchain(const class_state_finite &state) {
    if(state.measurements.entanglement_entropy_midchain) { return state.measurements.entanglement_entropy_midchain.value(); }
    tools::common::profile::get_default_prof()["t_ent"]->tic();
    auto &                   LC                      = state.midchain_bond();
    Eigen::Tensor<Scalar, 0> SE                      = -LC.square().contract(LC.square().log().eval(), idx({0}, {0}));
    state.measurements.entanglement_entropy_midchain = std::real(SE(0));
    tools::common::profile::get_default_prof()["t_ent"]->toc();
    return state.measurements.entanglement_entropy_midchain.value();
}

std::vector<double> tools::finite::measure::entanglement_entropies(const class_state_finite &state) {
    if(state.measurements.entanglement_entropies) { return state.measurements.entanglement_entropies.value(); }
    tools::common::profile::get_default_prof()["t_ent"]->tic();
    std::vector<double> entanglement_entropies;
    for(size_t pos = 0; pos < state.get_length(); pos++) {
        auto &                   L  = state.get_mps_site(pos).get_L();
        Eigen::Tensor<Scalar, 0> SE = -L.square().contract(L.square().log().eval(), idx({0}, {0}));
        entanglement_entropies.emplace_back(std::real(SE(0)));
        if(state.get_mps_site(pos).isCenter()) {
            auto &LC = state.get_mps_site(pos).get_LC();
            SE       = -LC.square().contract(LC.square().log().eval(), idx({0}, {0}));
            entanglement_entropies.emplace_back(std::real(SE(0)));
            state.measurements.entanglement_entropy_current = std::real(SE(0));
        }
    }
    state.measurements.entanglement_entropies = entanglement_entropies;
    tools::common::profile::get_default_prof()["t_ent"]->toc();
    return state.measurements.entanglement_entropies.value();
}

std::vector<double> tools::finite::measure::renyi_entropies(const class_state_finite &state, double q) {
    if(q == 1.0) return entanglement_entropies(state);
    if(q == 2.0 and state.measurements.renyi_2) return state.measurements.renyi_2.value();
    if(q == 3.0 and state.measurements.renyi_3) return state.measurements.renyi_3.value();
    if(q == 4.0 and state.measurements.renyi_4) return state.measurements.renyi_4.value();
    if(q == 100.0 and state.measurements.renyi_100) return state.measurements.renyi_100.value();
    tools::common::profile::get_default_prof()["t_ent"]->tic();
    std::vector<double> renyi_q;
    for(size_t pos = 0; pos < state.get_length(); pos++) {
        const auto &             L = state.get_mps_site(pos).get_L();
        Eigen::Tensor<Scalar, 0> RE;
        RE = (1.0 / 1.0 - q) * L.pow(2.0 * q).sum().log();
        renyi_q.emplace_back(std::real(RE(0)));
        if(state.get_mps_site(pos).isCenter()) {
            const auto &LC = state.get_mps_site(pos).get_LC();
            RE             = (1.0 / 1.0 - q) * LC.pow(2.0 * q).sum().log();
            renyi_q.emplace_back(std::real(RE(0)));
        }
    }
    tools::common::profile::get_default_prof()["t_ent"]->toc();
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
    if(state.measurements.spin_components) { return state.measurements.spin_components.value(); }
    double spin_x                      = measure::spin_component(state, qm::spinHalf::sx);
    double spin_y                      = measure::spin_component(state, qm::spinHalf::sy);
    double spin_z                      = measure::spin_component(state, qm::spinHalf::sz);
    state.measurements.spin_components = {spin_x, spin_y, spin_z};
    return state.measurements.spin_components.value();
}

double tools::finite::measure::spin_component(const class_state_finite &state, const Eigen::Matrix2cd &paulimatrix) {
    tools::common::profile::get_default_prof()["t_spn"]->tic();
    auto [mpo, L, R] = qm::mpo::pauli_mpo(paulimatrix);
    Eigen::Tensor<Scalar, 3> temp;
    for(auto && mps : state.mps_sites) {
        const auto &M = mps->get_M();
        temp.resize(M.dimension(2), M.dimension(2), mpo.dimension(1));
        temp.device(Textra::omp::getDevice()) = L.contract(M, idx({0}, {1})).contract(M.conjugate(), idx({0}, {1})).contract(mpo, idx({0, 1, 3}, {0, 2, 3}));
        L                                     = temp;
    }

    assert(L.dimensions() == R.dimensions());
    Eigen::Tensor<Scalar, 0> spin_tmp = L.contract(R, idx({0, 1, 2}, {0, 1, 2}));
    double                   spin     = std::real(spin_tmp(0));
    tools::common::profile::get_default_prof()["t_spn"]->toc();
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
    for(const auto &mps : state.mps_sites) {
        truncation_errors.emplace_back(mps->get_truncation_error());
        if(mps->isCenter()) truncation_errors.emplace_back(mps->get_truncation_error_LC());
    }
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
    double                   norm_chain = Textra::TensorVectorMap(chain).norm();
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
        tools::log->trace("Measuring energy on sites: model {} | edges {}", model.active_sites, edges.active_sites);
        tools::common::profile::get_default_prof()["t_ene"]->tic();
        double e_minus_ered = tools::common::contraction::expectation_value(state, mpo, env.L, env.R);
        tools::common::profile::get_default_prof()["t_ene"]->toc();
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
    //      "Actual energy" = (E - E_reduced) + E_reduced = (E)  + 0 = E
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
    // Depending on whether the mpo's are reduced or not we get different formulas.
    // If mpo's are reduced:
    //      Var H = <(H-E_red)^2> - <(H-E_red)>^2 = <H^2> - 2<H>E_red + E_red^2 - (<H> - E_red) ^2
    //                                            = H2    - 2*E*E_red + E_red^2 - E^2 + 2*E*E_red - E_red^2
    //                                            = H2    - E^2
    //      so Var H = <(H-E_red)^2> - energy_minus_energy_reduced^2 = H2 - ~0
    //      where H2 is computed with reduced mpo's. Note that ~0 is not exactly zero
    //      because E_red != E necessarily (though they are supposed to be very close)
    // Else, if E_red = 0 (i.e. not reduced) we get the usual formula:
    //      Var H = <(H - 0)^2> - <H - 0>^2 = H2 - E^2
    if(measurements != nullptr and measurements->energy_variance) return measurements->energy_variance.value();

    if constexpr(std::is_same_v<state_or_mps_type, class_state_finite>) {
        if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
            throw std::runtime_error(fmt::format("Could not compute energy variance: active sites are not equal: state {} | model {} | edges {}",
                                                 state.active_sites, model.active_sites, edges.active_sites));
        if(state.active_sites.empty()) throw std::runtime_error("Could not compute energy variance: active sites are empty");
        return tools::finite::measure::energy_variance(state.get_multisite_mps(), model, edges, measurements);
    } else {
        double energy = 0;
        if(model.is_reduced()) energy = tools::finite::measure::energy_minus_energy_reduced(state, model, edges, measurements);
        else
            energy = tools::finite::measure::energy(state, model, edges, measurements);
        double E2 = energy * energy;

        if(not num::all_equal(model.active_sites, edges.active_sites))
            throw std::runtime_error(
                fmt::format("Could not compute energy variance: active sites are not equal: model {} | edges {}", model.active_sites, edges.active_sites));

        const auto &mpo = model.get_multisite_mpo_squared();
        const auto &env = edges.get_multisite_var_blk();
        tools::log->trace("Measuring energy variance");
        if(state.dimension(0) != mpo.dimension(2))
            throw std::runtime_error(
                fmt::format("State and model have incompatible physical dimension: state dim {} | model dim {}", state.dimension(0), mpo.dimension(2)));
        tools::common::profile::get_default_prof()["t_var"]->tic();
        double H2 = tools::common::contraction::expectation_value(state, mpo, env.L, env.R);
        tools::common::profile::get_default_prof()["t_var"]->toc();
        double var = std::abs(H2 - E2);
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

double tools::finite::measure::energy_minus_energy_reduced(const class_tensors_finite &tensors) {
    return energy_minus_energy_reduced(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
}

double tools::finite::measure::energy(const class_tensors_finite &tensors) {
    if(tensors.measurements.energy) return tensors.measurements.energy.value();
    else {
        tensors.measurements.energy = tools::finite::measure::energy(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
        return tensors.measurements.energy.value();
    }
}

double tools::finite::measure::energy_per_site(const class_tensors_finite &tensors) {
    if(tensors.measurements.energy_per_site) return tensors.measurements.energy_per_site.value();
    else {
        tensors.measurements.energy_per_site = tools::finite::measure::energy_per_site(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
        return tensors.measurements.energy_per_site.value();
    }
}

double tools::finite::measure::energy_variance(const class_tensors_finite &tensors) {
    if(tensors.measurements.energy_variance) return tensors.measurements.energy_variance.value();
    else {
        tensors.measurements.energy_variance = tools::finite::measure::energy_variance(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
        return tensors.measurements.energy_variance.value();
    }
}

double tools::finite::measure::energy_variance_per_site(const class_tensors_finite &tensors) {
    if(tensors.measurements.energy_variance_per_site) return tensors.measurements.energy_variance_per_site.value();
    else {
        tensors.measurements.energy_variance_per_site =
            tools::finite::measure::energy_variance_per_site(*tensors.state, *tensors.model, *tensors.edges, &tensors.measurements);
        return tensors.measurements.energy_variance_per_site.value();
    }
}

double tools::finite::measure::energy_normalized(const class_tensors_finite &tensors, double emin, double emax) {
    return tools::finite::measure::energy_normalized(*tensors.state, *tensors.model, *tensors.edges, emin, emax);
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
