//
// Created by david on 2019-02-01.
//

#include <general/nmspc_tensor_extra.h>
#include <tensors/class_tensors_infinite.h>
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/contraction.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/infinite/measure.h>

void tools::infinite::measure::do_all_measurements(const class_tensors_infinite &tensors) {
    tensors.measurements.length                       = tools::infinite::measure::length(tensors);
    tensors.measurements.energy_mpo                   = tools::infinite::measure::energy_mpo(tensors);
    tensors.measurements.energy_per_site_mpo          = tools::infinite::measure::energy_per_site_mpo(tensors);
    tensors.measurements.energy_variance_mpo          = tools::infinite::measure::energy_variance_mpo(tensors);
    tensors.measurements.energy_per_site_ham          = tools::infinite::measure::energy_per_site_ham(tensors);
    tensors.measurements.energy_per_site_mom          = tools::infinite::measure::energy_per_site_mom(tensors);
    tensors.measurements.energy_variance_per_site_mpo = tools::infinite::measure::energy_variance_per_site_mpo(tensors);
    tensors.measurements.energy_variance_per_site_ham = tools::infinite::measure::energy_variance_per_site_ham(tensors);
    tensors.measurements.energy_variance_per_site_mom = tools::infinite::measure::energy_variance_per_site_mom(tensors);
}

void tools::infinite::measure::do_all_measurements(const class_state_infinite &state) {
    state.measurements.norm                 = tools::infinite::measure::norm(state);
    state.measurements.bond_dimension       = tools::infinite::measure::bond_dimension(state);
    state.measurements.entanglement_entropy = tools::infinite::measure::entanglement_entropy(state);
    state.measurements.truncation_error     = tools::infinite::measure::truncation_error(state);
}

size_t tools::infinite::measure::length(const class_tensors_infinite &tensors) { return tensors.edges->get_length(); }
size_t tools::infinite::measure::length(const class_edges_infinite &edges) { return edges.get_length(); }

double tools::infinite::measure::norm(const class_state_infinite &state) {
    if(state.measurements.norm) return state.measurements.norm.value();
    Eigen::Tensor<Scalar, 0> norm = state.get_2site_mps().contract(state.get_2site_mps().conjugate(), Textra::idx({0, 1, 2}, {0, 1, 2}));
    return std::abs(norm(0));
}

long tools::infinite::measure::bond_dimension(const class_state_infinite &state) {
    if(state.measurements.bond_dimension) return state.measurements.bond_dimension.value();
    state.measurements.bond_dimension = state.chiC();
    return state.measurements.bond_dimension.value();
}

double tools::infinite::measure::truncation_error(const class_state_infinite &state) {
    if(state.measurements.truncation_error) return state.measurements.truncation_error.value();
    state.measurements.truncation_error = state.get_truncation_error();
    return state.measurements.truncation_error.value();
}

double tools::infinite::measure::entanglement_entropy(const class_state_infinite &state) {
    if(state.measurements.entanglement_entropy) return state.measurements.entanglement_entropy.value();
    tools::common::profile::get_default_prof()["t_ent"]->tic();
    const auto &             LC = state.LC();
    Eigen::Tensor<Scalar, 0> SA = -LC.square().contract(LC.square().log().eval(), Textra::idx({0}, {0}));
    state.measurements.entanglement_entropy = std::real(SA(0));
    tools::common::profile::get_default_prof()["t_ent"]->toc();
    return state.measurements.entanglement_entropy.value();
}

template<typename state_or_mps_type>
double tools::infinite::measure::energy_minus_energy_reduced(const state_or_mps_type &state, const class_model_infinite &model,
                                                             const class_edges_infinite &edges) {
    if constexpr(std::is_same_v<state_or_mps_type, class_state_infinite>) {
        return tools::infinite::measure::energy_minus_energy_reduced(state.get_2site_mps(), model, edges);
    } else {
        tools::log->trace("Measuring energy mpo");
        const auto &mpo = model.get_2site_mpo();
        const auto &env = edges.get_ene_blk();
        tools::common::profile::get_default_prof()["t_ene"]->tic();
        double e_minus_ered = tools::common::contraction::expectation_value(state, mpo, env.L, env.R);
        tools::common::profile::get_default_prof()["t_ene"]->toc();
        return e_minus_ered;
    }
}

template double tools::infinite::measure::energy_minus_energy_reduced(const class_state_infinite &, const class_model_infinite &model,
                                                                      const class_edges_infinite &edges);
template double tools::infinite::measure::energy_minus_energy_reduced(const Eigen::Tensor<Scalar, 3> &, const class_model_infinite &model,
                                                                      const class_edges_infinite &edges);

template<typename state_or_mps_type>
double tools::infinite::measure::energy_mpo(const state_or_mps_type &state, const class_model_infinite &model, const class_edges_infinite &edges) {
    if constexpr(std::is_same_v<state_or_mps_type, class_state_infinite>) return tools::infinite::measure::energy_mpo(state.get_2site_mps(), model, edges);
    else
        return tools::infinite::measure::energy_minus_energy_reduced(state, model, edges) +
               model.get_energy_per_site_reduced() * static_cast<double>(edges.get_length());
}

template double tools::infinite::measure::energy_mpo(const class_state_infinite &, const class_model_infinite &model, const class_edges_infinite &edges);
template double tools::infinite::measure::energy_mpo(const Eigen::Tensor<Scalar, 3> &, const class_model_infinite &model, const class_edges_infinite &edges);

template<typename state_or_mps_type>
double tools::infinite::measure::energy_per_site_mpo(const state_or_mps_type &state, const class_model_infinite &model, const class_edges_infinite &edges) {
    tools::log->warn("energy_per_site_mpo: CHECK DIVISION");
    return tools::infinite::measure::energy_mpo(state, model, edges) / static_cast<double>(2);
}

template double tools::infinite::measure::energy_per_site_mpo(const class_state_infinite &, const class_model_infinite &model,
                                                              const class_edges_infinite &edges);
template double tools::infinite::measure::energy_per_site_mpo(const Eigen::Tensor<Scalar, 3> &, const class_model_infinite &model,
                                                              const class_edges_infinite &edges);

template<typename state_or_mps_type>
double tools::infinite::measure::energy_variance_mpo(const state_or_mps_type &state, const class_model_infinite &model, const class_edges_infinite &edges) {
    // Depending on whether the mpo's are reduced or not we get different formulas.
    // If mpo's are reduced:
    //      Var H = <(H-E_red)^2> - <(H-E_red)>^2 = <H^2> - 2<H>E_red + E_red^2 - (<H> - E_red) ^2
    //                                            = H2    - 2*E*E_red + E_red^2 - E^2 + 2*E*E_red - E_red^2
    //                                            = H2    - E^2
    //      so Var H = <(H-E_red)^2> - energy_minus_energy_reduced^2 = H2 - ~0
    //      where H2 is computed with reduced mpo's. Note that ~0 is not exactly zero
    //      because E_red != E necessarily (though they are supposed to be very close)
    // Else:
    //      Var H = <(H - 0)^2> - <H - 0>^2 = H2 - E^2
    if constexpr(std::is_same_v<state_or_mps_type, class_state_infinite>) {
        auto var = tools::infinite::measure::energy_variance_mpo(state.get_2site_mps(), model, edges);
        if(var < state.lowest_recorded_variance) state.lowest_recorded_variance = var;
        return var;
    } else {
        tools::log->trace("Measuring energy variance mpo");
        double energy = 0;
        if(model.is_reduced()) energy = tools::infinite::measure::energy_minus_energy_reduced(state, model, edges);
        else
            energy = tools::infinite::measure::energy_mpo(state, model, edges);
        double      E2  = energy * energy;
        const auto &mpo = model.get_2site_mpo();
        const auto &env = edges.get_var_blk();
        tools::log->trace("Measuring energy variance mpo");
        tools::common::profile::get_default_prof()["t_ene"]->tic();
        double H2 = tools::common::contraction::expectation_value(state, mpo, env.L, env.R);
        tools::common::profile::get_default_prof()["t_ene"]->toc();
        double var = std::abs(H2 - E2);
        return var;
    }
}

template double tools::infinite::measure::energy_variance_mpo(const class_state_infinite &, const class_model_infinite &model,
                                                              const class_edges_infinite &edges);
template double tools::infinite::measure::energy_variance_mpo(const Eigen::Tensor<Scalar, 3> &, const class_model_infinite &model,
                                                              const class_edges_infinite &edges);

template<typename state_or_mps_type>
double tools::infinite::measure::energy_variance_per_site_mpo(const state_or_mps_type &state, const class_model_infinite &model,
                                                              const class_edges_infinite &edges) {
    tools::log->warn("energy_per_site_mpo: CHECK DIVISION");
    return tools::infinite::measure::energy_variance_mpo(state, model, edges) / static_cast<double>(2);
}

template double tools::infinite::measure::energy_variance_per_site_mpo(const class_state_infinite &, const class_model_infinite &model,
                                                                       const class_edges_infinite &edges);
template double tools::infinite::measure::energy_variance_per_site_mpo(const Eigen::Tensor<Scalar, 3> &, const class_model_infinite &model,
                                                                       const class_edges_infinite &edges);

double tools::infinite::measure::energy_mpo(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_mpo) return tensors.measurements.energy_mpo.value();
    tensors.measurements.energy_mpo = tools::infinite::measure::energy_mpo(*tensors.state, *tensors.model, *tensors.edges);
    return tensors.measurements.energy_mpo.value();
}

double tools::infinite::measure::energy_per_site_mpo(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_per_site_mpo) return tensors.measurements.energy_per_site_mpo.value();
    tools::log->warn("energy_per_site_mpo: CHECK DIVISION");
    auto L                                   = tools::infinite::measure::length(tensors);
    tensors.measurements.energy_per_site_mpo = tools::infinite::measure::energy_mpo(tensors) / static_cast<double>(L);
    return tensors.measurements.energy_per_site_mpo.value();
}

double tools::infinite::measure::energy_variance_mpo(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_variance_mpo) return tensors.measurements.energy_variance_mpo.value();
    tensors.measurements.energy_variance_mpo = tools::infinite::measure::energy_variance_mpo(*tensors.state, *tensors.model, *tensors.edges);
    return tensors.measurements.energy_variance_mpo.value();
}

double tools::infinite::measure::energy_variance_per_site_mpo(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_variance_per_site_mpo) return tensors.measurements.energy_variance_per_site_mpo.value();
    auto L                                            = tools::infinite::measure::length(tensors);
    tensors.measurements.energy_variance_per_site_mpo = tools::infinite::measure::energy_variance_mpo(tensors) / static_cast<double>(L);
    return tensors.measurements.energy_variance_per_site_mpo.value();
}

double tools::infinite::measure::energy_mpo(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_infinite &tensors) {
    return tools::infinite::measure::energy_mpo(mps, *tensors.model, *tensors.edges);
}
double tools::infinite::measure::energy_per_site_mpo(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_infinite &tensors) {
    return tools::infinite::measure::energy_per_site_mpo(mps, *tensors.model, *tensors.edges);
}
double tools::infinite::measure::energy_variance_mpo(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_infinite &tensors) {
    return tools::infinite::measure::energy_variance_mpo(mps, *tensors.model, *tensors.edges);
}
double tools::infinite::measure::energy_variance_per_site_mpo(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_infinite &tensors) {
    return tools::infinite::measure::energy_variance_per_site_mpo(mps, *tensors.model, *tensors.edges);
}
