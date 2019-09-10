//
// Created by david on 2019-09-10.
//
#include <tools/nmspc_tools.h>
#include <state/class_finite_state.h>
#include <simulation/nmspc_settings.h>
#include <spdlog/fmt/bundled/ranges.h>
using namespace Textra;
using Scalar = class_finite_state::Scalar;

class_finite_state tools::finite::measure::reduced::get_state_with_energy_reduced_mpo (const class_finite_state & state){
    tools::log->trace("Generating reduced energy state");

    auto state_reduced = state;
    state_reduced.unset_measurements();
    state_reduced.clear_cache();
    auto   multitheta        = state_reduced.get_multitheta();
    double energy_per_site   = tools::finite::measure::multisite::energy_per_site         (state_reduced,multitheta);
//    double variance_mul      = tools::finite::measure::multisite::energy_variance_per_site(state_reduced,multitheta);
//    tools::log->trace("State before reduction: energy per site {} | log10 variance per site: multisite = {}",energy_per_site,std::log10(variance_mul));

    for (auto &mpo : state_reduced.MPO_L) mpo->set_reduced_energy(energy_per_site);
    for (auto &mpo : state_reduced.MPO_R) mpo->set_reduced_energy(energy_per_site);
    tools::finite::mps::rebuild_environments(state_reduced);
    state_reduced.unset_measurements();
    state_reduced.clear_cache();
//    energy_per_site   = tools::finite::measure::multisite::energy_per_site(state_reduced,multitheta);
//    variance_mul      = tools::finite::measure::multisite::energy_variance_per_site(state_reduced,multitheta);
//    tools::log->trace("State after reduction:  energy per site {} | log10 variance per site: multisite = {}",energy_per_site,std::log10(variance_mul));
    return state_reduced;

}

double tools::finite::measure::reduced::internal::significant_digits(double H2, double E2){
    double max_digits    = std::numeric_limits<double>::max_digits10;
    double lost_bits     = -std::log2(1.0 - std::abs(std::min(H2,E2)/std::max(H2,E2)));
    double lost_digits   = std::log10(std::pow(2.0,lost_bits));
//    tools::log->trace("Significant digits: {}",std::floor(max_digits - lost_digits));
    return digits = std::floor(max_digits - lost_digits);
}



double tools::finite::measure::reduced::energy_variance(const class_finite_state &state, const Eigen::Tensor<Scalar,3> & multitheta){
    auto state_reduced = tools::finite::measure::reduced::get_state_with_energy_reduced_mpo(state);
    double E_reduced   = tools::finite::measure::multisite::energy(state_reduced,multitheta);
    double H2_reduced  = tools::finite::measure::multisite::energy_variance(state_reduced,multitheta);
    internal::significant_digits(H2_reduced, E_reduced*E_reduced);
    return std::abs(H2_reduced - E_reduced*E_reduced);
}


double tools::finite::measure::reduced::energy_variance_per_site(const class_finite_state &state,const Eigen::Tensor<Scalar,3> & multitheta){
    return reduced::energy_variance(state,multitheta)/state.get_length();
}

double tools::finite::measure::reduced::energy_variance(const class_finite_state &state){
    return reduced::energy_variance(state,state.get_multitheta());
}


double tools::finite::measure::reduced::energy_variance_per_site(const class_finite_state &state){
    return reduced::energy_variance(state)/state.get_length();
}

