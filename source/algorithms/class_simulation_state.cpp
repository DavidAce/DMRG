//
// Created by david on 2019-02-14.
//

#include "class_simulation_state.h"


void class_simulation_state::clear(){*this = class_simulation_state();}


std::ostream & operator << (std::ostream& os, const class_simulation_state & sim_state){
    return os
            << std::string("iteration                      : ")  << sim_state.iteration                                         << '\n'
            << std::string("step                           : ")  << sim_state.step                                              << '\n'
            << std::string("position                       : ")  << sim_state.position                                          << '\n'
            << std::string("chi_temp                       : ")  << sim_state.chi_temp                                          << '\n'
            << std::string("chi_max                        : ")  << sim_state.chi_max                                           << '\n'
            << std::string("min_sweeps                     : ")  << sim_state.min_sweeps                                        << '\n'
            << std::string("energy_min                     : ")  << sim_state.energy_min                                        << '\n'
            << std::string("energy_max                     : ")  << sim_state.energy_max                                        << '\n'
            << std::string("energy_target                  : ")  << sim_state.energy_target                                     << '\n'
            << std::string("energy_ubound                  : ")  << sim_state.energy_ubound                                     << '\n'
            << std::string("energy_lbound                  : ")  << sim_state.energy_lbound                                     << '\n'
            << std::string("energy_now                     : ")  << sim_state.energy_now                                        << '\n'
            << std::string("energy_dens                    : ")  << sim_state.energy_dens                                       << '\n'
            << std::string("energy_dens_target             : ")  << sim_state.energy_dens_target                                << '\n'
            << std::string("energy_dens_window             : ")  << sim_state.energy_dens_window                                << '\n'
            << std::string("phys_time                      : ")  << sim_state.phys_time                                         << '\n'
            << std::string("wall_time                      : ")  << sim_state.wall_time                                         << '\n'
            << std::string("simu_time                      : ")  << sim_state.simu_time                                         << '\n'
            << std::string("delta_t                        : ")  << sim_state.delta_t                                           << '\n'
            << std::string("time_step_has_converged        : ")  << std::boolalpha << sim_state.time_step_has_converged         << '\n'
            << std::string("simulation_has_converged       : ")  << std::boolalpha << sim_state.simulation_has_converged        << '\n'
            << std::string("simulation_has_saturated       : ")  << std::boolalpha << sim_state.simulation_has_saturated        << '\n'
            << std::string("simulation_has_to_stop         : ")  << std::boolalpha << sim_state.simulation_has_to_stop          << '\n'
            << std::string("bond_dimension_has_reached_max : ")  << std::boolalpha << sim_state.bond_dimension_has_reached_max  << '\n'
            << std::string("entanglement_has_converged     : ")  << std::boolalpha << sim_state.entanglement_has_converged      << '\n'
            << std::string("entanglement_has_saturated     : ")  << std::boolalpha << sim_state.entanglement_has_saturated      << '\n'
            << std::string("variance_mpo_has_converged     : ")  << std::boolalpha << sim_state.variance_mpo_has_converged      << '\n'
            << std::string("variance_mpo_has_saturated     : ")  << std::boolalpha << sim_state.variance_mpo_has_saturated      << '\n'
            << std::string("variance_ham_has_converged     : ")  << std::boolalpha << sim_state.variance_ham_has_converged      << '\n'
            << std::string("variance_ham_has_saturated     : ")  << std::boolalpha << sim_state.variance_ham_has_saturated      << '\n'
            << std::string("variance_mom_has_converged     : ")  << std::boolalpha << sim_state.variance_mom_has_converged      << '\n'
            << std::string("variance_mom_has_saturated     : ")  << std::boolalpha << sim_state.variance_mom_has_saturated      << '\n'
            << std::string("variance_mpo_saturated_for     : ")  << sim_state.variance_mpo_saturated_for                        << '\n'
            << std::string("variance_ham_saturated_for     : ")  << sim_state.variance_ham_saturated_for                        << '\n'
            << std::string("variance_mom_saturated_for     : ")  << sim_state.variance_mom_saturated_for                        << '\n';

}