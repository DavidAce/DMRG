//
// Created by david on 2019-02-14.
//

#include "class_simulation_status.h"


void class_simulation_status::clear(){*this = class_simulation_status();}

//void class_simulation_status::get_all() {
//
//    return std::array
//
//
//}


std::ostream & operator << (std::ostream& os, const class_simulation_status & sim_status){
    return os
            << std::string("iteration                      : ") << sim_status.iteration << '\n'
            << std::string("step                           : ") << sim_status.step << '\n'
            << std::string("position                       : ") << sim_status.position << '\n'
            << std::string("moves                          : ") << sim_status.moves << '\n'
            << std::string("num_resets                     : ") << sim_status.num_resets << '\n'
            << std::string("chi_max                        : ") << sim_status.chi_max << '\n'
            << std::string("chi_lim                        : ") << sim_status.chi_lim << '\n'
            << std::string("min_sweeps                     : ")  << sim_status.min_sweeps                                        << '\n'
            << std::string("energy_min                     : ")  << sim_status.energy_min                                        << '\n'
            << std::string("energy_max                     : ")  << sim_status.energy_max                                        << '\n'
            << std::string("energy_target                  : ")  << sim_status.energy_target                                     << '\n'
            << std::string("energy_ubound                  : ")  << sim_status.energy_ubound                                     << '\n'
            << std::string("energy_lbound                  : ")  << sim_status.energy_lbound                                     << '\n'
            << std::string("energy_dens                    : ")  << sim_status.energy_dens                                       << '\n'
            << std::string("energy_dens_target             : ")  << sim_status.energy_dens_target                                << '\n'
            << std::string("energy_dens_window             : ")  << sim_status.energy_dens_window                                << '\n'
            << std::string("phys_time                      : ") << sim_status.phys_time << '\n'
            << std::string("wall_time                      : ") << sim_status.wall_time << '\n'
            << std::string("simu_time                      : ") << sim_status.simu_time << '\n'
            << std::string("delta_t                        : ") << sim_status.delta_t << '\n'
            << std::string("time_step_has_converged        : ") << std::boolalpha << sim_status.time_step_has_converged         << '\n'
            << std::string("simulation_has_converged       : ") << std::boolalpha << sim_status.simulation_has_converged        << '\n'
            << std::string("simulation_has_saturated       : ") << std::boolalpha << sim_status.simulation_has_saturated        << '\n'
            << std::string("simulation_has_succeeded       : ") << std::boolalpha << sim_status.simulation_has_succeeded        << '\n'
            << std::string("simulation_has_got_stuck       : ") << std::boolalpha << sim_status.simulation_has_got_stuck        << '\n'
            << std::string("simulation_has_stuck_for       : ") << sim_status.simulation_has_stuck_for                          << '\n'
            << std::string("simulation_has_to_stop         : ") << std::boolalpha << sim_status.simulation_has_to_stop          << '\n'
            << std::string("chi_lim_has_reached_chi_max    : ") << std::boolalpha << sim_status.chi_lim_has_reached_chi_max     << '\n'
            << std::string("entanglement_has_converged     : ") << std::boolalpha << sim_status.entanglement_has_converged      << '\n'
            << std::string("entanglement_has_saturated     : ") << std::boolalpha << sim_status.entanglement_has_saturated      << '\n'
            << std::string("variance_mpo_has_converged     : ") << std::boolalpha << sim_status.variance_mpo_has_converged      << '\n'
            << std::string("variance_mpo_has_saturated     : ") << std::boolalpha << sim_status.variance_mpo_has_saturated      << '\n'
            << std::string("variance_ham_has_converged     : ") << std::boolalpha << sim_status.variance_ham_has_converged      << '\n'
            << std::string("variance_ham_has_saturated     : ") << std::boolalpha << sim_status.variance_ham_has_saturated      << '\n'
            << std::string("variance_mom_has_converged     : ") << std::boolalpha << sim_status.variance_mom_has_converged      << '\n'
            << std::string("variance_mom_has_saturated     : ") << std::boolalpha << sim_status.variance_mom_has_saturated      << '\n'
            << std::string("entanglement_saturated_for     : ") << sim_status.entanglement_saturated_for                        << '\n'
            << std::string("variance_mpo_saturated_for     : ") << sim_status.variance_mpo_saturated_for                        << '\n'
            << std::string("variance_ham_saturated_for     : ") << sim_status.variance_ham_saturated_for                        << '\n'
            << std::string("variance_mom_saturated_for     : ") << sim_status.variance_mom_saturated_for                        << '\n';

}