//
// Created by david on 2019-02-14.
//

#ifndef DMRG_CLASS_SIMULATION_STATE_H
#define DMRG_CLASS_SIMULATION_STATE_H

#include <memory>
#include <string>
#include <iostream>

class class_simulation_state{
public:
    // Common variables
    int    iteration                      = 0; //In idmrg and itebd: iterations, in fdmrg and xdmrg: full sweeps along the chain.
    int    step                           = 0; //In fdmrg and xdmrg: how many individual moves along the chain.
    int    position                       = 0;
    long   chi_temp                       = 8;
    long   chi_max                        = 8;
    int    min_sweeps                     = 2 ;
    double energy_min                     = 0;
    double energy_max                     = 0;
    double energy_target                  = 0;
    double energy_ubound                  = 0;
    double energy_lbound                  = 0;
    double energy_now                     = 0;
    double energy_dens                    = 0;
    double phys_time                      = 0;
    double delta_t                        = 0; //Make sure this one gets initialized to delta_t0!
    bool   time_step_has_converged        = false;
    bool   simulation_has_converged       = false;
    bool   simulation_has_to_stop         = false;
    bool   bond_dimension_has_reached_max = false;
    bool   entanglement_has_converged     = false;
    bool   entanglement_has_saturated     = false;
    bool   variance_mpo_has_converged     = false;
    bool   variance_mpo_has_saturated     = false;
    bool   variance_ham_has_converged     = false;
    bool   variance_ham_has_saturated     = false;
    bool   variance_mom_has_converged     = false;
    bool   variance_mom_has_saturated     = false;
    int    variance_mpo_saturated_for     = 0;
    int    variance_ham_saturated_for     = 0;
    int    variance_mom_saturated_for     = 0;


    void clear(){*this = class_simulation_state();}
    friend std::ostream& operator <<(std::ostream& os, class_simulation_state const& sim_state)
    {


        return os << std::string("iteration                      : ")  << sim_state.iteration                                         << '\n'
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
                  << std::string("phys_time                      : ")  << sim_state.phys_time                                         << '\n'
                  << std::string("delta_t                        : ")  << sim_state.delta_t                                           << '\n'
                  << std::string("time_step_has_converged        : ")  << std::boolalpha << sim_state.time_step_has_converged         << '\n'
                  << std::string("simulation_has_converged       : ")  << std::boolalpha << sim_state.simulation_has_converged        << '\n'
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


};



#endif //DMRG_CLASS_SIMULATION_STATE_H
