//
// Created by david on 2019-02-14.
//

#ifndef DMRG_CLASS_SIMULATION_STATE_H
#define DMRG_CLASS_SIMULATION_STATE_H

#include <memory>

class class_simulation_state{
public:
    // Common variables
    int    iteration                      = 0; //In idmrg and itebd: iterations, in fdmrg and xdmrg: full sweeps along the chain.
    int    step                           = 0; //In fdmrg and xdmrg: how many individual moves along the chain.
    int    position                       = 0;
    long   chi_temp                       = 16;
    int    min_sweeps                     = 2 ;
    bool   simulation_has_converged       = false;
    bool   simulation_has_to_stop         = false;
    bool   bond_dimension_has_reached_max = false;
    bool   entanglement_has_converged     = false;
    bool   entanglement_has_saturated     = false;
    int    variance_mpo_saturated_for     = 0;
    bool   variance_mpo_has_converged     = false;
    bool   variance_mpo_has_saturated     = false;
    int    variance_ham_saturated_for     = 0;
    bool   variance_ham_has_converged     = false;
    bool   variance_ham_has_saturated     = false;
    int    variance_mom_saturated_for     = 0;
    bool   variance_mom_has_converged     = false;
    bool   variance_mom_has_saturated     = false;

    double energy_min     = 0;
    double energy_max     = 0;
    double energy_target  = 0;
    double energy_ubound  = 0;
    double energy_lbound  = 0;
    double energy_now     = 0;

    double phys_time               = 0;
    double delta_t                 = 0; //Make sure this one gets initialized to delta_t0!
    bool   time_step_has_converged = false;

    void clear(){*this = class_simulation_state();}
};



#endif //DMRG_CLASS_SIMULATION_STATE_H
