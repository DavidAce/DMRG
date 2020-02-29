//
// Created by david on 2019-02-14.
//

#pragma once


#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <hdf5.h>
#include <hdf5_hl.h>

class class_simulation_status{
    public:
    // common variables
    size_t iteration                      = 0; //In idmrg and itebd: iterations, in fdmrg and xdmrg: full sweeps along the chain.
    size_t step                           = 0; //How many dmrg steps have been taken (each step may cover multiple sites)
    size_t position                       = 0;
    size_t moves                          = 0; //In fdmrg and xdmrg: how many individual moves along the chain.
    size_t num_resets                     = 0;
    size_t num_states                     = 0;  /*!< xDMRG can produce several states per disorder realization. This counts states produced */
    size_t min_sweeps                     = 0 ;
    long   chi_max                        = 0;
    long   chi_lim                        = 0;
    double energy_min                     = 0;
    double energy_max                     = 0;
    double energy_target                  = 0;
    double energy_ubound                  = 0;
    double energy_lbound                  = 0;
    double energy_dens                    = 0;
    double energy_dens_target             = 0;
    double energy_dens_window             = 0;
    double phys_time                      = 0;
    double wall_time                      = 0;
    double simu_time                      = 0;
    double delta_t                        = 0; //Make sure this one gets initialized to delta_t0!
    size_t simulation_has_stuck_for       = 0;
    size_t entanglement_saturated_for     = 0;
    size_t variance_mpo_saturated_for     = 0;
    size_t variance_ham_saturated_for     = 0;
    size_t variance_mom_saturated_for     = 0;
    bool   simulation_has_converged       = false;
    bool   simulation_has_saturated       = false;
    bool   simulation_has_succeeded       = false;
    bool   simulation_has_got_stuck       = false;
    bool   simulation_has_to_stop         = false;
    bool   chi_lim_has_reached_chi_max    = false;
    bool   entanglement_has_converged     = false;
    bool   entanglement_has_saturated     = false;
    bool   variance_mpo_has_converged     = false;
    bool   variance_mpo_has_saturated     = false;
    bool   variance_ham_has_converged     = false;
    bool   variance_ham_has_saturated     = false;
    bool   variance_mom_has_converged     = false;
    bool   variance_mom_has_saturated     = false;
    bool   time_step_has_converged        = false;

    void clear();
    friend std::ostream& operator <<(std::ostream& os, const class_simulation_status & sim_status);
};


