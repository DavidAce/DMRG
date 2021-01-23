//
// Created by david on 2019-02-14.
//

#pragma once
#include <cstddef>
#include <complex>
class class_algorithm_status {
    public:
    // common variables
    size_t               iter                        = 0; // In idmrg and itebd: iterations, in fdmrg and xdmrg: full sweeps along the chain.
    size_t               step                        = 0; // How many dmrg steps have been taken (each step may cover multiple sites)
    long                 position                    = 0;
    int                  direction                   = 1;
    size_t               num_resets                  = 0;
    size_t               min_iters                   = 0;
    long                 chi_lim_max                 = 0; /*!< Maximum allowable bond dimension during an algorithm run */
    long                 chi_lim_init                = 0; /*!< Initial limit on bond dimension when an algorithm starts */
    long                 chi_lim                     = 0; /*!< Current limit on bond dimension, can be increased dynamically */
    double               energy_min_per_site         = 0;
    double               energy_max_per_site         = 0;
    double               energy_tgt_per_site         = 0;
    double               energy_ulim_per_site        = 0;
    double               energy_llim_per_site        = 0;
    double               energy_dens                 = 0;
    double               energy_dens_target          = 0;
    double               energy_dens_window          = 0;
    double               energy_variance_lowest      = 1;
    double               phys_time                   = 0;
    double               wall_time                   = 0;
    double               algo_time                   = 0;
    std::complex<double> delta_t                     = 0; // Note this is complex!! Make sure this one gets initialized to delta_t0!
    size_t               algorithm_has_stuck_for     = 0;
    size_t               entanglement_saturated_for  = 0;
    size_t               variance_mpo_saturated_for  = 0;
    size_t               variance_ham_saturated_for  = 0;
    size_t               variance_mom_saturated_for  = 0;
    bool                 algorithm_has_finished      = false;
    bool                 algorithm_has_converged     = false;
    bool                 algorithm_has_saturated     = false;
    bool                 algorithm_has_succeeded     = false;
    bool                 algorithm_has_got_stuck     = false;
    bool                 algorithm_has_to_stop       = false;
    bool                 chi_lim_has_reached_chi_max = false;
    bool                 entanglement_has_converged  = false;
    bool                 entanglement_has_saturated  = false;
    bool                 variance_mpo_has_converged  = false;
    bool                 variance_mpo_has_saturated  = false;
    bool                 variance_ham_has_converged  = false;
    bool                 variance_ham_has_saturated  = false;
    bool                 variance_mom_has_converged  = false;
    bool                 variance_mom_has_saturated  = false;
    bool                 time_step_has_converged     = false;
    void                 clear() { *this = class_algorithm_status(); }
    void                 reset() {
        // Keeps some data for simulations that follow
        auto status = *this;
        clear();
        energy_min_per_site  = status.energy_min_per_site;
        energy_max_per_site  = status.energy_max_per_site;
        min_iters            = status.min_iters;
        chi_lim_max          = status.chi_lim_max;
        chi_lim_init         = status.chi_lim_init;
        chi_lim              = status.chi_lim;
        energy_min_per_site  = status.energy_min_per_site;
        energy_max_per_site  = status.energy_max_per_site;
        energy_tgt_per_site  = status.energy_tgt_per_site;
        energy_ulim_per_site = status.energy_ulim_per_site;
        energy_llim_per_site = status.energy_llim_per_site;
        energy_dens          = status.energy_dens;
        energy_dens_target   = status.energy_dens_target;
        energy_dens_window   = status.energy_dens_window;
    }
};
