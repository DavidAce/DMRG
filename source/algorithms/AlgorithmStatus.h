#pragma once
#include <complex>
#include <config/enums.h>
#include <cstddef>

class AlgorithmStatus {
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
    size_t               energy_variance_max_digits  = 0;
    double               energy_variance_prec_limit  = 0;
    double               sub_expansion_alpha         = 0; /*!< subspace expansion factor alpha */
    double               sub_expansion_variance      = 0; /*!< lowest variance when alpha was last updated */
    size_t               sub_expansion_step          = 0; /*!< step when alpha was last updated */
    double               phys_time                   = 0;
    double               wall_time                   = 0;
    double               algo_time                   = 0;
    std::complex<double> delta_t                     = 0; // Note this is complex!! Make sure this one gets initialized to delta_t0!
    AlgorithmType        algo_type                   = AlgorithmType::ANY;
    AlgorithmStop        algo_stop                   = AlgorithmStop::NONE;
    bool                 algorithm_has_finished      = false;
    bool                 algorithm_has_succeeded     = false;
    bool                 algorithm_has_to_stop       = false;
    size_t               algorithm_has_stuck_for     = 0;
    size_t               algorithm_saturated_for     = 0;
    size_t               algorithm_converged_for     = 0;
    size_t               entanglement_converged_for  = 0;
    size_t               entanglement_saturated_for  = 0;
    size_t               variance_mpo_converged_for  = 0;
    size_t               variance_mpo_saturated_for  = 0;
    size_t               variance_ham_converged_for  = 0;
    size_t               variance_ham_saturated_for  = 0;
    size_t               variance_mom_converged_for  = 0;
    size_t               variance_mom_saturated_for  = 0;
    bool                 chi_lim_has_reached_chi_max = false;
    bool                 spin_parity_has_converged   = false;
    bool                 time_step_has_converged     = false;
    bool                 fes_is_running              = false;
    void                 clear() {
        AlgorithmType algo_type_ = algo_type; // Retain the algorithm type
        *this                    = AlgorithmStatus();
        algo_type                = algo_type_;
    }
    void reset() {
        // Keeps some data for simulations that follow
        auto status = *this;
        clear();
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
        algo_type            = status.algo_type;
    }
    [[nodiscard]] std::string_view algo_type_sv() const { return enum2sv(algo_type); }
    [[nodiscard]] std::string      algo_type_str() const { return std::string(algo_type_sv()); }
    [[nodiscard]] std::string_view algo_stop_sv() const { return enum2sv(algo_stop); }
    [[nodiscard]] std::string      algo_stop_str() const { return std::string(algo_stop_sv()); }
};
