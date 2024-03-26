#pragma once
#include "config/enums.h"
#include "math/float.h"
#include <complex>
#include <cstddef>
#include <h5pp/details/h5ppVstr.h>
#include <h5pp/details/h5ppFstr.h>

class AlgorithmStatus {
    public:
    // common variables
    size_t               iter                       = 0; // In idmrg and itebd: iterations, in fdmrg and xdmrg: sweeps along the chain.
    size_t               step                       = 0; // How many dmrg steps have been taken (each step may cover multiple sites)
    long                 position                   = 0;
    int                  direction                  = 1;
    StorageEvent         event                      = StorageEvent::NONE;
    size_t               num_resets                 = 0;
    size_t               min_iters                  = 0;
    long                 bond_lim                   = 0; /*!< Current limit on bond dimension, can be increased dynamically */
    long                 bond_max                   = 0; /*!< Maximum allowable bond dimension during an algorithm run */
    long                 bond_init                  = 0; /*!< Initial limit on bond dimension when an algorithm starts */
    double               trnc_lim                   = 0; /*!< Current truncation error limit */
    double               trnc_min                   = 0; /*!< Minimum truncation error limit for this simulation */
    double               trnc_init                  = 0; /*!< Initial max limit on truncation error */
    double               energy_min                 = std::numeric_limits<double>::quiet_NaN();
    double               energy_max                 = std::numeric_limits<double>::quiet_NaN();
    double               energy_tgt                 = std::numeric_limits<double>::quiet_NaN();
    double               energy_dens                = std::numeric_limits<double>::quiet_NaN();
    double               energy_dens_target         = std::numeric_limits<double>::quiet_NaN();
    double               energy_variance_lowest     = 1;
    size_t               energy_variance_max_digits = 0;
    double               energy_variance_prec_limit = 0;
    double               env_expansion_alpha        = 0; /*!< subspace expansion factor alpha */
    double               env_expansion_variance     = 0; /*!< lowest variance when alpha was last updated */
    size_t               env_expansion_step         = 0; /*!< step when alpha was last updated */
    h5pp::fstr_t<64>     phys_time                  = {};
    double               wall_time                  = 0;
    double               algo_time                  = 0;
    h5pp::fstr_t<128>    delta_t                    = {}; // Note this is complex!! Make sure this one gets initialized to delta_t0!
    AlgorithmType        algo_type                  = AlgorithmType::ANY;
    AlgorithmStop        algo_stop                  = AlgorithmStop::NONE;
    bool                 algorithm_has_finished     = false;
    bool                 algorithm_has_succeeded    = false;
    bool                 algorithm_has_to_stop      = false;
    size_t               algorithm_has_stuck_for    = 0;
    size_t               algorithm_saturated_for    = 0;
    size_t               algorithm_converged_for    = 0;
    size_t               entanglement_converged_for = 0;
    size_t               entanglement_saturated_for = 0;
    size_t               variance_mpo_converged_for = 0;
    size_t               variance_mpo_saturated_for = 0;
    size_t               variance_ham_converged_for = 0;
    size_t               variance_ham_saturated_for = 0;
    size_t               variance_mom_converged_for = 0;
    size_t               variance_mom_saturated_for = 0;
    bool                 bond_limit_has_reached_max = false;
    bool                 trnc_limit_has_reached_min = false;
    bool                 spin_parity_has_converged  = false;
    bool                 time_step_has_converged    = false;
    bool                 fes_is_running             = false;

    void                           clear();
    void                           reset();
    [[nodiscard]] std::string_view algo_type_sv() const;
    [[nodiscard]] std::string_view algo_stop_sv() const;
    [[nodiscard]] bool             operator==(const AlgorithmStatus &s) const;
};
