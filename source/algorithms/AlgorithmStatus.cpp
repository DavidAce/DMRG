#include "AlgorithmStatus.h"

void AlgorithmStatus::clear() {
    AlgorithmType algo_type_ = algo_type; // Retain the algorithm type
    *this                    = AlgorithmStatus();
    algo_type                = algo_type_;
}
void AlgorithmStatus::reset() {
    // Keeps some data for simulations that follow
    auto status = *this;
    clear();
    min_iters          = status.min_iters;
    bond_max           = status.bond_max;
    bond_init          = status.bond_init;
    bond_lim           = status.bond_lim;
    trnc_min           = status.trnc_min;
    trnc_init          = status.trnc_init;
    trnc_lim           = status.trnc_lim;
    energy_min         = status.energy_min;
    energy_max         = status.energy_max;
    energy_tgt         = status.energy_tgt;
    energy_ulim        = status.energy_ulim;
    energy_llim        = status.energy_llim;
    energy_dens        = status.energy_dens;
    energy_dens_target = status.energy_dens_target;
    energy_dens_window = status.energy_dens_window;
    algo_type          = status.algo_type;
}
std::string_view AlgorithmStatus::algo_type_sv() const { return enum2sv(algo_type); }
std::string      AlgorithmStatus::algo_type_str() const { return std::string(algo_type_sv()); }
std::string_view AlgorithmStatus::algo_stop_sv() const { return enum2sv(algo_stop); }
std::string      AlgorithmStatus::algo_stop_str() const { return std::string(algo_stop_sv()); }

bool AlgorithmStatus::operator==(const AlgorithmStatus &s) const {
    return
        /* clang-format off */
        this->iter                          == s.iter and
        this->step                          == s.step and
        this->position                      == s.position and
        this->direction                     == s.direction and
        this->event                         == s.event and
        this->num_resets                    == s.num_resets and
        this->min_iters                     == s.min_iters and
        this->bond_max                      == s.bond_max and
        this->bond_init                     == s.bond_init and
        this->bond_lim                      == s.bond_lim and
        this->trnc_min                      == s.trnc_min and
        this->trnc_init                     == s.trnc_init and
        this->trnc_lim                      == s.trnc_lim and
        this->energy_min                    == s.energy_min and
        this->energy_max                    == s.energy_max and
        this->energy_tgt                    == s.energy_tgt and
        this->energy_ulim                   == s.energy_ulim and
        this->energy_llim                   == s.energy_llim and
        this->energy_dens                   == s.energy_dens and
        this->energy_dens_target            == s.energy_dens_target and
        this->energy_dens_window            == s.energy_dens_window and
        this->energy_variance_lowest        == s.energy_variance_lowest and
        this->energy_variance_max_digits    == s.energy_variance_max_digits and
        this->energy_variance_prec_limit    == s.energy_variance_prec_limit and
        this->env_expansion_alpha           == s.env_expansion_alpha and
        this->env_expansion_variance        == s.env_expansion_variance and
        this->env_expansion_step            == s.env_expansion_step and
        this->phys_time                     == s.phys_time and
        this->wall_time                     == s.wall_time and
        this->algo_time                     == s.algo_time and
        this->delta_t                       == s.delta_t and
        this->algo_type                     == s.algo_type and
        this->algo_stop                     == s.algo_stop and
        this->algorithm_has_finished        == s.algorithm_has_finished and
        this->algorithm_has_succeeded       == s.algorithm_has_succeeded and
        this->algorithm_has_to_stop         == s.algorithm_has_to_stop and
        this->algorithm_has_stuck_for       == s.algorithm_has_stuck_for and
        this->algorithm_saturated_for       == s.algorithm_saturated_for and
        this->algorithm_converged_for       == s.algorithm_converged_for and
        this->entanglement_converged_for    == s.entanglement_converged_for and
        this->entanglement_saturated_for    == s.entanglement_saturated_for and
        this->variance_mpo_converged_for    == s.variance_mpo_converged_for and
        this->variance_mpo_saturated_for    == s.variance_mpo_saturated_for and
        this->variance_ham_converged_for    == s.variance_ham_converged_for and
        this->variance_ham_saturated_for    == s.variance_ham_saturated_for and
        this->variance_mom_converged_for    == s.variance_mom_converged_for and
        this->variance_mom_saturated_for    == s.variance_mom_saturated_for and
        this->bond_limit_has_reached_max    == s.bond_limit_has_reached_max and
        this->trnc_limit_has_reached_min    == s.trnc_limit_has_reached_min and
        this->spin_parity_has_converged     == s.spin_parity_has_converged and
        this->time_step_has_converged       == s.time_step_has_converged and
        this->fes_is_running                == s.fes_is_running;
    /* clang-format on */
}