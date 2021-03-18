//
// Created by david on 2018-02-09.
//

#pragma once

#include "class_algorithm_finite.h"
#include <deque>
#include <physics/class_quantum_gates.h>

/*!
 * \brief Class that runs the excited-state DMRG algorithm.
 */

class class_state_finite;
class class_xdmrg : public class_algorithm_finite {
    private:
    double energy_window_growth_factor = 1.0;
    std::optional<double> variance_before_step = std::nullopt;
    std::optional<double> alpha_expansion_prev = std::nullopt; // Remember this value;
    std::vector<qm::Gate> ham_gates_1body;
    std::vector<qm::Gate> ham_gates_2body;
    std::vector<qm::Gate> ham_gates_3body;
    std::pair<std::vector<qm::Gate>,std::vector<qm::Gate>> time_gates_1site;
    std::pair<std::vector<qm::Gate>,std::vector<qm::Gate>> time_gates_2site;
    std::pair<std::vector<qm::Gate>,std::vector<qm::Gate>> time_gates_3site;


    struct OptConf {
        OptMode               optMode           = OptMode::VARIANCE;
        OptSpace              optSpace          = OptSpace::DIRECT;
        OptType               optType           = OptType::CPLX;
        OptInit               optInit           = OptInit::CURRENT_STATE;
        size_t                max_sites         = 2;
        size_t                min_sites         = 1;
        long                  max_problem_size  = 0;
        long                  problem_size      = 0;
        bool                  second_chance     = true;
        std::optional<double> alpha_expansion   = std::nullopt;
        std::array<long, 3>   problem_dims;
        std::vector<size_t>   chosen_sites;
        std::string           label;
    };
    std::vector<OptConf> get_opt_conf_list();

    public:
    // Inherit the constructor of class_algorithm_base
    using class_algorithm_finite::class_algorithm_finite;
    explicit class_xdmrg(std::shared_ptr<h5pp::File> h5ppFile_);
    void   find_energy_range();
    void   init_energy_limits(std::optional<double> energy_density_target = std::nullopt, std::optional<double> energy_density_window = std::nullopt);
    void   single_xDMRG_step(std::vector<class_xdmrg::OptConf> conf = {});
    void   randomize_into_state_in_energy_window(ResetReason reason, StateInit state_type, std::optional<std::string> sector = std::nullopt);
    void   run_task_list(std::deque<xdmrg_task> &task_list);
    void   run_preprocessing() final;
    void   resume() final;
    void   run_default_task_list() final;
    void   run_algorithm() final;
    void   check_convergence() final;
    void   update_time_step();
    void   create_hamiltonian_gates();
    void   create_time_evolution_gates();
    bool   cfg_algorithm_is_on() final;
    long   cfg_chi_lim_max() final;
    size_t cfg_print_freq() final;
    bool   cfg_chi_lim_grow() final;
    long   cfg_chi_lim_init() final;
    bool   cfg_store_wave_function() final;
};
