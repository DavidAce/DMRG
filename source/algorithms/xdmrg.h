#pragma once

#include "AlgorithmFinite.h"
#include <deque>
#include <qm/gate.h>

namespace tools::finite::opt {
    class opt_mps;
    struct OptMeta;
}
/*!
 * \brief Class that runs the excited-state DMRG algorithm.
 */

class StateFinite;
class xdmrg : public AlgorithmFinite {
    using OptMeta = tools::finite::opt::OptMeta;

    private:
    double                                                  energy_window_growth_factor = 1.0;
    std::optional<double>                                   variance_before_step        = std::nullopt;
    std::vector<qm::Gate>                                   ham_gates_1body;
    std::vector<qm::Gate>                                   ham_gates_2body;
    std::vector<qm::Gate>                                   ham_gates_3body;
    std::pair<std::vector<qm::Gate>, std::vector<qm::Gate>> time_gates_1site;
    std::pair<std::vector<qm::Gate>, std::vector<qm::Gate>> time_gates_2site;
    std::pair<std::vector<qm::Gate>, std::vector<qm::Gate>> time_gates_3site;

    [[nodiscard]] OptMeta     get_opt_meta();

    public:
    // Inherit the constructor of class_algorithm_base
    using AlgorithmFinite::AlgorithmFinite;
    explicit xdmrg(std::shared_ptr<h5pp::File> h5ppFile_);
    void     find_energy_range();
    void     init_energy_limits(std::optional<double> energy_density_target = std::nullopt, std::optional<double> energy_density_window = std::nullopt);
    void     initialize_state_in_energy_window(ResetReason reason, StateInit state_type, const std::optional<std::string> &sector = std::nullopt);
    void     run_task_list(std::deque<xdmrg_task> &task_list);
    void     run_preprocessing() final;
    void     run_default_task_list() final;
    void     run_algorithm() final;
    void     run_fes_analysis() final;
    void     update_state() final;
    void     resume() final;
    void     check_convergence() final;
    void     update_time_step();
    void     create_hamiltonian_gates();
    void     create_time_evolution_gates();
};
