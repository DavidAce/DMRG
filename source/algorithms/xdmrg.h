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
    double                energy_window_growth_factor = 1.0;
    std::optional<double> variance_before_step        = std::nullopt;

    [[nodiscard]] OptMeta get_opt_meta();

    public:
    // Inherit the constructor of class_algorithm_base
    using AlgorithmFinite::AlgorithmFinite;
    explicit xdmrg(std::shared_ptr<h5pp::File> h5ppFile_);
    void     find_energy_range();
    void     init_energy_target(std::optional<double> energy_density_target = std::nullopt);
    void     run_task_list(std::deque<xdmrg_task> &task_list);
    void     run_preprocessing() final;
    void     run_default_task_list() final;
    void     run_algorithm() final;
    void     update_state() final;
    void     resume() final;
    void     check_convergence() final;
    void     update_time_step();
    void     set_energy_shift_mpo() final;
};
