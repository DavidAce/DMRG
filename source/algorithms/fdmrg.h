#pragma once
#include "AlgorithmFinite.h"
#include <deque>
class class_h5table_measurements_finite;
namespace tools::finite::opt {
    struct OptMeta;
}
/*!
// * \brief Class that runs the finite DMRG algorithm.
 */

class StateFinite;
class fdmrg : public AlgorithmFinite {
    private:
    std::string_view get_state_name() const;
    public:
    using OptMeta = tools::finite::opt::OptMeta;
    // Inherit the constructor of class_algorithm_base
    using AlgorithmFinite::AlgorithmFinite;
    fdmrg();
    explicit fdmrg(std::shared_ptr<h5pp::File> h5file_);
    OptRitz ritz = OptRitz::SR;
    void    resume() final;
    void    run_task_list(std::deque<fdmrg_task> &task_list);
    void    run_default_task_list() final;
    void    run_preprocessing() final;
    void    run_algorithm() final;
    void    run_fes_analysis() final;
    void    update_state() final;
    void    check_convergence() final;
};
