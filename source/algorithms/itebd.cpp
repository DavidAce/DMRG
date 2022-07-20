#include "itebd.h"
#include "config/settings.h"
#include "qm/time.h"
#include "tensors/model/ModelInfinite.h"
#include "tensors/site/mpo/MpoSite.h"
#include "tensors/state/StateInfinite.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/infinite/opt.h"

itebd::itebd(std::shared_ptr<h5pp::File> h5ppFile_) : AlgorithmInfinite(std::move(h5ppFile_), AlgorithmType::iTEBD) {
    tools::log->trace("Constructing class_itebd");
}

void itebd::run_preprocessing() {
    tools::log->info("Running {} preprocessing", status.algo_type_sv());
    auto t_pre = tid::tic_scope("pre");
    init_bond_dimension_limits();
    init_truncation_error_limits();
    randomize_model(); // First use of random!
    auto t_init    = tid::tic_scope("init");
    status.delta_t = std::complex<double>(settings::itebd::time_step_init_real, settings::itebd::time_step_init_imag);
    h_evn          = tensors.model->get_2site_ham_AB();
    h_odd          = tensors.model->get_2site_ham_BA();

    unitary_time_evolving_operators = qm::time::get_twosite_time_evolution_operators(status.delta_t, settings::itebd::suzuki_order, h_evn, h_odd);
    tools::log->info("Finished {} preprocessing", status.algo_type_sv());
}

void itebd::run_simulation() {
    tools::log->info("Starting {} simulation", status.algo_type_sv());
    auto t_run = tid::tic_scope("run");
    while(status.iter < settings::itebd::max_iters and status.algorithm_converged_for == 0) {
        single_TEBD_step();
        check_convergence();
        print_status_update();
        write_to_file();

        status.iter++;
        status.wall_time = tid::get_unscoped("t_tot").get_time();
        status.algo_time = t_run->get_time();
    }
}

void itebd::run_postprocessing() {
    auto t_pos = tid::tic_scope("post");
    print_status_full();
}

void itebd::single_TEBD_step() {
    /*!
     * \fn single_TEBD_step(class_superblock &state)
     * \brief infinite Time evolving block decimation.
     */
    auto t_step = tid::tic_scope("step");
    for(auto &U : unitary_time_evolving_operators) {
        Eigen::Tensor<Scalar, 3> twosite_tensor = tools::infinite::opt::time_evolve_state(*tensors.state, U);
        tensors.merge_twosite_tensor(twosite_tensor, status.bond_lim);
        if(&U != &unitary_time_evolving_operators.back()) { tensors.state->swap_AB(); }
    }
    status.phys_time += std::abs(status.delta_t);
    tensors.clear_measurements();
}

void itebd::check_convergence() {
    auto t_con = tid::tic_scope("conv");
    check_convergence_entg_entropy();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    update_bond_dimension_limit();
    check_convergence_time_step();
    if(status.entanglement_converged_for > 0 and status.variance_ham_converged_for > 0 and status.variance_mom_converged_for > 0 and
       status.bond_limit_has_reached_max and status.time_step_has_converged) {
        status.algorithm_converged_for++;
    } else
        status.algorithm_converged_for = 0;
}

void itebd::check_convergence_time_step() {
    if(std::abs(status.delta_t) <= settings::itebd::time_step_min) {
        status.time_step_has_converged = true;
    } else if(status.bond_limit_has_reached_max and status.entanglement_converged_for > 0) {
        status.delta_t                  = std::max(settings::itebd::time_step_min, std::abs(status.delta_t) * 0.5);
        unitary_time_evolving_operators = qm::time::get_twosite_time_evolution_operators(status.delta_t, settings::itebd::suzuki_order, h_evn, h_odd);
        //        state->H->update_evolution_step_size(-status.delta_t, settings::itebd::suzuki_order);
        clear_convergence_status();
    }
}
