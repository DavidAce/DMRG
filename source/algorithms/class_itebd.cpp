//
// Created by david on 2018-01-18.
//
#include "class_itebd.h"
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <physics/nmspc_quantum_mechanics.h>
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/infinite/opt.h>

using namespace std;
using namespace Textra;

class_itebd::class_itebd(std::shared_ptr<h5pp::File> h5ppFile_) : class_algorithm_infinite(std::move(h5ppFile_), AlgorithmType::iTEBD) {
    tools::log->trace("Constructing class_itebd");
}

void class_itebd::run_preprocessing() {
    tools::log->info("Running {} preprocessing", algo_name);
    auto t_pre = tools::common::profile::prof[algo_type]["t_pre"]->tic_token();
    init_bond_dimension_limits();
    randomize_model(); // First use of random!
    status.delta_t = std::complex<double>(settings::itebd::time_step_init_real, settings::itebd::time_step_init_imag);
    h_evn          = tensors.model->get_2site_ham_AB();
    h_odd          = tensors.model->get_2site_ham_BA();

    unitary_time_evolving_operators = qm::timeEvolution::get_twosite_time_evolution_operators(status.delta_t, settings::itebd::suzuki_order, h_evn, h_odd);
    tools::log->info("Finished {} preprocessing", algo_name);
}

void class_itebd::run_simulation() {
    tools::log->info("Starting {} simulation", algo_name);
    auto t_sim = tools::common::profile::prof[algo_type]["t_sim"]->tic_token();
    while(status.iter < settings::itebd::max_iters and not status.algorithm_has_converged) {
        single_TEBD_step();
        write_to_file();
        copy_from_tmp();
        print_status_update();
        check_convergence();
        status.iter++;
    }
}

void class_itebd::run_postprocessing() {
    auto t_pos = tools::common::profile::prof[algo_type]["t_pos"]->tic_token();
    print_status_full();
    tools::common::profile::print_profiling();
}

void class_itebd::single_TEBD_step() {
    /*!
     * \fn single_TEBD_step(class_superblock &state)
     * \brief infinite Time evolving block decimation.
     */
    for(auto &U : unitary_time_evolving_operators) {
        Eigen::Tensor<Scalar, 3> twosite_tensor = tools::infinite::opt::time_evolve_state(*tensors.state, U);
        tensors.merge_twosite_tensor(twosite_tensor, status.chi_lim);
        if(&U != &unitary_time_evolving_operators.back()) { tensors.state->swap_AB(); }
    }
    status.phys_time += std::abs(status.delta_t);
    tensors.clear_measurements();
    status.wall_time = tools::common::profile::t_tot->get_measured_time();
    status.algo_time = tools::common::profile::prof[algo_type]["t_sim"]->get_measured_time();
}

void class_itebd::check_convergence() {
    auto t_con = tools::common::profile::prof[algo_type]["t_con"]->tic_token();
    check_convergence_entg_entropy();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    update_bond_dimension_limit();
    check_convergence_time_step();
    if(status.entanglement_has_converged and status.variance_ham_has_converged and status.variance_mom_has_converged and status.chi_lim_has_reached_chi_max and
       status.time_step_has_converged) {
        status.algorithm_has_converged = true;
    }
}

void class_itebd::check_convergence_time_step() {
    if(std::abs(status.delta_t) <= settings::itebd::time_step_min) {
        status.time_step_has_converged = true;
    } else if(status.chi_lim_has_reached_chi_max and status.entanglement_has_converged) {
        status.delta_t                  = std::max(settings::itebd::time_step_min, std::abs(status.delta_t) * 0.5);
        unitary_time_evolving_operators = qm::timeEvolution::get_twosite_time_evolution_operators(status.delta_t, settings::itebd::suzuki_order, h_evn, h_odd);
        //        state->H->update_evolution_step_size(-status.delta_t, settings::itebd::suzuki_order);
        clear_convergence_status();
    }
}


bool class_itebd::cfg_algorithm_is_on() { return settings::itebd::on; }
long class_itebd::cfg_chi_lim_max() { return settings::itebd::chi_lim_max; }
size_t class_itebd::cfg_print_freq() { return settings::itebd::print_freq; }
bool   class_itebd::cfg_chi_lim_grow() { return settings::itebd::chi_lim_grow; }
long   class_itebd::cfg_chi_lim_init() { return settings::itebd::chi_lim_init; }
