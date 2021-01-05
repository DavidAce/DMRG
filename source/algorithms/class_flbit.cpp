//
// Created by david on 2018-01-31.
//

#include "class_flbit.h"
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <physics/nmspc_quantum_mechanics.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/io.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <tools/finite/opt.h>
#include <tools/finite/print.h>
#include <unsupported/Eigen/CXX11/Tensor>

#include <iostream>
#include <math/num.h>
class_flbit::class_flbit(std::shared_ptr<h5pp::File> h5pp_file_) : class_algorithm_finite(std::move(h5pp_file_), AlgorithmType::fLBIT) {
    tools::log->trace("Constructing class_flbit");
    tensors.state->set_name("state_real");
}

void class_flbit::resume() {
    // Resume can imply many things
    // 1) Resume a simulation which terminated prematurely
    // 2) Resume a previously successful simulation. This may be desireable if the config
    //    wants something that is not present in the file.
    //      a) A certain number of states
    //      b) A state inside of a particular energy window
    //      c) The ground or "roof" states
    // To guide the behavior, we check the setting ResumePolicy.

    auto state_prefix = tools::common::io::h5resume::find_resumable_state(*h5pp_file, algo_type);
    if(state_prefix.empty()) throw std::runtime_error("Could not resume: no valid state candidates found for resume");
    tools::log->info("Resuming state [{}]", state_prefix);
    tools::finite::io::h5resume::load_simulation(*h5pp_file, state_prefix, tensors, status);
    clear_convergence_status();
    // Our first task is to decide on a state name for the newly loaded state
    // The simplest is to inferr it from the state prefix itself
    auto name = tools::common::io::h5resume::extract_state_name(state_prefix);

    // Initialize a custom task list
    std::list<flbit_task> task_list;

    if(not status.algorithm_has_finished) {
        // This could be a checkpoint state
        // Simply "continue" the algorithm until convergence
        if(name.find("lbit") != std::string::npos)
            task_list.emplace_back(flbit_task::TIME_EVOLVE);
        else
            throw std::runtime_error(fmt::format("Unrecognized state name for flbit: [{}]", name));
        task_list.emplace_back(flbit_task::POST_DEFAULT);
        run_task_list(task_list);
    }
    // If we reached this point the current state has finished for one reason or another.
    // TODO: We may still have some more things to do, e.g. the config may be asking for more states
}

void class_flbit::run_task_list(std::list<flbit_task> &task_list) {
    while(not task_list.empty()) {
        auto task = task_list.front();
        switch(task) {
            case flbit_task::INIT_RANDOMIZE_MODEL: randomize_model(); break;
            case flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE: randomize_state(ResetReason::INIT, StateInit::RANDOM_PRODUCT_STATE); break;
            case flbit_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE: randomize_state(ResetReason::INIT, StateInit::RANDOM_ENTANGLED_STATE); break;
            case flbit_task::INIT_BOND_DIM_LIMITS: init_bond_dimension_limits(); break;
            case flbit_task::INIT_WRITE_MODEL: write_to_file(StorageReason::MODEL); break;
            case flbit_task::INIT_CLEAR_STATUS: status.clear(); break;
            case flbit_task::INIT_DEFAULT: run_preprocessing(); break;
            case flbit_task::TIME_EVOLVE:
                tensors.state->set_name("state_real");
                run_algorithm();
                break;
            case flbit_task::POST_WRITE_RESULT: write_to_file(StorageReason::FINISHED); break;
            case flbit_task::POST_PRINT_RESULT: print_status_full(); break;
            case flbit_task::POST_PRINT_PROFILING: tools::common::profile::print_profiling(algo_type); break;
            case flbit_task::POST_DEFAULT: run_postprocessing(); break;
            case flbit_task::PROF_RESET: tools::common::profile::reset_profiling(algo_type); break;
        }
        task_list.pop_front();
    }
}

void class_flbit::run_default_task_list() {
    std::list<flbit_task> default_task_list = {
        flbit_task::INIT_DEFAULT,
        flbit_task::TIME_EVOLVE,
        flbit_task::POST_DEFAULT,
    };

    run_task_list(default_task_list);
    if(not default_task_list.empty()) {
        for(auto &task : default_task_list) tools::log->critical("Unfinished task: {}", enum2str(task));
        throw std::runtime_error("Simulation ended with unfinished tasks");
    }
}

void class_flbit::run_preprocessing() {
    tools::log->info("Running {} preprocessing", algo_name);
    tools::common::profile::prof[algo_type]["t_pre"]->tic();
    status.clear();
    randomize_model(); // First use of random!
    init_bond_dimension_limits();

    // Create a state in the l-bit basis
    randomize_state(ResetReason::INIT, settings::strategy::initial_state);
    while(not tensors.state->position_is_any_edge()) tensors.move_center_point(status.chi_lim);
    if(not tensors.state->position_is_any_edge()){tools::finite::print::dimensions(tensors); throw std::logic_error("Put the state on an edge!");}

    tools::finite::print::model(*tensors.model);
    status.delta_t = std::complex<double>(settings::flbit::time_step_init_real,settings::flbit::time_step_init_imag);
    if(settings::model::model_size <= 10){
        // Create a copy of the state as a full state vector for ED comparison
        auto list_Lsite = num::range<size_t>(0,settings::model::model_size,1);
        Upsi_ed = tools::finite::measure::mps_wavefn(*tensors.state);
        tools::log->info("<Ψ_ed|Ψ_ed>   : {:.16f}",Textra::TensorVectorMap(Upsi_ed).norm());
        ham_gates_Lsite.emplace_back(qm::Gate(tensors.model->get_multisite_ham(list_Lsite,{1,2,3}), list_Lsite));
        time_gates_Lsite = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_Lsite);
    }

    // Create the hamiltonian gates with n-site terms
    auto                    list_1site = num::range<size_t>(0,settings::model::model_size-0,1);
    auto                    list_2site = num::range<size_t>(0,settings::model::model_size-1,1);
    auto                    list_3site = num::range<size_t>(0,settings::model::model_size-2,1);
    for(auto pos : list_1site) ham_gates_1body.emplace_back(qm::Gate(tensors.model->get_multisite_ham({pos}, {1}), {pos}));
    for(auto pos : list_2site) ham_gates_2body.emplace_back(qm::Gate(tensors.model->get_multisite_ham({pos,pos+1}, {2}), {pos,pos+1}));
    for(auto pos : list_3site) ham_gates_3body.emplace_back(qm::Gate(tensors.model->get_multisite_ham({pos,pos+1,pos+2},{3}), {pos,pos+1,pos+2}));

    for(auto &&ham : ham_gates_1body) std::cout << "ham 1body:\n" << ham.op << std::endl;
    for(auto &&ham : ham_gates_2body) std::cout << "ham 2body:\n" << ham.op << std::endl;
    for(auto &&ham : ham_gates_3body) std::cout << "ham 3body:\n" << ham.op << std::endl;
    for(auto &&ham : ham_gates_1body) if(Textra::TensorMatrixMap(ham.op).isZero()) throw std::runtime_error("Ham1 is all zeros");
    for(auto &&ham : ham_gates_2body) if(Textra::TensorMatrixMap(ham.op).isZero()) throw std::runtime_error("Ham2 is all zeros");
    for(auto &&ham : ham_gates_3body) if(Textra::TensorMatrixMap(ham.op).isZero()) throw std::runtime_error("Ham3 is all zeros");

    // Get the time evolution operators
    time_gates_1site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_1body);
    time_gates_2site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_2body);
    time_gates_3site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_3body);

    unitary_gates_2site_layer0 = qm::lbit::get_unitary_2gate_layer(settings::model::model_size,0.1);
    unitary_gates_2site_layer1 = qm::lbit::get_unitary_2gate_layer(settings::model::model_size,0.1);
    unitary_gates_2site_layer2 = qm::lbit::get_unitary_2gate_layer(settings::model::model_size,0.1);
    unitary_gates_2site_layer3 = qm::lbit::get_unitary_2gate_layer(settings::model::model_size,0.1);


    if(not tensors.state->position_is_any_edge()) throw std::logic_error("Put the state on an edge!");

    // Generate the corresponding state in lbit basis
    transform_to_lbit_basis();

    tools::common::profile::prof[algo_type]["t_pre"]->toc();
    tools::log->info("Finished {} preprocessing", algo_name);
}

void class_flbit::run_algorithm() {
    if(tensors.state->get_name().empty()) tensors.state->set_name("state_real");
    if(not state_lbit) throw std::logic_error("state_lbit == nullptr: Set the state in lbit basis before running the flbit algorithm");
    tools::log->info("Starting {} algorithm with model [{}] for state [{}]", algo_name, enum2str(settings::model::model_type),  tensors.state->get_name());
    tools::common::profile::prof[algo_type]["t_sim"]->tic();
    if(not tensors.state->position_is_any_edge()) throw std::logic_error("Put the state on an edge!");
    while(true) {
        single_flbit_step();
        check_convergence();
        print_status_update();
        print_profiling_lap();
        write_to_file();

        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);

        // It's important not to perform the last move, so we break now: that last state would not get optimized
        if(stop_reason != StopReason::NONE) break;
        update_time_step();
//        update_bond_dimension_limit(); // Will update bond dimension if the state precision is being limited by bond dimension
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", algo_name, tensors.state->get_name(), enum2str(stop_reason));
    status.algorithm_has_finished = true;
    tools::common::profile::prof[algo_type]["t_sim"]->toc();
}

void class_flbit::single_flbit_step() {
    /*!
     * \fn void single_DMRG_step(std::string ritz)
     */
    tools::log->trace("Starting single flbit time-step");
    tensors.activate_sites({0,1});
    tensors.clear_measurements();
    tensors.clear_cache();
    state_lbit->clear_cache();
    state_lbit->clear_measurements();

    // Time evolve here
    tools::common::profile::prof[algo_type]["t_evo"]->tic();
    tools::finite::mps::apply_gates(*state_lbit, time_gates_1site, false, status.chi_lim);
    tools::finite::mps::apply_gates(*state_lbit, time_gates_2site, false, status.chi_lim);
    tools::finite::mps::apply_gates(*state_lbit, time_gates_3site, false, status.chi_lim);
    tools::common::profile::prof[algo_type]["t_evo"]->toc();

    transform_to_real_basis();

    if(settings::model::model_size <= 10 and settings::debug){
        Eigen::Tensor<Scalar,1> Upsi_mps = tools::finite::measure::mps_wavefn(*tensors.state);
        Eigen::Tensor<Scalar,1> Upsi_tmp = time_gates_Lsite[0].op.contract(Upsi_ed, Textra::idx({1},{0}));
        Upsi_ed = Upsi_tmp * std::exp(std::arg(Upsi_tmp(0)) * Scalar(0,-1));
        Upsi_mps = Upsi_mps * std::exp(std::arg(Upsi_mps(0)) * Scalar(0,-1));
        Eigen::Tensor<Scalar,0> overlap = Upsi_ed.conjugate().contract(Upsi_mps, Textra::idx({0},{0}));
        tools::log->info("<UΨ_tmp|UΨ_tmp> : {:.16f}",Textra::TensorVectorMap(Upsi_ed).norm());
        tools::log->info("<UΨ_ed|UΨ_ed>   : {:.16f}",Textra::TensorVectorMap(Upsi_ed).norm());
        tools::log->info("<UΨ_mps|UΨ_mps> : {:.16f}",Textra::TensorVectorMap(Upsi_mps).norm());
        tools::log->info("<UΨ_ed|UΨ_mps>  : {:.16f}{:+.16f}i",std::real(overlap(0)),std::imag(overlap(0)));
    }
    tensors.clear_measurements();
    tensors.clear_cache();
    tensors.rebuild_edges_ene();
    status.iter     += 1;
    status.step     += settings::model::model_size;
    status.phys_time += std::abs(status.delta_t);
    status.wall_time = tools::common::profile::t_tot->get_measured_time();
    status.algo_time = tools::common::profile::prof[algo_type]["t_sim"]->get_measured_time();

}

void class_flbit::update_time_step(){
    if(std::abs(status.delta_t) * settings::flbit::time_step_growth_factor >= settings::flbit::time_step_max_size) return;
    static double time_last_update = 0;
    if(std::abs(status.delta_t) == 0.0) status.delta_t = std::complex<double>(settings::flbit::time_step_init_real, settings::flbit::time_step_init_imag);
    double time_until_update = time_last_update + std::abs(status.delta_t) * settings::flbit::time_step_growth_factor*10 - status.phys_time;
    tools::log->trace("time_last_update : {}",time_last_update);
    tools::log->trace("time_until_update: {}",time_until_update);
    tools::log->trace("phys_time        : {}",status.phys_time);
    tools::log->trace("delta_t          : {:.16f}{:+.16f}i",status.delta_t.real(),status.delta_t.imag());
    if(time_until_update <= 1e-10){
        time_last_update = status.phys_time;
        status.delta_t *= settings::flbit::time_step_growth_factor;
        time_gates_1site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_1body);
        time_gates_2site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_2body);
        time_gates_3site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_3body);
        if(settings::model::model_size <= 10 and settings::debug){
            time_gates_Lsite = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_Lsite);
        }
    }
}


void class_flbit::check_convergence() {
    tools::common::profile::prof[algo_type]["t_con"]->tic();
    if(tensors.position_is_any_edge()) {
        check_convergence_entg_entropy();
    }

    status.algorithm_has_saturated = status.entanglement_saturated_for >= min_saturation_iters;
    status.algorithm_has_converged = status.entanglement_has_converged;
    status.algorithm_has_succeeded = status.algorithm_has_saturated and status.algorithm_has_converged;
    status.algorithm_has_got_stuck = status.algorithm_has_saturated and not status.algorithm_has_converged;

    if(tensors.state->position_is_any_edge()) status.algorithm_has_stuck_for = status.algorithm_has_got_stuck ? status.algorithm_has_stuck_for + 1 : 0;
    status.algorithm_has_to_stop = status.algorithm_has_stuck_for >= max_stuck_iters;

    if(tensors.state->position_is_any_edge()) {
        tools::log->debug("Simulation report: converged {} | saturated {} | succeeded {} | stuck {} for {} iters | has to stop {}",
                          status.algorithm_has_converged, status.algorithm_has_saturated, status.algorithm_has_succeeded, status.algorithm_has_got_stuck,
                          status.algorithm_has_stuck_for, status.algorithm_has_to_stop);
    }
    stop_reason = StopReason::NONE;
    if(tensors.position_is_any_edge() and status.iter > settings::flbit::min_iters) {
        if(status.iter >= settings::flbit::max_iters) stop_reason = StopReason::MAX_ITERS;
        if(status.phys_time >= settings::flbit::time_limit) stop_reason = StopReason::SUCCEEDED;
//        if(status.algorithm_has_succeeded) stop_reason = StopReason::SUCCEEDED;
//        if(status.algorithm_has_to_stop) stop_reason = StopReason::SATURATED;
        if(status.num_resets > settings::strategy::max_resets) stop_reason = StopReason::MAX_RESET;
    }

    tools::common::profile::prof[algo_type]["t_con"]->toc();
}

void class_flbit::transform_to_real_basis(){
    tools::common::profile::prof[algo_type]["t_map"]->tic();
    tensors.state = std::make_unique<class_state_finite>(*state_lbit);
    tensors.state->set_name("state_real");
    tools::log->info("Transforming {} to {}", state_lbit->get_name(), tensors.state->get_name());
    tensors.clear_measurements();
    tensors.clear_cache();
    tools::finite::mps::apply_gates(*tensors.state,unitary_gates_2site_layer0, false, status.chi_lim);
    tools::finite::mps::apply_gates(*tensors.state,unitary_gates_2site_layer1, false, status.chi_lim);
    tools::finite::mps::apply_gates(*tensors.state,unitary_gates_2site_layer2, false, status.chi_lim);
    tools::finite::mps::apply_gates(*tensors.state,unitary_gates_2site_layer3, false, status.chi_lim);
    auto has_normalized = tools::finite::mps::normalize_state(*tensors.state, status.chi_lim, settings::precision::svd_threshold, NormPolicy::IFNEEDED);
    if constexpr(settings::debug)
        if(has_normalized and tools::log->level() == spdlog::level::trace){
            tools::common::profile::prof[algo_type]["t_dbg"]->tic();
            Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "  [", "]");
            tools::log->trace("After normalization");
            for(auto &&mps : tensors.state->mps_sites)
                std::cout << "M(" << mps->get_position() << ") dims [" << mps->spin_dim() << "," << mps->get_chiL() << "," << mps->get_chiR() << "]:\n"
                          << Textra::TensorMatrixMap(mps->get_M_bare(), mps->spin_dim(), mps->get_chiL() * mps->get_chiR()).format(CleanFmt) << std::endl;
            tools::common::profile::prof[algo_type]["t_dbg"]->toc();
        }
    tools::common::profile::prof[algo_type]["t_map"]->toc();
}

void class_flbit::transform_to_lbit_basis(){
    state_lbit = std::make_unique<class_state_finite>(*tensors.state);
    state_lbit->set_name("state_lbit");
    tools::log->info("Transforming {} to {}", tensors.state->get_name(),state_lbit->get_name());
    state_lbit->clear_cache();
    state_lbit->clear_measurements();
    tensors.clear_measurements();
    tensors.clear_cache();
    tools::common::profile::prof[algo_type]["t_map"]->tic();
    tools::finite::mps::apply_gates(*state_lbit,unitary_gates_2site_layer0, true, status.chi_lim);
    tools::finite::mps::apply_gates(*state_lbit,unitary_gates_2site_layer1, true, status.chi_lim);
    tools::finite::mps::apply_gates(*state_lbit,unitary_gates_2site_layer2, true, status.chi_lim);
    tools::finite::mps::apply_gates(*state_lbit,unitary_gates_2site_layer3, true, status.chi_lim);
    tools::common::profile::prof[algo_type]["t_map"]->toc();
    tools::common::profile::prof[AlgorithmType::ANY]["t_map_norm"]->tic();
    auto has_normalized = tools::finite::mps::normalize_state(*state_lbit, status.chi_lim, settings::precision::svd_threshold, NormPolicy::IFNEEDED);
    tools::common::profile::prof[AlgorithmType::ANY]["t_map_norm"]->toc();
    if constexpr(settings::debug)
        if(has_normalized and tools::log->level() == spdlog::level::trace){
            Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "  [", "]");
            tools::log->trace("After normalization");
            for(auto &&mps : state_lbit->mps_sites)
                std::cout << "M(" << mps->get_position() << ") dims [" << mps->spin_dim() << "," << mps->get_chiL() << "," << mps->get_chiR() << "]:\n"
                          << Textra::TensorMatrixMap(mps->get_M_bare(), mps->spin_dim(), mps->get_chiL() * mps->get_chiR()).format(CleanFmt) << std::endl;
        }
}

void class_flbit::write_to_file(StorageReason storage_reason, std::optional<CopyPolicy> copy_file) {
    tools::common::profile::prof[AlgorithmType::ANY]["t_write_h5pp"]->tic();
    tensors.clear_cache();
    tensors.clear_measurements();
    class_algorithm_finite::write_to_file(storage_reason, *tensors.state, copy_file);
    tensors.clear_cache();
    tensors.clear_measurements();
    class_algorithm_finite::write_to_file(storage_reason, *state_lbit, copy_file);
    tensors.clear_cache();
    tensors.clear_measurements();
    tools::common::profile::prof[AlgorithmType::ANY]["t_write_h5pp"]->toc();

}


bool   class_flbit::cfg_algorithm_is_on() { return settings::flbit::on; }
long   class_flbit::cfg_chi_lim_max() { return settings::flbit::chi_lim_max; }
size_t class_flbit::cfg_print_freq() { return settings::flbit::print_freq; }
bool   class_flbit::cfg_chi_lim_grow() { return settings::flbit::chi_lim_grow; }
long   class_flbit::cfg_chi_lim_init() { return settings::flbit::chi_lim_init; }
bool   class_flbit::cfg_store_wave_function() { return settings::flbit::store_wavefn; }



void class_flbit::print_status_update() {
    if(num::mod(status.iter, cfg_print_freq()) != 0) return;
    if(cfg_print_freq() == 0) return;
    tools::common::profile::prof[algo_type]["t_out"]->tic();
    std::string report;
    report += fmt::format("{:<} ", tensors.state->get_name());
    report += fmt::format("iter:{:<4} ", status.iter);
    report += fmt::format("step:{:<5} ", status.step);
    report += fmt::format("L:{} ", tensors.get_length());
    if(tensors.active_sites.empty()) report += fmt::format("l:{:<2} ", tensors.get_position());
    else if(tensors.state->get_direction() > 0)
        report += fmt::format("l:[{:>2}-{:<2}] ", tensors.active_sites.front(), tensors.active_sites.back());
    else if(tensors.state->get_direction() < 0)
        report += fmt::format("l:[{:>2}-{:<2}] ", tensors.active_sites.back(), tensors.active_sites.front());
    report += fmt::format("E/L:{:<20.16f} ", tools::finite::measure::energy_per_site(tensors));
    if(algo_type == AlgorithmType::xDMRG) { report += fmt::format("ε:{:<6.4f} ", status.energy_dens); }
    report += fmt::format("Sₑ(L/2):{:<10.8f} ", tools::finite::measure::entanglement_entropy_midchain(*tensors.state));
//    report += fmt::format("log₁₀σ²E:{:<10.6f} [{:<10.6f}] ", std::log10(tools::finite::measure::energy_variance(tensors)),
//                          std::log10(status.energy_variance_lowest));
    report +=
        fmt::format("χ:{:<3}|{:<3}|{:<3} ", cfg_chi_lim_max(), status.chi_lim, tools::finite::measure::bond_dimension_midchain(*tensors.state));

    report += fmt::format("log₁₀trnc:{:<8.4f} ", std::log10(tensors.state->get_truncation_error_midchain()));
    report += fmt::format("stk:{:<1} ", status.algorithm_has_stuck_for);
    report += fmt::format("sat:[σ² {:<1} Sₑ {:<1}] ", status.variance_mpo_saturated_for, status.entanglement_saturated_for);
    report += fmt::format("wtime:{:<} ",fmt::format("{:>6.2f}s",tools::common::profile::t_tot->get_measured_time()));
    report += fmt::format("ptime:{:<} ",fmt::format("{:>6.2f}s",status.phys_time));
    report += fmt::format("mem[rss {:<.1f}|peak {:<.1f}|vm {:<.1f}]MB ", tools::common::profile::mem_rss_in_mb(), tools::common::profile::mem_hwm_in_mb(),
                          tools::common::profile::mem_vm_in_mb());
    tools::log->info(report);
    tools::common::profile::prof[algo_type]["t_out"]->toc();
}