#include "flbit.h"
#include <config/settings.h>
#include <debug/exceptions.h>
#include <general/iter.h>
#include <h5pp/h5pp.h>
#include <io/fmt.h>
#include <math/num.h>
#include <math/tenx.h>
#include <qm/lbit.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateFinite.h>
#include <tid/tid.h>
#include <tools/common/h5.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/h5.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <tools/finite/opt.h>
#include <tools/finite/print.h>
#include <unsupported/Eigen/CXX11/Tensor>

flbit::flbit(std::shared_ptr<h5pp::File> h5pp_file_) : AlgorithmFinite(std::move(h5pp_file_), AlgorithmType::fLBIT) {
    tools::log->trace("Constructing class_flbit");
    tensors.state->set_name("state_real");
}

void flbit::resume() {
    // Resume can imply many things
    // 1) Resume a simulation which terminated prematurely
    // 2) Resume a previously successful simulation. This may be desireable if the config
    //    wants something that is not present in the file.
    //      a) A certain number of states
    //      b) A state inside of a particular energy window
    //      c) The ground or "roof" states
    // To guide the behavior, we check the setting ResumePolicy.

    auto state_prefix = tools::common::h5::resume::find_resumable_state(*h5pp_file, status.algo_type, "state_real");
    if(state_prefix.empty()) throw except::state_error("Could not resume: no valid state candidates found for resume");
    tools::log->info("Resuming state [{}]", state_prefix);
    tools::finite::h5::load::simulation(*h5pp_file, state_prefix, tensors, status, status.algo_type);

    // Load the unitaries
    unitary_gates_2site_layers.clear();
    std::string grouppath = "/fLBIT/model/unitary_gates";
    if(not h5pp_file->linkExists(grouppath)) throw std::runtime_error(fmt::format("Missing link: {}", grouppath));
    auto num_layers = h5pp_file->readAttribute<size_t>("num_layers", grouppath);
    if(num_layers != settings::model::lbit::u_layer)
        throw std::runtime_error(fmt::format("Mismatch in number of layers: file {} != cfg {}", num_layers != settings::model::lbit::u_layer));
    unitary_gates_2site_layers.resize(num_layers);
    for(auto &&[idx_layer, layer] : iter::enumerate(unitary_gates_2site_layers)) {
        std::string layerpath = fmt::format("{}/layer_{}", grouppath, idx_layer);
        auto        num_gates = h5pp_file->readAttribute<size_t>("num_gates", layerpath);
        layer.resize(num_gates);
        for(auto &&[idx_gate, u] : iter::enumerate(layer)) {
            std::string gatepath = fmt::format("{}/u_{}", layerpath, idx_gate);
            auto        op       = h5pp_file->readDataset<Eigen::Tensor<Scalar, 2>>(gatepath);
            auto        pos      = h5pp_file->readAttribute<std::vector<size_t>>("pos", gatepath);
            auto        dim      = h5pp_file->readAttribute<std::vector<long>>("dim", gatepath);
            u                    = qm::Gate(op, pos, dim);
        }
    }
    clear_convergence_status();
    tensors.move_center_point_to_edge(status.chi_lim);
    tensors.rebuild_edges();

    // Our first task is to decide on a state name for the newly loaded state
    // The simplest is to inferr it from the state prefix itself
    auto name = tools::common::h5::resume::extract_state_name(state_prefix);

    // Initialize a custom task list
    std::deque<flbit_task> task_list;

    if(status.algorithm_has_finished) {
        task_list.emplace_back(flbit_task::POST_PRINT_RESULT);
    } else {
        // This could be a savepoint state
        // Simply "continue" the algorithm until convergence
        if(name.find("real") != std::string::npos) {
            task_list.emplace_back(flbit_task::INIT_TIME);
            task_list.emplace_back(flbit_task::INIT_GATES);
            task_list.emplace_back(flbit_task::TRANSFORM_TO_LBIT);
            task_list.emplace_back(flbit_task::TIME_EVOLVE);
        } else
            throw std::runtime_error(fmt::format("Unrecognized state name for flbit: [{}]", name));
        task_list.emplace_back(flbit_task::POST_DEFAULT);
    }
    run_task_list(task_list);
    // If we reached this point the current state has finished for one reason or another.
    // TODO: We may still have some more things to do, e.g. the config may be asking for more states
}

void flbit::run_task_list(std::deque<flbit_task> &task_list) {
    while(not task_list.empty()) {
        auto task = task_list.front();
        switch(task) {
            case flbit_task::INIT_RANDOMIZE_MODEL: randomize_model(); break;
            case flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE: randomize_state(ResetReason::INIT, StateInit::RANDOM_PRODUCT_STATE); break;
            case flbit_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE: randomize_state(ResetReason::INIT, StateInit::RANDOM_ENTANGLED_STATE); break;
            case flbit_task::INIT_BOND_DIM_LIMITS: init_bond_dimension_limits(); break;
            case flbit_task::INIT_WRITE_MODEL: write_to_file(StorageReason::MODEL); break;
            case flbit_task::INIT_CLEAR_STATUS: status.clear(); break;
            case flbit_task::INIT_CLEAR_CONVERGENCE: clear_convergence_status(); break;
            case flbit_task::INIT_DEFAULT: run_preprocessing(); break;
            case flbit_task::INIT_TIME: {
                create_time_points();
                update_time_step();
                break;
            }
            case flbit_task::INIT_GATES: {
                create_hamiltonian_gates();
                create_time_evolution_gates();
                create_lbit_transform_gates();
                break;
            }
            case flbit_task::TIME_EVOLVE:
                tensors.state->set_name("state_real");
                run_algorithm();
                break;
            case flbit_task::TRANSFORM_TO_LBIT: transform_to_lbit_basis(); break;
            case flbit_task::TRANSFORM_TO_REAL: transform_to_real_basis(); break;
            case flbit_task::POST_WRITE_RESULT: write_to_file(StorageReason::FINISHED); break;
            case flbit_task::POST_PRINT_RESULT: print_status_full(); break;
            case flbit_task::POST_PRINT_PROFILING: tools::common::profile::print_profiling(status); break;
            case flbit_task::POST_DEFAULT: run_postprocessing(); break;
            case flbit_task::PROF_RESET: tid::reset("fLBIT"); break;
        }
        task_list.pop_front();
    }
}

void flbit::run_default_task_list() {
    std::deque<flbit_task> default_task_list = {
        flbit_task::INIT_DEFAULT,
        flbit_task::TIME_EVOLVE,
        flbit_task::POST_DEFAULT,
    };

    run_task_list(default_task_list);
    if(not default_task_list.empty()) {
        for(auto &task : default_task_list) tools::log->critical("Unfinished task: {}", enum2sv(task));
        throw std::runtime_error("Simulation ended with unfinished tasks");
    }
}

void flbit::run_preprocessing() {
    tools::log->info("Running {} preprocessing", status.algo_type_sv());
    status.clear();
    randomize_model(); // First use of random!
    init_bond_dimension_limits();

    // Create a state in the l-bit basis
    randomize_state(ResetReason::INIT, settings::strategy::initial_state);
    tensors.move_center_point_to_edge(status.chi_lim);
    tools::finite::print::model(*tensors.model);
    create_time_points();
    update_time_step();
    if(settings::model::model_size <= 10) {
        // Create a copy of the state as a full state vector for ED comparison
        auto list_Lsite = num::range<size_t>(0, settings::model::model_size, 1);
        Upsi_ed         = tools::finite::measure::mps_wavefn(*tensors.state);
        tools::log->info("<Ψ_ed|Ψ_ed>   : {:.16f}", tenx::VectorMap(Upsi_ed).norm());
        ham_gates_Lsite.emplace_back(qm::Gate(tensors.model->get_multisite_ham(list_Lsite, {1, 2, 3}), list_Lsite, tensors.state->get_spin_dims(list_Lsite)));
        time_gates_Lsite = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_Lsite);
    }

    create_hamiltonian_gates();
    create_time_evolution_gates();
    create_lbit_transform_gates();

    if(not tensors.position_is_inward_edge()) throw std::logic_error("Put the state on an edge!");

    // Generate the corresponding state in lbit basis
    transform_to_lbit_basis();

    write_to_file(StorageReason::MODEL, CopyPolicy::TRY);
    tools::log->info("Finished {} preprocessing", status.algo_type_sv());
}

void flbit::run_algorithm() {
    if(tensors.state->get_name().empty()) tensors.state->set_name("state_real");
    if(not state_lbit) transform_to_lbit_basis();
    if(time_points.empty()) create_time_points();
    if(std::abs(status.delta_t) == 0) update_time_step();
    tools::log->info("Starting {} algorithm with model [{}] for state [{}]", status.algo_type_sv(), enum2sv(settings::model::model_type),
                     tensors.state->get_name());
    auto t_run = tid::tic_scope("run");
    if(not tensors.position_is_inward_edge()) throw std::logic_error("Put the state on an edge!");
    while(true) {
        single_flbit_step();
        check_convergence();
        write_to_file(StorageReason::SAVEPOINT);
        write_to_file(StorageReason::CHECKPOINT);
        print_status_update();
        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);

        // It's important not to perform the last move, so we break now: that last state would not get optimized
        if(status.algo_stop != AlgorithmStop::NONE) break;
        update_time_step();
        status.wall_time = tid::get_unscoped("t_tot").get_time();
        status.algo_time = t_run->get_time();
        t_run->start_lap();
        tools::common::profile::print_profiling();
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", status.algo_type_sv(), tensors.state->get_name(), status.algo_stop_sv());
    status.algorithm_has_finished = true;
}

void flbit::single_flbit_step() {
    /*!
     * \fn void single_DMRG_step(std::string ritz)
     */
    tools::log->trace("Starting single flbit time-step");
    if(not state_lbit) throw std::logic_error("state_lbit == nullptr: Set the state in lbit basis before running an flbit step");

    // Time evolve here
    auto t_step = tid::tic_scope("step");
    auto t_evo  = tid::tic_scope("time_evo");
    tools::finite::mps::apply_gates(*state_lbit, time_gates_1site, false, status.chi_lim);
    tools::finite::mps::apply_gates(*state_lbit, time_gates_2site, false, status.chi_lim);
    tools::finite::mps::apply_gates(*state_lbit, time_gates_3site, false, status.chi_lim);
    t_evo.toc();
    transform_to_real_basis();

    if(settings::model::model_size <= 10 and settings::debug) {
        Eigen::Tensor<Scalar, 1> Upsi_mps = tools::finite::measure::mps_wavefn(*tensors.state);
        Eigen::Tensor<Scalar, 1> Upsi_tmp = time_gates_Lsite[0].op.contract(Upsi_ed, tenx::idx({1}, {0}));
        Upsi_ed                           = Upsi_tmp * std::exp(std::arg(Upsi_tmp(0)) * Scalar(0, -1));
        Upsi_mps                          = Upsi_mps * std::exp(std::arg(Upsi_mps(0)) * Scalar(0, -1));
        Eigen::Tensor<Scalar, 0> overlap  = Upsi_ed.conjugate().contract(Upsi_mps, tenx::idx({0}, {0}));
        tools::log->info("<UΨ_tmp|UΨ_tmp> : {:.16f}", tenx::VectorMap(Upsi_ed).norm());
        tools::log->info("<UΨ_ed|UΨ_ed>   : {:.16f}", tenx::VectorMap(Upsi_ed).norm());
        tools::log->info("<UΨ_mps|UΨ_mps> : {:.16f}", tenx::VectorMap(Upsi_mps).norm());
        tools::log->info("<UΨ_ed|UΨ_mps>  : {:.16f}{:+.16f}i", std::real(overlap(0)), std::imag(overlap(0)));
    }
    tensors.clear_measurements();
    tensors.clear_cache();
    tensors.rebuild_edges_ene();
    status.iter += 1;
    status.step += settings::model::model_size;
    status.position  = tensors.get_position<long>();
    status.direction = tensors.state->get_direction();
    status.phys_time += std::abs(status.delta_t);
}

void flbit::update_time_step() {
    tools::log->trace("Updating time step");
    auto t_updtstep = tid::tic_scope("update_time_step");
    if(time_points.empty()) create_time_points();
    auto time_point_idx0 = std::clamp(status.iter + 0, 0ul, time_points.size() - 1);
    auto time_point_idx1 = std::clamp(status.iter + 1, 0ul, time_points.size() - 1);
    if(time_point_idx0 == time_point_idx1) {
        status.algo_stop = AlgorithmStop::SUCCEEDED;
        return;
    }
    if(time_point_idx0 > time_point_idx1) throw std::logic_error(fmt::format("Time order error: idx0 ({}) > idx1 ({})", time_point_idx0, time_point_idx1));
    if(std::abs(std::abs(time_points[time_point_idx0]) - status.phys_time) > 1e-10)
        tools::log->warn("Physical time is currently {:.8e} | should be {:.8e}", status.phys_time, std::abs(time_points[time_point_idx0]));
    status.delta_t = time_points[time_point_idx1] - time_points[time_point_idx0];
    tools::log->trace("Time step iter {} = {:.8e}", status.iter, std::abs(status.delta_t));
    if(std::abs(status.delta_t) == 0) throw std::logic_error("Expected nonzero delta_t after time step update");
    time_gates_1site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_1body);
    time_gates_2site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_2body);
    time_gates_3site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_3body);
    if constexpr(settings::debug)
        if(settings::model::model_size <= 10) time_gates_Lsite = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_Lsite);
}

void flbit::check_convergence() {
    if(not tensors.position_is_inward_edge()) return;
    auto t_con = tid::tic_scope("conv");
    check_convergence_entg_entropy();
    if(status.entanglement_saturated_for > 0)
        status.algorithm_saturated_for++;
    else
        status.algorithm_saturated_for = 0;

    status.algorithm_converged_for = status.iter + 1 - std::min(settings::flbit::time_num_steps, status.iter + 1);
    status.algorithm_has_succeeded = status.algorithm_converged_for > 0;
    if(status.algorithm_saturated_for > min_saturation_iters and status.algorithm_converged_for == 0)
        status.algorithm_has_stuck_for++;
    else
        status.algorithm_has_stuck_for = 0;

    status.algorithm_has_to_stop = status.algorithm_has_stuck_for >= max_stuck_iters;

    tools::log->debug("Simulation report: converged {} | saturated {} | stuck {} | succeeded {} | has to stop {}", status.algorithm_converged_for,
                      status.algorithm_saturated_for, status.algorithm_has_stuck_for, status.algorithm_has_succeeded, status.algorithm_has_to_stop);
    status.algo_stop = AlgorithmStop::NONE;
    if(status.iter >= settings::flbit::min_iters) {
        if(status.iter >= settings::flbit::max_iters) status.algo_stop = AlgorithmStop::MAX_ITERS;
        if(status.iter >= settings::flbit::time_num_steps) status.algo_stop = AlgorithmStop::SUCCEEDED;
        if(status.algorithm_has_to_stop) status.algo_stop = AlgorithmStop::SATURATED;
        if(status.phys_time >= std::abs(std::complex<double>(settings::flbit::time_final_real, settings::flbit::time_final_imag)))
            status.algo_stop = AlgorithmStop::SUCCEEDED;
        if(status.num_resets > settings::strategy::max_resets) status.algo_stop = AlgorithmStop::MAX_RESET;
    }
}

void flbit::create_time_points() {
    tools::log->info("Creating time points");
    auto t_crt      = tid::tic_scope("create_time_points");
    auto time_start = std::complex<double>(settings::flbit::time_start_real, settings::flbit::time_start_imag);
    auto time_final = std::complex<double>(settings::flbit::time_final_real, settings::flbit::time_final_imag);
    auto time_diff  = time_start - time_final;
    // Check that there will be some time evolution
    if(std::abs(time_diff) == 0) throw std::logic_error("time_start - time_final == 0");
    // Check that the time limits are purely real or imaginary!
    bool time_is_real = std::abs(time_diff.real()) > 0;
    bool time_is_imag = std::abs(time_diff.imag()) > 0;
    if(time_is_real and time_is_imag)
        throw std::logic_error(fmt::format("time_start and time_final must both be either purely real or imaginary. Got:\n"
                                           "time_start = {:.8f}{:+.8f}\n"
                                           "time_final = {:.8f}{:+.8f}",
                                           time_start.real(), time_start.imag(), time_final.real(), time_final.imag()));
    time_points.reserve(settings::flbit::time_num_steps);
    if(time_is_real)
        for(const auto &t : num::LogSpaced(settings::flbit::time_num_steps, std::real(time_start), std::real(time_final))) {
            if(std::isinf(t) or std::isnan(t)) throw std::runtime_error(fmt::format("Invalid time point: {}", t));
            time_points.emplace_back(t);
        }
    else
        for(const auto &t : num::LogSpaced(settings::flbit::time_num_steps, std::imag(time_start), std::imag(time_final))) {
            if(std::isinf(t) or std::isnan(t)) throw std::runtime_error(fmt::format("Invalid time point: {}", t));
            time_points.emplace_back(std::complex<double>(0, t));
        }
    status.phys_time = std::abs(time_points[std::min(time_points.size() - 1, status.iter)]);
    tools::log->trace("Created time points: {}", time_points);
    tools::log->trace("Current physical time {:.8e}", status.phys_time);
}

void flbit::create_hamiltonian_gates() {
    tools::log->info("Creating Hamiltonian gates");
    auto t_hamgates = tid::tic_scope("hamgates");
    ham_gates_1body.clear();
    ham_gates_2body.clear();
    ham_gates_3body.clear();
    // Create the hamiltonian gates with n-site terms
    auto list_1site = num::range<size_t>(0, settings::model::model_size - 0, 1);
    auto list_2site = num::range<size_t>(0, settings::model::model_size - 1, 1);
    auto list_3site = num::range<size_t>(0, settings::model::model_size - 2, 1);
    for(auto pos : list_1site) {
        auto range = 1ul;                                  // Number of site indices
        auto sites = num::range<size_t>(pos, pos + range); // A list of site indices
        auto nbody = std::vector<size_t>{1};               // A list of included nbody interaction terms (1: on-site terms, 2: pairwise, and so on)
        auto spins = tensors.state->get_spin_dims(sites);  // A list of spin dimensions for each site (should all be 2 for two-level systems)
        tools::log->info("Generating {}-body hamiltonian on sites {}", nbody, sites);
        ham_gates_1body.emplace_back(qm::Gate(tensors.model->get_multisite_ham({pos}, nbody), sites, spins));
    }
    for(auto pos : list_2site) {
        auto range = settings::model::lbit::J2_span; // Max distance |i-j| to the furthest interacting site
        auto sites = num::range<size_t>(pos, std::clamp<size_t>(pos + range + 1, 0ul, settings::model::model_size)); // +1 to include last site
        auto nbody = std::vector<size_t>{2};
        auto spins = tensors.state->get_spin_dims(sites);
        tools::log->info("Generating {}-body hamiltonian on sites {}", nbody, sites);
        ham_gates_2body.emplace_back(qm::Gate(tensors.model->get_multisite_ham(sites, nbody), sites, spins));
    }
    for(auto pos : list_3site) {
        auto range = 2ul; // Distance to next-nearest neighbor when 3 sites interact
        auto sites = num::range<size_t>(pos, std::clamp<size_t>(pos + range + 1, 0ul, settings::model::model_size)); // +1 to include last site
        auto nbody = std::vector<size_t>{3};
        auto spins = tensors.state->get_spin_dims(sites);
        tools::log->info("Generating {}-body hamiltonian on sites {}", nbody, sites);
        ham_gates_3body.emplace_back(qm::Gate(tensors.model->get_multisite_ham(sites, nbody), sites, spins));
    }
    for(const auto &[idx, ham] : iter::enumerate(ham_gates_1body)) {
        if(tenx::MatrixMap(ham.op).isZero()) {
            throw std::runtime_error(fmt::format("ham1[{}] is all zeros", idx));
            tools::log->warn("Ham1 is all zeros");
        }
    }
    for(const auto &[idx, ham] : iter::enumerate(ham_gates_2body)) {
        if(tenx::MatrixMap(ham.op).isZero()) {
            throw std::runtime_error(fmt::format("ham2[{}] is all zeros", idx));
            tools::log->warn("Ham2 is all zeros");
        }
    }
    for(const auto &[idx, ham] : iter::enumerate(ham_gates_3body)) {
        if(tenx::MatrixMap(ham.op).isZero()) {
            throw std::runtime_error(fmt::format("ham3[{}] is all zeros", idx));
            tools::log->warn("Ham3 is all zeros");
        }
    }
    //    for(const auto &ham : ham_gates_2body)
    //        if(tenx::MatrixMap(ham.op).isZero()) tools::log->warn("Ham2 is all zeros");
    //    for(const auto &ham : ham_gates_3body)
    //        if(tenx::MatrixMap(ham.op).isZero()) tools::log->warn("Ham3 is all zeros");
}
void flbit::create_time_evolution_gates() {
    // Create the time evolution operators
    if(time_points.empty()) create_time_points();
    if(ham_gates_1body.empty() or ham_gates_2body.empty() or ham_gates_3body.empty()) create_hamiltonian_gates();
    if(std::abs(status.delta_t) == 0) update_time_step();
    tools::log->info("Creating time evolution gates");
    time_gates_1site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_1body);
    time_gates_2site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_2body);
    time_gates_3site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_3body);
}

void flbit::create_lbit_transform_gates() {
    if(unitary_gates_2site_layers.size() == settings::model::lbit::u_layer) return;
    tools::log->info("Creating {} layers of 2-site unitary gates", settings::model::lbit::u_layer);
    unitary_gates_2site_layers.clear();
    for(size_t idx = 0; idx < settings::model::lbit::u_layer; idx++)
        unitary_gates_2site_layers.emplace_back(qm::lbit::get_unitary_2gate_layer(settings::model::model_size, settings::model::lbit::f_mixer));
}

void flbit::transform_to_real_basis() {
    if(unitary_gates_2site_layers.size() != settings::model::lbit::u_layer) create_lbit_transform_gates();
    auto t_map    = tid::tic_scope("l2r");
    tensors.state = std::make_unique<StateFinite>(*state_lbit);
    tensors.state->set_name("state_real");
    tools::log->debug("Transforming {} to {} using {} unitary layers", state_lbit->get_name(), tensors.state->get_name(), unitary_gates_2site_layers.size());
    for(const auto &layer : unitary_gates_2site_layers) tools::finite::mps::apply_gates(*tensors.state, layer, false, status.chi_lim);

    tensors.clear_measurements();
    tensors.clear_cache();
    [[maybe_unused]] auto has_normalized = tools::finite::mps::normalize_state(*tensors.state, status.chi_lim, std::nullopt, NormPolicy::IFNEEDED);
    status.position                      = tensors.get_position<long>();
    status.direction                     = tensors.state->get_direction();

    if constexpr(settings::debug) {
        auto t_dbg = tid::tic_scope("debug");
        // Double check the transform operation
        // Check that the transform backwards is equal to to the original state
        auto state_lbit_debug = *tensors.state;
        for(auto &layer : iter::reverse(unitary_gates_2site_layers)) tools::finite::mps::apply_gates(state_lbit_debug, layer, true, status.chi_lim);
        auto overlap = tools::finite::ops::overlap(*state_lbit, state_lbit_debug);
        tools::log->info("Debug overlap: {:.16f}", overlap);
        if(std::abs(overlap - 1) > 1e-10) throw std::runtime_error(fmt::format("State overlap after transform back from real is not 1: Got {:.16f}", overlap));
    }
}

void flbit::transform_to_lbit_basis() {
    if(unitary_gates_2site_layers.size() != settings::model::lbit::u_layer) create_lbit_transform_gates();
    auto t_map = tid::tic_scope("r2l");
    state_lbit = std::make_unique<StateFinite>(*tensors.state);
    state_lbit->set_name("state_lbit");

    tools::log->debug("Transforming {} to {}", tensors.state->get_name(), state_lbit->get_name());
    state_lbit->clear_cache();
    state_lbit->clear_measurements();
    for(auto &layer : iter::reverse(unitary_gates_2site_layers)) tools::finite::mps::apply_gates(*state_lbit, layer, true, status.chi_lim);
    [[maybe_unused]] auto has_normalized = tools::finite::mps::normalize_state(*state_lbit, status.chi_lim, std::nullopt, NormPolicy::IFNEEDED);
    if constexpr(settings::debug) {
        //        auto t_dbg = tid::tic_scope("debug");
        // Double check the transform operation
        // Check that the transform backwards is equal to to the original state
        auto state_real_debug = *state_lbit;
        for(auto &layer : unitary_gates_2site_layers) tools::finite::mps::apply_gates(state_real_debug, layer, false, status.chi_lim);
        auto overlap = tools::finite::ops::overlap(*tensors.state, state_real_debug);
        tools::log->info("Debug overlap: {:.16f}", overlap);
        if(std::abs(overlap - 1) > 1e-10) throw std::runtime_error(fmt::format("State overlap after transform back from lbit is not 1: Got {:.16f}", overlap));
    }
}

void flbit::write_to_file(StorageReason storage_reason, std::optional<CopyPolicy> copy_file) {
    tensors.clear_cache();
    tensors.clear_measurements();
    AlgorithmFinite::write_to_file(storage_reason, *tensors.state, *tensors.model, *tensors.edges, copy_file);

    // Save the unitaries once
    if(storage_reason == StorageReason::MODEL) {
        auto        t_h5      = tid::tic_scope("h5");
        auto        t_model   = tid::tic_scope("MODEL");
        auto        t_gates   = tid::tic_scope("gates");
        std::string grouppath = "/fLBIT/model/unitary_gates";
        if(h5pp_file->linkExists(grouppath)) return;
        for(const auto &[idx_layer, layer] : iter::enumerate(unitary_gates_2site_layers)) {
            std::string layerpath = fmt::format("{}/layer_{}", grouppath, idx_layer);
            for(const auto &[idx_gate, u] : iter::enumerate(layer)) {
                std::string gatepath = fmt::format("{}/u_{}", layerpath, idx_gate);
                h5pp_file->writeDataset(u.op, gatepath);
                h5pp_file->writeAttribute(u.dim, "dim", gatepath);
                h5pp_file->writeAttribute(u.pos, "pos", gatepath);
            }
            h5pp_file->writeAttribute(static_cast<size_t>(idx_layer), "idx_layer", layerpath);
            h5pp_file->writeAttribute(layer.size(), "num_gates", layerpath);
        }
        h5pp_file->writeAttribute(unitary_gates_2site_layers.size(), "num_layers", grouppath);
    }

    // Save the lbit analysis once
    if(storage_reason == StorageReason::MODEL) {
        auto t_h5       = tid::tic_scope("h5");
        auto t_model    = tid::tic_scope("MODEL");
        auto t_analysis = tid::tic_scope("analysis");
        if(h5pp_file->linkExists("/fLBIT/analysis")) return;
        std::vector<size_t> urange;
        std::vector<double> frange;
        size_t              sample = 1;
        if(settings::flbit::compute_lbit_stats) {
            sample = 50;
            urange = num::range<size_t>(1, 4);
            frange = num::range<double>(0, 0.8, 0.05);
        } else if(settings::flbit::compute_lbit_length) {
            urange = {settings::model::lbit::u_layer};
            frange = {settings::model::lbit::f_mixer};
        }
        if(not urange.empty() and not frange.empty()) {
            auto [cls_avg, sse_avg, decay, lioms] = qm::lbit::get_lbit_analysis(urange, frange, tensors.get_length(), sample);
            h5pp_file->writeDataset(cls_avg, "/fLBIT/analysis/cls_avg");
            h5pp_file->writeDataset(sse_avg, "/fLBIT/analysis/sse_avg");
            h5pp_file->writeDataset(decay, "/fLBIT/analysis/decay");
            h5pp_file->writeDataset(lioms, "/fLBIT/analysis/lioms");
            h5pp_file->writeAttribute(urange, "u_depth", "/fLBIT/analysis/cls_avg");
            h5pp_file->writeAttribute(urange, "u_depth", "/fLBIT/analysis/sse_avg");
            h5pp_file->writeAttribute(urange, "u_depth", "/fLBIT/analysis/decay");
            h5pp_file->writeAttribute(urange, "u_depth", "/fLBIT/analysis/lioms");
            h5pp_file->writeAttribute(frange, "f_mixer", "/fLBIT/analysis/cls_avg");
            h5pp_file->writeAttribute(frange, "f_mixer", "/fLBIT/analysis/sse_avg");
            h5pp_file->writeAttribute(frange, "f_mixer", "/fLBIT/analysis/decay");
            h5pp_file->writeAttribute(frange, "f_mixer", "/fLBIT/analysis/lioms");
            h5pp_file->writeAttribute(sample, "samples", "/fLBIT/analysis/cls_avg");
            h5pp_file->writeAttribute(sample, "samples", "/fLBIT/analysis/sse_avg");
            h5pp_file->writeAttribute(sample, "samples", "/fLBIT/analysis/decay");
            h5pp_file->writeAttribute(sample, "samples", "/fLBIT/analysis/lioms");
        }
    }
}

void flbit::print_status_update() {
    if(num::mod(status.iter, settings::print_freq(status.algo_type)) != 0) return;
    if(settings::print_freq(status.algo_type) == 0) return;
    auto        t_print = tid::tic_scope("print");
    std::string report;
    report += fmt::format("{:<} ", tensors.state->get_name());
    report += fmt::format("iter:{:<4} ", status.iter);
    report += fmt::format("step:{:<5} ", status.step);
    report += fmt::format("L:{} ", tensors.get_length());
    if(tensors.active_sites.empty())
        report += fmt::format("l:{:<2} ", tensors.get_position());
    else if(tensors.state->get_direction() > 0)
        report += fmt::format("l:[{:>2}-{:<2}] ", tensors.active_sites.front(), tensors.active_sites.back());
    else if(tensors.state->get_direction() < 0)
        report += fmt::format("l:[{:>2}-{:<2}] ", tensors.active_sites.back(), tensors.active_sites.front());
    report += fmt::format("E/L:{:<20.16f} ", tools::finite::measure::energy_per_site(tensors));
    report += fmt::format("Sₑ(L/2):{:<10.8f} ", tools::finite::measure::entanglement_entropy_midchain(*tensors.state));
    if(tensors.state->measurements.number_entropy_midchain) // This one is expensive
        report += fmt::format("Sₙ(L/2):{:<10.8f} ", tensors.state->measurements.number_entropy_midchain.value());

    if(state_lbit->measurements.number_entropy_midchain) // This one is expensive
        report += fmt::format("Sₙ(L/2) lbit:{:<10.8f} ", state_lbit->measurements.number_entropy_midchain.value());

    report += fmt::format("χ:{:<3}|{:<3}|{:<3} ", settings::chi_lim_max(status.algo_type), status.chi_lim,
                          tools::finite::measure::bond_dimension_midchain(*tensors.state));
    report += fmt::format("log₁₀trnc:{:<8.4f} ", std::log10(tensors.state->get_truncation_error_midchain()));
    report += fmt::format("wtime:{:<} ", fmt::format("{:>6.2f}s", tid::get_unscoped("t_tot").get_time()));
    report += fmt::format("ptime:{:<} ", fmt::format("{:+>8.2e}s", status.phys_time));
    report += fmt::format("itime:{:<} ", fmt::format("{:>4.2f}s", tid::get_unscoped("fLBIT")["run"].get_lap()));
    report += fmt::format("mem[rss {:<.1f}|peak {:<.1f}|vm {:<.1f}]MB ", tools::common::profile::mem_rss_in_mb(), tools::common::profile::mem_hwm_in_mb(),
                          tools::common::profile::mem_vm_in_mb());
    tools::log->info(report);
}