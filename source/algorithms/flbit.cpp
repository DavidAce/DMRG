#include "flbit.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "debug/info.h"
#include "general/iter.h"
#include "io/fmt.h"
#include "math/num.h"
#include "math/tenx.h"
#include "qm/lbit.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include "tools/common/prof.h"
#include "tools/finite/h5.h"
#include "tools/finite/measure.h"
#include "tools/finite/mps.h"
#include "tools/finite/ops.h"
#include "tools/finite/print.h"
#include <h5pp/h5pp.h>
#include <unsupported/Eigen/CXX11/Tensor>

flbit::flbit(std::shared_ptr<h5pp::File> h5file_) : AlgorithmFinite(std::move(h5file_), AlgorithmType::fLBIT) {
    tools::log->trace("Constructing class_flbit");
    tensors.state->set_name("state_real");
}

void flbit::resume() {
    // Resume can imply many things
    // 1) Resume a simulation which terminated prematurely
    // 2) Resume a previously successful simulation. This may be desireable if the config
    //    wants something that is not present in the file.
    //      a) A certain number of states
    //      b) A state inside a particular energy window
    //      c) The ground or "roof" states
    // To guide the behavior, we check the setting ResumePolicy.

    auto state_name = settings::storage::file_resume_name;
    if(state_name.empty()) state_name = "state_real";

    auto state_prefixes = tools::common::h5::resume::find_state_prefixes(*h5file, status.algo_type, state_name);
    if(state_prefixes.empty()) throw except::state_error("no resumable states were found");
    for(const auto &state_prefix : state_prefixes) {
        tools::log->info("Resuming state [{}]", state_prefix);
        try {
            tools::finite::h5::load::simulation(*h5file, state_prefix, tensors, status, status.algo_type);
        } catch(const except::load_error &le) { continue; }

        // Our first task is to decide on a state name for the newly loaded state
        // The simplest is to inferr it from the state prefix itself
        auto name = tools::common::h5::resume::extract_state_name(state_prefix);
        tensors.state->set_name(name);
        clear_convergence_status();

        // Create an initial state in the real basis
        auto bond_lim = settings::get_bond_init(status.algo_type);
        switch(settings::strategy::initial_state) {
            case StateInit::PRODUCT_STATE_NEEL:
            case StateInit::RANDOM_PRODUCT_STATE_ZEROMAG: break;
            default:
                tools::log->warn("Expected initial_state: PRODUCT_STATE_NEEL|RANDOM_PRODUCT_STATE_ZEROMAG. Got {}", enum2sv(settings::strategy::initial_state));
        }

        if(settings::strategy::initial_axis != "+z") tools::log->warn("Expected initial_axis == +z. Got {}", settings::strategy::initial_axis);

        tensors.randomize_state(ResetReason::INIT, settings::strategy::initial_state, StateInitType::REAL, settings::strategy::initial_axis,
                                settings::strategy::use_eigenspinors, settings::input::bitfield, bond_lim);

        tensors.move_center_point_to_inward_edge();

        // Load the unitaries
        unitary_gates_2site_layers.clear();
        std::string grouppath = "/fLBIT/model/unitary_gates";
        if(not h5file->linkExists(grouppath)) throw except::runtime_error("Missing link: {}", grouppath);
        auto num_layers = h5file->readAttribute<size_t>(grouppath, "num_layers");
        if(num_layers != settings::model::lbit::u_layer)
            throw except::runtime_error("Mismatch in number of layers: file {} != cfg {}", num_layers != settings::model::lbit::u_layer);
        unitary_gates_2site_layers.resize(num_layers);
        for(auto &&[idx_layer, layer] : iter::enumerate(unitary_gates_2site_layers)) {
            std::string layerpath = fmt::format("{}/layer_{}", grouppath, idx_layer);
            auto        num_gates = h5file->readAttribute<size_t>(layerpath, "num_gates");
            layer.resize(num_gates);
            for(auto &&[idx_gate, u] : iter::enumerate(layer)) {
                std::string gatepath = fmt::format("{}/u_{}", layerpath, idx_gate);
                auto        op       = h5file->readDataset<Eigen::Tensor<Scalar, 2>>(gatepath);
                auto        pos      = h5file->readAttribute<std::vector<size_t>>(gatepath, "pos");
                auto        dim      = h5file->readAttribute<std::vector<long>>(gatepath, "dim");
                u                    = qm::Gate(op, pos, dim);
            }
        }

        // Generate the corresponding state in lbit basis
        transform_to_lbit_basis();

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
                throw except::runtime_error("Unrecognized state name for flbit: [{}]", name);
            task_list.emplace_back(flbit_task::POST_DEFAULT);
        }
        run_task_list(task_list);
    }

    // If we reached this point the current state has finished for one reason or another.
    // TODO: We may still have some more things to do, e.g. the config may be asking for more states
}

void flbit::run_task_list(std::deque<flbit_task> &task_list) {
    while(not task_list.empty()) {
        auto task = task_list.front();
        switch(task) {
            case flbit_task::INIT_RANDOMIZE_MODEL: randomize_model(); break;
            case flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE: randomize_state(ResetReason::INIT, StateInit::RANDOM_PRODUCT_STATE); break;
            case flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE_ZEROMAG: randomize_state(ResetReason::INIT, StateInit::RANDOM_PRODUCT_STATE_ZEROMAG); break;
            case flbit_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE: randomize_state(ResetReason::INIT, StateInit::RANDOM_ENTANGLED_STATE); break;
            case flbit_task::INIT_BOND_LIMITS: init_bond_dimension_limits(); break;
            case flbit_task::INIT_TRNC_LIMITS: init_truncation_error_limits(); break;
            case flbit_task::INIT_WRITE_MODEL: write_to_file(StorageEvent::MODEL); break;
            case flbit_task::INIT_CLEAR_STATUS: status.clear(); break;
            case flbit_task::INIT_CLEAR_CONVERGENCE: clear_convergence_status(); break;
            case flbit_task::INIT_DEFAULT: run_preprocessing(); break;
            case flbit_task::INIT_TIME: {
                create_time_points();
                update_time_step();
                break;
            }
            case flbit_task::INIT_GATES: {
                if(settings::flbit::use_swap_gates) {
                    create_hamiltonian_swap_gates();
                    update_time_evolution_swap_gates();
                } else {
                    create_hamiltonian_gates();
                    update_time_evolution_gates();
                }
                create_unitary_circuit_gates();
                break;
            }
            case flbit_task::TIME_EVOLVE:
                tensors.state->set_name("state_real");
                run_algorithm();
                break;
            case flbit_task::TRANSFORM_TO_LBIT: transform_to_lbit_basis(); break;
            case flbit_task::TRANSFORM_TO_REAL: transform_to_real_basis(); break;
            case flbit_task::POST_WRITE_RESULT: write_to_file(StorageEvent::LAST_STATE); break;
            case flbit_task::POST_PRINT_RESULT: print_status_full(); break;
            case flbit_task::POST_PRINT_TIMERS: tools::common::timer::print_timers(); break;
            case flbit_task::POST_FES_ANALYSIS: run_fes_analysis(); break;
            case flbit_task::POST_DEFAULT: run_postprocessing(); break;
            case flbit_task::TIMER_RESET: tid::reset("fLBIT"); break;
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
        throw except::runtime_error("Simulation ended with unfinished tasks");
    }
}

void flbit::run_preprocessing() {
    tools::log->info("Running {} preprocessing", status.algo_type_sv());
    auto t_pre = tid::tic_scope("pre");
    status.clear();
    randomize_model(); // First use of random!
    init_bond_dimension_limits();
    init_truncation_error_limits();

    // Create an initial state in the real basis
    randomize_state(ResetReason::INIT, settings::strategy::initial_state);
    tensors.move_center_point_to_inward_edge();
    tools::finite::print::model(*tensors.model);
    create_time_points();
    update_time_step();
    if(settings::model::model_size <= 6) {
        // Create a copy of the state as a full state vector for ED comparison
        auto list_Lsite = num::range<size_t>(0, settings::model::model_size, 1);
        Upsi_ed         = tools::finite::measure::mps_wavefn(*tensors.state);
        tools::log->info("<Ψ_ed|Ψ_ed>   : {:.16f}", tenx::VectorMap(Upsi_ed).norm());
        auto nbody = std::vector<size_t>{1, 2, 3};
        ham_gates_Lsite.emplace_back(qm::Gate(tensors.model->get_multisite_ham(list_Lsite, nbody), list_Lsite, tensors.state->get_spin_dims(list_Lsite)));
        time_gates_Lsite = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_Lsite, settings::flbit::time_gate_id_threshold);
    }

    if(settings::flbit::use_swap_gates) {
        create_hamiltonian_swap_gates();
        update_time_evolution_swap_gates();
    } else {
        create_hamiltonian_gates();
        update_time_evolution_gates();
    }

    create_unitary_circuit_gates();

    if(not tensors.position_is_inward_edge()) throw except::logic_error("Put the state on an edge!");

    // Generate the corresponding state in lbit basis
    transform_to_lbit_basis();
    write_to_file(StorageEvent::MODEL, CopyPolicy::TRY);
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
    if(not tensors.position_is_inward_edge()) throw except::logic_error("Put the state on an edge!");
    while(true) {
        update_state();
        check_convergence();
        write_to_file(StorageEvent::ITER_STATE);
        print_status();
        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);

        // It's important not to perform the last move, so we break now: that last state would not get optimized
        if(status.algo_stop != AlgorithmStop::NONE) break;
        if(settings::flbit::use_swap_gates) {
            update_time_evolution_swap_gates();
            // #pragma message "no need to update the normal gates"
            //             update_time_evolution_gates();

        } else
            update_time_evolution_gates();

        status.wall_time = tid::get_unscoped("t_tot").get_time();
        status.algo_time = t_run->get_time();
        t_run->start_lap();
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", status.algo_type_sv(), tensors.state->get_name(), status.algo_stop_sv());
    status.algorithm_has_finished = true;
}

void flbit::run_fes_analysis() {
    if(settings::strategy::fes_rate == 0) return;
    tools::log->warn("FES is not yet implemented for flbit");
}

void flbit::update_state() {
    /*!
     * \fn void update_state()
     */
    tools::log->debug("Starting fLBIT: iter {} | Δt = ({:.2e}, {:.2e})", status.iter, std::real(status.delta_t), std::imag(status.delta_t));
    if(not state_lbit) throw except::logic_error("state_lbit == nullptr: Set the state in lbit basis before running an flbit step");
    if(not state_lbit_init) {
        state_lbit_init = std::make_unique<StateFinite>(*state_lbit);
        tools::finite::mps::normalize_state(*state_lbit_init, std::nullopt, NormPolicy::ALWAYS);
    }
    *state_lbit = *state_lbit_init;

    // Time evolve from 0 to time_point[iter] here
    auto t_step  = tid::tic_scope("step");
    auto t_evo   = tid::tic_scope("time_evo");
    auto svd_cfg = svd::config(status.bond_lim, status.trnc_lim);

    if(settings::flbit::use_swap_gates) {
        tools::log->debug("Applying time evolution swap gates Δt = ({:.2e}, {:.2e})", std::real(status.delta_t), std::imag(status.delta_t));
        tools::finite::mps::apply_swap_gates(*state_lbit, time_swap_gates_1site, false, GateMove::AUTO, svd_cfg);
        tools::finite::mps::apply_swap_gates(*state_lbit, time_swap_gates_2site, false, GateMove::AUTO, svd_cfg);
        tools::finite::mps::apply_swap_gates(*state_lbit, time_swap_gates_3site, false, GateMove::AUTO, svd_cfg);
    } else {
        tools::log->debug("Applying time evolution gates Δt = ({:.2e}, {:.2e})", std::real(status.delta_t), std::imag(status.delta_t));
        tools::finite::mps::apply_gates(*state_lbit, time_gates_1site, false, GateMove::AUTO, svd_cfg);
        tools::finite::mps::apply_gates(*state_lbit, time_gates_2site, false, GateMove::AUTO, svd_cfg);
        tools::finite::mps::apply_gates(*state_lbit, time_gates_3site, false, GateMove::AUTO, svd_cfg);
    }
    tools::finite::mps::normalize_state(*state_lbit, std::nullopt, NormPolicy::IFNEEDED);

    t_evo.toc();
    transform_to_real_basis();

    if constexpr(settings::debug) {
        if(tools::log->level() <= spdlog::level::debug) {
            tools::log->debug("{} bond dimensions: {}", state_lbit->get_name(), tools::finite::measure::bond_dimensions(*state_lbit));
            tools::log->debug("{} bond dimensions: {}", tensors.state->get_name(), tools::finite::measure::bond_dimensions(*tensors.state));
        }
        for(const auto &mps : state_lbit->mps_sites) mps->assert_normalized();
        if(settings::model::model_size <= 6) {
            if(Upsi_ed.dimension(0) != time_gates_Lsite[0].op.dimension(1))
                throw except::logic_error("Upsi_ed may not have been initialized: Upsi_ed: {}", Upsi_ed.dimensions());
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
    }

    tensors.clear_measurements();
    tensors.clear_cache();
    status.phys_time = std::abs(time_points[std::min(status.iter, time_points.size() - 1)]);

    status.iter += 1;
    status.step += settings::model::model_size;
    status.position  = tensors.get_position<long>();
    status.direction = tensors.state->get_direction();
}

void flbit::update_time_step() {
    if(time_points.empty()) create_time_points();
    tools::log->trace("Updating time step");
    auto t_updtstep = tid::tic_scope("update_time_step");
    if(status.iter >= time_points.size()) {
        status.delta_t   = time_points.back();
        status.algo_stop = AlgorithmStop::SUCCESS;
        return;
    }
    status.delta_t = time_points[status.iter];
    if(std::abs(status.delta_t) == 0) throw except::logic_error("Expected nonzero delta_t after time step update");
    tools::log->debug("Time step iter {} | Δt = {} | t = {:8.2e}", status.iter, status.delta_t, status.phys_time);
}

void flbit::check_convergence() {
    if(not tensors.position_is_inward_edge()) return;
    auto t_con = tid::tic_scope("check_conv");
    check_convergence_entg_entropy();
    if(status.entanglement_saturated_for > 0)
        status.algorithm_saturated_for++;
    else
        status.algorithm_saturated_for = 0;

    status.algorithm_converged_for = status.iter + 1 - std::min(settings::flbit::time_num_steps, status.iter + 1);
    status.algorithm_has_succeeded = status.algorithm_converged_for > 0;
    if(status.algorithm_saturated_for > settings::strategy::min_saturation_iters and status.algorithm_converged_for == 0)
        status.algorithm_has_stuck_for++;
    else
        status.algorithm_has_stuck_for = 0;

    status.algorithm_has_to_stop = false; // Never stop due to saturation

    tools::log->debug("Simulation report: converged {} | saturated {} | stuck {} | succeeded {} | has to stop {}", status.algorithm_converged_for,
                      status.algorithm_saturated_for, status.algorithm_has_stuck_for, status.algorithm_has_succeeded, status.algorithm_has_to_stop);
    status.algo_stop = AlgorithmStop::NONE;
    if(status.iter >= settings::flbit::min_iters) {
        if(status.iter >= settings::flbit::max_iters) status.algo_stop = AlgorithmStop::MAX_ITERS;
        if(status.iter >= settings::flbit::time_num_steps) status.algo_stop = AlgorithmStop::SUCCESS;
        if(status.algorithm_has_to_stop) status.algo_stop = AlgorithmStop::SATURATED;
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
    if(std::abs(time_diff) == 0) throw except::logic_error("time_start - time_final == 0");
    // Check that the time limits are purely real or imaginary!
    bool time_is_real = std::abs(time_diff.real()) > 0;
    bool time_is_imag = std::abs(time_diff.imag()) > 0;
    if(time_is_real and time_is_imag)
        throw except::logic_error("time_start and time_final must both be either purely real or imaginary. Got:\n"
                                  "time_start = {:.8f}{:+.8f}\n"
                                  "time_final = {:.8f}{:+.8f}",
                                  time_start.real(), time_start.imag(), time_final.real(), time_final.imag());
    time_points.reserve(settings::flbit::time_num_steps);
    if(time_is_real)
        for(const auto &t : num::LogSpaced(settings::flbit::time_num_steps, std::real(time_start), std::real(time_final))) {
            if(std::isinf(t) or std::isnan(t)) throw except::runtime_error("Invalid time point: {}", t);
            time_points.emplace_back(t);
        }
    else
        for(const auto &t : num::LogSpaced(settings::flbit::time_num_steps, std::imag(time_start), std::imag(time_final))) {
            if(std::isinf(t) or std::isnan(t)) throw except::runtime_error("Invalid time point: {}", t);
            time_points.emplace_back(std::complex<double>(0, t));
        }
    tools::log->trace(FMT_STRING("Created {} time points:\n{}"), time_points.size(), time_points);
    // Sanity check
    if(time_points.size() != settings::flbit::time_num_steps)
        throw except::logic_error("Got time_points.size():[{}] != settings::flbit::time_num_steps:[{}]", time_points.size(), settings::flbit::time_num_steps);
}

void flbit::create_hamiltonian_gates() {
    if(ready_hamiltonian_gates) throw except::logic_error("Hamiltonian gates have already been constructed");
    tools::log->info("Creating Hamiltonian gates");
    auto t_hamgates = tid::tic_scope("create_ham_gates");
    ham_gates_1body.clear();
    ham_gates_2body.clear();
    ham_gates_3body.clear();

    // Create the hamiltonian gates with n-site terms
    //
    // Note for 2-body terms
    // We want MPO's with only 2-body terms of the Hamiltonian, for distances |i-j| <= J2_span sites.
    // Let's assume we have L==8 and J2_span == 3, then we apply time-evolution mpo's
    // L     :  0,1,2,3,4,5,6,7
    // mpo[0]: [0,1,2,3]
    // mpo[1]:   [1,2,3,4]
    // mpo[2]:     [2,3,4,5]
    // mpo[3]:       [3,4,5,6]
    // mpo[4]:         [4,5,6,7]
    // Note that the interaction on sites {1,2} gets applied twice, and {3,4} thrice.
    // To compensate for this, pass a "0" to nbody, which tells the mpo-generator to divide J[i,j] by
    // the number of times it is applied. This ensures that each interaction is time-evolved correctly.

    auto L          = settings::model::model_size;
    auto list_1site = num::range<size_t>(0, L - 0, 1);
    auto list_2site = num::range<size_t>(0, L - 1, 1);
    auto list_3site = num::range<size_t>(0, L - 2, 1);
    for(auto pos : list_1site) {
        auto sites = std::vector<size_t>{pos};            // A list of site indices
        auto nbody = std::vector<size_t>{1};              // A list of included nbody interaction terms (1: on-site terms, 2: pairwise, and so on)
        auto spins = tensors.state->get_spin_dims(sites); // A list of spin dimensions for each site (should all be 2 for two-level systems)
        tools::log->info("Generating {}-body hamiltonian on sites {}", nbody, sites);
        ham_gates_1body.emplace_back(qm::Gate(tensors.model->get_multisite_ham({pos}, nbody), sites, spins));
    }

    auto J2_ctof = std::min(settings::model::lbit::J2_span, L - 1); // Max distance |i-j| to the furthest interacting site L-1
    for(auto posL : list_2site) {
        auto posR = posL + J2_ctof;
        if(J2_ctof == 0) break;
        if(posR >= L) break;
        auto sites = num::range<size_t>(posL, posR + 1); // +1 to include last site, which is otherwise not included in range
        auto nbody = std::vector<size_t>{0, 2};          // zero is a flag to enable compensation for double-counting
        auto spins = tensors.state->get_spin_dims(sites);
        tools::log->info("Generating {}-body hamiltonian on sites {}", nbody, sites);
        ham_gates_2body.emplace_back(qm::Gate(tensors.model->get_multisite_ham(sites, nbody), sites, spins));
    }

    for(auto posL : list_3site) {
        auto range = 2ul; // Distance to next-nearest neighbor when 3 sites interact
        auto posR  = posL + range;
        if(posR >= L) break;
        auto sites = num::range<size_t>(posL, posR + 1); // +1 to include last site
        auto nbody = std::vector<size_t>{3};
        auto spins = tensors.state->get_spin_dims(sites);
        tools::log->info("Generating {}-body hamiltonian on sites {}", nbody, sites);
        ham_gates_3body.emplace_back(qm::Gate(tensors.model->get_multisite_ham(sites, nbody), sites, spins));
    }

    for(const auto &[idx, ham] : iter::enumerate(ham_gates_1body))
        if(tenx::isZero(ham.op)) throw except::runtime_error("ham1[{}] is all zeros", idx);

    for(const auto &[idx, ham] : iter::enumerate(ham_gates_2body))
        if(tenx::isZero(ham.op)) throw except::runtime_error("ham2[{}] is all zeros", idx);

    for(const auto &[idx, ham] : iter::enumerate(ham_gates_3body))
        if(tenx::isZero(ham.op)) throw except::runtime_error("ham3[{}] is all zeros", idx);
    ready_hamiltonian_gates = true;
}

void flbit::create_hamiltonian_swap_gates() {
    if(ready_hamiltonian_swap_gates) throw except::logic_error("Hamiltonian swap gates have already been constructed");
    tools::log->info("Creating Hamiltonian swap gates");
    auto t_swaphamgates = tid::tic_scope("create_ham_gates_swap");
    ham_swap_gates_1body.clear();
    ham_swap_gates_2body.clear();
    ham_swap_gates_3body.clear();

    // Create the hamiltonian 2-site swap gates with n-site range
    auto L          = settings::model::model_size;
    auto list_1site = num::range<size_t>(0, L - 0, 1);
    auto list_2site = num::range<size_t>(0, L - 1, 1);
    auto list_3site = num::range<size_t>(0, L - 2, 1);
    for(auto pos : list_1site) {
        auto range = 0ul;                                      // Number of site indices
        auto sites = num::range<size_t>(pos, pos + range + 1); // A list of site indices. +1 to include last site (in this case the only site)
        auto nbody = std::vector<size_t>{1};                   // A list of included nbody interaction terms (1: on-site terms, 2: pairwise, and so on)
        auto spins = tensors.state->get_spin_dims(sites);      // A list of spin dimensions for each site (should all be 2 for two-level systems)
        tools::log->info("Generating {}-body hamiltonian on sites {}", nbody, sites);
        ham_swap_gates_1body.emplace_back(qm::SwapGate(tensors.model->get_multisite_ham({pos}, nbody), sites, spins));
    }
    auto J2_ctof = std::min(settings::model::lbit::J2_span, L - 1); // Max distance |i-j| to the furthest interacting site L-1
    for(auto posL : list_2site) {
        auto maxR = std::min<size_t>(posL + J2_ctof, L - 1);
        if(maxR == posL) continue;
        for(auto posR : num::range<size_t>(posL, maxR + 1)) { // maxR+1 to include the last site in range
            if(posL == posR) continue;
            if(posL >= L) throw except::logic_error("posL {} >= L {}", posL, L);
            if(posR >= L) throw except::logic_error("posR {} >= L {}", posR, L);
            auto sites = std::vector<size_t>{posL, posR}; // +1 to include last site
            auto nbody = std::vector<size_t>{2};
            auto spins = tensors.state->get_spin_dims(sites);
            tools::log->info("Generating {}-body hamiltonian on sites {}", nbody, sites);
            auto ham_swap_gate = qm::SwapGate(tensors.model->get_multisite_ham(sites, nbody), sites, spins);
            // Accept all swap gates even if all elements are near zero on gates for remote sites,
            // since at large times t these can become relevant again by exp(-itH)
            ham_swap_gates_2body.emplace_back(ham_swap_gate);
        }
    }
    // Ignore Hamiltonians with entries smaller than J2_zero: the timescale is too small to resolve them.

    for(auto pos : list_3site) {
        auto range = 2ul;                                                                  // Distance to next-nearest neighbor when 3 sites interact
        auto sites = num::range<size_t>(pos, std::clamp<size_t>(pos + range + 1, 0ul, L)); // +1 to include last site
        auto nbody = std::vector<size_t>{3};
        auto spins = tensors.state->get_spin_dims(sites);
        tools::log->info("Generating {}-body hamiltonian on sites {}", nbody, sites);
        ham_swap_gates_3body.emplace_back(qm::SwapGate(tensors.model->get_multisite_ham(sites, nbody), sites, spins));
    }
    for(const auto &[idx, ham] : iter::enumerate(ham_swap_gates_1body)) {
        if(tenx::isZero(ham.op)) throw except::runtime_error("ham1[{}] is all zeros", idx);
    }
    for(const auto &[idx, ham] : iter::enumerate(ham_swap_gates_2body)) {
        if(tenx::isZero(ham.op)) tools::log->debug("hamiltonian 2-body swap gate {} for sites {} is a zero-matrix", idx, ham.pos);
    }
    for(const auto &[idx, ham] : iter::enumerate(ham_swap_gates_3body))
        if(tenx::isZero(ham.op)) throw except::runtime_error("ham3[{}] is all zeros", idx);

    ready_hamiltonian_swap_gates = true;
}

void flbit::update_time_evolution_gates() {
    // Create the time evolution operators
    if(time_points.empty()) create_time_points();
    if(not ready_hamiltonian_gates) throw except::logic_error("Hamiltonian gates have not been constructed");
    update_time_step();
    tools::log->debug("Updating time evolution gates to iter {} | Δt = ({:.2e}, {:.2e})", status.iter, std::real(status.delta_t), std::imag(status.delta_t));
    time_gates_1site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_1body, settings::flbit::time_gate_id_threshold);
    time_gates_2site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_2body, settings::flbit::time_gate_id_threshold);
    time_gates_3site = qm::lbit::get_time_evolution_gates(status.delta_t, ham_gates_3body, settings::flbit::time_gate_id_threshold);
}

void flbit::update_time_evolution_swap_gates() {
    // Create the time evolution operators
    if(time_points.empty()) create_time_points();
    if(not ready_hamiltonian_swap_gates) throw except::logic_error("Hamiltonian swap gates have not been constructed");
    update_time_step();
    tools::log->debug("Updating time evolution swap gates to iter {} | Δt = ({:.2e}, {:.2e})", status.iter, std::real(status.delta_t),
                      std::imag(status.delta_t));
    time_swap_gates_1site = qm::lbit::get_time_evolution_swap_gates(status.delta_t, ham_swap_gates_1body, settings::flbit::time_gate_id_threshold);
    time_swap_gates_2site = qm::lbit::get_time_evolution_swap_gates(status.delta_t, ham_swap_gates_2body, settings::flbit::time_gate_id_threshold);
    time_swap_gates_3site = qm::lbit::get_time_evolution_swap_gates(status.delta_t, ham_swap_gates_3body, settings::flbit::time_gate_id_threshold);
}

void flbit::create_unitary_circuit_gates() {
    if(unitary_gates_2site_layers.size() == settings::model::lbit::u_layer) return;
    tools::log->info("Creating {} layers of 2-site unitary gates", settings::model::lbit::u_layer);
    unitary_gates_2site_layers.clear();
    std::vector<double> fields;
    for(const auto &field : tensors.model->get_parameter("J1_rand")) fields.emplace_back(std::any_cast<double>(field));
    tools::log->info("Fields: {}", fields);
    for(size_t idx = 0; idx < settings::model::lbit::u_layer; idx++)
        unitary_gates_2site_layers.emplace_back(
            qm::lbit::get_unitary_2gate_layer_blocked(settings::model::model_size, settings::model::lbit::f_mixer, fields, settings::model::lbit::J1_wdth));
}

void flbit::transform_to_real_basis() {
    if(unitary_gates_2site_layers.size() != settings::model::lbit::u_layer) create_unitary_circuit_gates();
    auto t_map    = tid::tic_scope("l2r");
    tensors.state = std::make_unique<StateFinite>(*state_lbit);
    tensors.state->set_name("state_real");
    tools::log->debug("Transforming {} to {} using {} unitary layers", state_lbit->get_name(), tensors.state->get_name(), unitary_gates_2site_layers.size());
    for(const auto &layer : unitary_gates_2site_layers)
        tools::finite::mps::apply_gates(*tensors.state, layer, false, GateMove::OFF, svd::config(status.bond_lim, status.trnc_lim)); // L16: true 29 | false
    for(const auto &layer : unitary_gates_2site_layers)
        for(const auto &u : layer) u.unmark_as_used();

    tensors.clear_measurements();
    tensors.clear_cache();
    tools::finite::mps::normalize_state(*tensors.state, std::nullopt, NormPolicy::IFNEEDED);
    status.position  = tensors.get_position<long>();
    status.direction = tensors.state->get_direction();

    if constexpr(settings::debug) {
        auto t_dbg = tid::tic_scope("debug");
        if(tools::log->level() <= spdlog::level::debug) {
            tools::log->debug("{} bond dimensions: {}", state_lbit->get_name(), tools::finite::measure::bond_dimensions(*state_lbit));
            tools::log->debug("{} bond dimensions: {}", tensors.state->get_name(), tools::finite::measure::bond_dimensions(*tensors.state));
        }
        // Check normalization
        for(const auto &mps : state_lbit->mps_sites) mps->assert_normalized();

        // Double-check the transform operation
        // Check that the transform backwards is equal to to the original state
        auto state_lbit_debug = *tensors.state;
        for(const auto &layer : iter::reverse(unitary_gates_2site_layers))
            tools::finite::mps::apply_gates(state_lbit_debug, layer, true, GateMove::AUTO, svd::config(status.bond_lim, status.trnc_lim));
        for(const auto &layer : iter::reverse(unitary_gates_2site_layers))
            for(const auto &u : layer) u.unmark_as_used();
        auto overlap = tools::finite::ops::overlap(*state_lbit, state_lbit_debug);
        tools::log->info("Debug overlap: {:.16f}", overlap);
        if(std::abs(overlap - 1) > 10 * status.trnc_lim)
            throw except::runtime_error("State overlap after transform back from real is not 1: Got {:.16f}", overlap);
    }
}

void flbit::transform_to_lbit_basis() {
    if(unitary_gates_2site_layers.size() != settings::model::lbit::u_layer) create_unitary_circuit_gates();
    auto t_map = tid::tic_scope("r2l");
    state_lbit = std::make_unique<StateFinite>(*tensors.state);
    state_lbit->set_name("state_lbit");

    tools::log->debug("Transforming {} to {}", tensors.state->get_name(), state_lbit->get_name());
    state_lbit->clear_cache();
    state_lbit->clear_measurements();
    for(const auto &layer : iter::reverse(unitary_gates_2site_layers))
        tools::finite::mps::apply_gates(*state_lbit, layer, true, GateMove::AUTO,
                                        svd::config(status.bond_lim, status.trnc_lim)); // L16: true 28 | false 29 svds
    for(const auto &layer : iter::reverse(unitary_gates_2site_layers))
        for(const auto &u : layer) u.unmark_as_used();
    tools::finite::mps::normalize_state(*state_lbit, std::nullopt, NormPolicy::IFNEEDED);

    if constexpr(settings::debug) {
        auto t_dbg = tid::tic_scope("debug");
        if(tools::log->level() <= spdlog::level::debug) {
            tools::log->debug("{} bond dimensions: {}", state_lbit->get_name(), tools::finite::measure::bond_dimensions(*state_lbit));
            tools::log->debug("{} bond dimensions: {}", tensors.state->get_name(), tools::finite::measure::bond_dimensions(*tensors.state));
        }
        // Check normalization
        for(const auto &mps : state_lbit->mps_sites) mps->assert_normalized();

        // Double-check the that transform operation backwards is equal to the original state
        auto state_real_debug = *state_lbit;
        for(auto &layer : unitary_gates_2site_layers)
            for(auto &g : layer) g.unmark_as_used();
        for(const auto &layer : unitary_gates_2site_layers)
            tools::finite::mps::apply_gates(state_real_debug, layer, false, GateMove::AUTO, svd::config(status.bond_lim, status.trnc_lim));
        for(const auto &layer : unitary_gates_2site_layers)
            for(const auto &u : layer) u.unmark_as_used();
        auto overlap = tools::finite::ops::overlap(*tensors.state, state_real_debug);
        tools::log->info("Debug overlap: {:.16f}", overlap);
        if(std::abs(overlap - 1) > 10 * status.trnc_lim)
            throw except::runtime_error("State overlap after transform back from lbit is not 1: Got {:.16f}", overlap);
    }
}

void flbit::write_to_file(StorageEvent storage_event, CopyPolicy copy_policy) {
    AlgorithmFinite::write_to_file(*tensors.state, *tensors.model, *tensors.edges, storage_event, copy_policy);
    // Save the unitaries once
    if(storage_event == StorageEvent::MODEL) {
        auto        t_h5      = tid::tic_scope("h5");
        auto        t_event   = tid::tic_scope(enum2sv(storage_event), tid::highest);
        auto        t_gates   = tid::tic_scope("gates");
        std::string grouppath = "/fLBIT/model/unitary_gates";
        if(h5file->linkExists(grouppath)) return;
        for(const auto &[idx_layer, layer] : iter::enumerate(unitary_gates_2site_layers)) {
            std::string layerpath = fmt::format("{}/layer_{}", grouppath, idx_layer);
            for(const auto &[idx_gate, u] : iter::enumerate(layer)) {
                std::string gatepath = fmt::format("{}/u_{}", layerpath, idx_gate);
                h5file->writeDataset(u.op, gatepath);
                h5file->writeAttribute(u.dim, gatepath, "dim");
                h5file->writeAttribute(u.pos, gatepath, "pos");
            }
            h5file->writeAttribute(static_cast<size_t>(idx_layer), layerpath, "idx_layer");
            h5file->writeAttribute(layer.size(), layerpath, "num_gates");
        }
        h5file->writeAttribute(unitary_gates_2site_layers.size(), grouppath, "num_layers");
    }

    // Save the lbit analysis once
    if(storage_event == StorageEvent::MODEL) {
        auto t_h5    = tid::tic_scope("h5");
        auto t_event = tid::tic_scope(enum2sv(storage_event), tid::highest);
        if(h5file->linkExists("/fLBIT/model/lbits")) return;
        auto urange = std::vector<size_t>{settings::model::lbit::u_layer};
        auto frange = std::vector<double>{settings::model::lbit::f_mixer};
        auto sample = settings::flbit::compute_lbit_stats;
        if(sample > 1) {
            urange = {8};   // num::range<size_t>(8, 9);
            frange = {0.5}; // num::range<double>(0.25, 0.50, 0.05);
        }
        if(sample > 0) {
            std::vector<double> fields;
            for(const auto &field : tensors.model->get_parameter("J1_rand")) fields.emplace_back(std::any_cast<double>(field));
            auto ugate_props = qm::lbit::UnitaryGateProperties{settings::model::lbit::ugate_type, fields, settings::model::lbit::J1_wdth};
            tools::log->info("Computing the lbit characteristic length-scale with {} gates", enum2sv(settings::model::lbit::ugate_type));
            auto lbitSA = qm::lbit::get_lbit_support_analysis(urange, frange, sample, tensors.get_length(), ugate_props, false);

            h5file->writeDataset(lbitSA.cls_avg, "/fLBIT/model/lbits/cls_avg");
            h5file->writeDataset(lbitSA.sse_avg, "/fLBIT/model/lbits/sse_avg");
            h5file->writeDataset(lbitSA.decay, "/fLBIT/model/lbits/decay");
            h5file->writeAttribute(urange, "/fLBIT/model/lbits", "u_layer");
            h5file->writeAttribute(frange, "/fLBIT/model/lbits", "f_mixer");
            h5file->writeAttribute(sample, "/fLBIT/model/lbits", "samples");
            if(settings::storage::storage_level_model != StorageLevel::NONE) {
                h5file->writeDataset(lbitSA.supps, "/fLBIT/model/lbits/data", H5D_CHUNKED);
                h5file->writeDataset(lbitSA.pupps, "/fLBIT/model/lbits/pata", H5D_CHUNKED);
                h5file->writeAttribute("The operator support matrix O(i,j) = (1/2^L) Tr(tau_i^z sigma_j^z)", "/fLBIT/model/lbits/data", "description");
                h5file->writeAttribute("The operator support matrix with shifted columns O(i,j) --> O(i,|i-j|)", "/fLBIT/model/lbits/pata", "description");
                h5file->writeAttribute(std::vector<std::string>{"f_mixer", "u_layer", "sample", "i", "j"}, "/fLBIT/model/lbits/data", "dims");
            }
        }
    }
}

void flbit::print_status() {
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
    //    report += fmt::format("E/L:{:<20.16f} ", tools::finite::measure::energy_per_site(tensors));
    report += fmt::format("ε:{:<8.2e} ", tensors.state->get_truncation_error_midchain());
    report += fmt::format("Sₑ(L/2):{:<10.8f} ", tools::finite::measure::entanglement_entropy_midchain(*tensors.state));
    if(tensors.state->measurements.number_entropy_midchain) // This one is expensive
        report += fmt::format("Sₙ(L/2):{:<10.8f} ", tensors.state->measurements.number_entropy_midchain.value());

    if(state_lbit->measurements.number_entropy_midchain) // This one is expensive
        report += fmt::format("Sₙ(L/2) lbit:{:<10.8f} ", state_lbit->measurements.number_entropy_midchain.value());

    report += fmt::format("χ:{:<3}|{:<3}|{:<3} ", settings::get_bond_max(status.algo_type), status.bond_lim,
                          tools::finite::measure::bond_dimension_midchain(*tensors.state));
    report += fmt::format("wtime:{:<} ", fmt::format("{:>6.2f}s", tid::get_unscoped("t_tot").get_time()));
    report += fmt::format("ptime:{:<} ", fmt::format("{:+>8.2e}s", status.phys_time));
    report += fmt::format("itime:{:<} ", fmt::format("{:>4.2f}s", tid::get_unscoped("fLBIT")["run"].get_lap()));
    report += fmt::format("mem[rss {:<.1f}|peak {:<.1f}|vm {:<.1f}]MB ", debug::mem_rss_in_mb(), debug::mem_hwm_in_mb(), debug::mem_vm_in_mb());
    tools::log->info(report);
}