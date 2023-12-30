#include "flbit.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "debug/info.h"
#include "general/iter.h"
#include "io/fmt.h"
#include "math/float.h"
#include "math/linalg/tensor.h"
#include "math/num.h"
#include "math/tenx.h"
#include "qm/lbit.h"
#include "qm/spin.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include "tools/common/prof.h"
#include "tools/common/split.h"
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
        tools::finite::h5::load::simulation(*h5file, state_prefix, tensors, status, status.algo_type);

        // Our first task is to decide on a state name for the newly loaded state
        // The simplest is to inferr it from the state prefix itself
        auto name = tools::common::h5::resume::extract_state_name(state_prefix);
        tensors.state->set_name(name);
        clear_convergence_status();

        // Create an initial state in the real basis
        auto bond_lim = settings::get_bond_init(status.algo_type);
        switch(settings::strategy::initial_state) {
            case StateInit::PRODUCT_STATE_DOMAIN_WALL:
            case StateInit::PRODUCT_STATE_PATTERN:
            case StateInit::PRODUCT_STATE_NEEL:
            case StateInit::PRODUCT_STATE_NEEL_SHUFFLED:
            case StateInit::PRODUCT_STATE_NEEL_DISLOCATED: break;
            default:
                tools::log->warn("Expected initial_state: "
                                 "PRODUCT_STATE_DOMAIN_WALL|"
                                 "PRODUCT_STATE_PATTERN|"
                                 "PRODUCT_STATE_NEEL|"
                                 "PRODUCT_STATE_NEEL_SHUFFLED|"
                                 "PRODUCT_STATE_NEEL_DISLOCATED"
                                 ". Got {}",
                                 enum2sv(settings::strategy::initial_state));
        }

        if(settings::strategy::initial_axis.find("z") == std::string::npos)
            tools::log->warn("Expected initial_axis == z. Got {}", settings::strategy::initial_axis);

        tensors.initialize_state(ResetReason::INIT, settings::strategy::initial_state, StateInitType::REAL, settings::strategy::initial_axis,
                                 settings::strategy::use_eigenspinors, bond_lim, settings::strategy::initial_pattern);

        tensors.move_center_point_to_inward_edge();

        // Load the unitary circuit
        unitary_gates_2site_layers = qm::lbit::read_unitary_2site_gate_layers(*h5file, "/fLBIT/model/unitary_circuit");

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
    using namespace settings::strategy;
    while(not task_list.empty()) {
        auto task = task_list.front();
        switch(task) {
            case flbit_task::INIT_RANDOMIZE_MODEL: initialize_model(); break;
            case flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE: initialize_state(ResetReason::INIT, StateInit::RANDOM_PRODUCT_STATE); break;
            case flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE_NEEL_SHUFFLED:
                initialize_state(ResetReason::INIT, StateInit::PRODUCT_STATE_NEEL_SHUFFLED);
                break;
            case flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE_NEEL_DISLOCATED:
                initialize_state(ResetReason::INIT, StateInit::PRODUCT_STATE_NEEL_DISLOCATED);
                break;
            case flbit_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE: initialize_state(ResetReason::INIT, StateInit::RANDOM_ENTANGLED_STATE); break;
            case flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE_PATTERN: initialize_state(ResetReason::INIT, StateInit::PRODUCT_STATE_PATTERN); break;
            case flbit_task::INIT_RANDOMIZE_INTO_MIDCHAIN_SINGLET_NEEL_STATE:
                initialize_state(ResetReason::INIT, StateInit::MIDCHAIN_SINGLET_NEEL_STATE);
                break;
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
                create_unitary_circuit_gates();
                create_hamiltonian_gates();
                update_time_evolution_gates();
                break;
            }
            case flbit_task::TIME_EVOLVE:
                tensors.state->set_name("state_real");
                run_algorithm();
                break;
            case flbit_task::TRANSFORM_TO_LBIT: transform_to_lbit_basis(); break;
            case flbit_task::TRANSFORM_TO_REAL: transform_to_real_basis(); break;
            case flbit_task::POST_WRITE_RESULT: write_to_file(StorageEvent::FINISHED); break;
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
    initialize_model(); // First use of random!
    init_bond_dimension_limits();
    init_truncation_error_limits();

    // Create an initial state in the real basis
    tensors.state->set_name("state_real");
    initialize_state(ResetReason::INIT, settings::strategy::initial_state);
    tensors.move_center_point_to_inward_edge();
    state_real_init = std::make_unique<StateFinite>(*tensors.state);
    tools::finite::print::model(*tensors.model);

    create_unitary_circuit_gates();

    // Store the model
    write_to_file(StorageEvent::MODEL, CopyPolicy::TRY);

    create_time_points();
    update_time_step();
    create_hamiltonian_gates();
    update_time_evolution_gates();

    if(not tensors.position_is_inward_edge()) throw except::logic_error("Put the state on an edge!");

    // Generate the corresponding state in lbit basis
    transform_to_lbit_basis();
    tools::log->info("Finished {} preprocessing", status.algo_type_sv());
}

void flbit::run_algorithm() {
    if(tensors.state->get_name().empty()) tensors.state->set_name("state_real");
    if(not state_lbit) transform_to_lbit_basis();
    if(time_points.empty()) create_time_points();
    if(cmp_t(status.delta_t.to_floating_point<cplx_t>(), 0.0)) update_time_step();
    tools::log->info("Starting {} algorithm with model [{}] for state [{}]", status.algo_type_sv(), enum2sv(settings::model::model_type),
                     tensors.state->get_name());
    auto t_run = tid::tic_scope("run");
    if(not tensors.position_is_inward_edge()) throw except::logic_error("Put the state on an edge!");
    while(true) {
        update_state();
        check_convergence();
        write_to_file(StorageEvent::ITERATION);
        print_status();
        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);

        // It's important not to perform the last move, so we break now: that last state would not get optimized
        if(status.algo_stop != AlgorithmStop::NONE) break;
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
    auto delta_t = status.delta_t.to_floating_point<cplx_t>();
    tools::log->debug("Starting fLBIT: iter {} | Δt = ({:.2e}, {:.2e})", status.iter, f128_t(std::real(delta_t)), f128_t(std::imag(delta_t)));
    if(not state_lbit) throw except::logic_error("state_lbit == nullptr: Set the state in lbit basis before running an flbit step");
    if(not state_lbit_init) {
        state_lbit_init = std::make_unique<StateFinite>(*state_lbit);
        tools::finite::mps::normalize_state(*state_lbit_init, std::nullopt, NormPolicy::ALWAYS);
    }
    *state_lbit = *state_lbit_init;

    // Time evolve from 0 to time_point[iter] here
    auto t_step = tid::tic_scope("step");
    time_evolve_lbit_state();
    transform_to_real_basis();

    tensors.clear_measurements();
    tensors.clear_cache();
    status.phys_time = abs_t(time_points[std::min(status.iter, time_points.size() - 1)]);

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
    if(settings::flbit::time_scale == TimeScale::LOGSPACED)
        if(cmp_t(status.delta_t.to_floating_point<cplx_t>(), 0.0)) throw except::logic_error("Expected nonzero delta_t after time step update");
    tools::log->debug("Time step iter {} | Δt = {} | t = {:8.2e}", status.iter, status.delta_t.to_floating_point<cplx>(),
                      status.phys_time.to_floating_point<real>());
}

void flbit::check_convergence() {
    if(not tensors.position_is_inward_edge()) return;
    auto t_con = tid::tic_scope("check_conv");
    //    check_convergence_entg_entropy();
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
    auto   t_crt = tid::tic_scope("create_time_points");
    cplx_t time_start(settings::flbit::time_start_real, settings::flbit::time_start_imag);
    cplx_t time_final(settings::flbit::time_final_real, settings::flbit::time_final_imag);
    tools::log->info("Creating time points ({},{}) -> ({},{})", settings::flbit::time_start_real, settings::flbit::time_start_imag,
                     settings::flbit::time_final_real, settings::flbit::time_final_imag);

    cplx_t time_diff = time_start - time_final;
    // Check that there will be some time evolution
    if(cmp_t(time_diff, 0.0)) throw except::logic_error("time_start - time_final == 0");
    // Check that the time limits are purely real or imaginary!
    bool time_is_real = abs_t(time_diff.real()) > 0;
    bool time_is_imag = abs_t(time_diff.imag()) > 0;
    if(time_is_real and time_is_imag)
        throw except::logic_error("time_start and time_final must both be either purely real or imaginary. Got:\n"
                                  "time_start = {:.8f}{:+.8f}\n"
                                  "time_final = {:.8f}{:+.8f}",
                                  f128_t(time_start.real()), f128_t(time_start.imag()), f128_t(time_final.real()), f128_t(time_final.imag()));
    time_points.reserve(settings::flbit::time_num_steps);
    if(settings::flbit::time_scale == TimeScale::LOGSPACED) {
        if(time_is_real) {
            for(const auto &t : num::LogSpaced(settings::flbit::time_num_steps, time_start.real(), time_final.real())) { time_points.emplace_back(t); }
        } else {
            for(const auto &t : num::LogSpaced(settings::flbit::time_num_steps, time_start.imag(), time_final.imag())) {
                time_points.emplace_back(cplx_t(0.0, t));
            }
        }
    } else if(settings::flbit::time_scale == TimeScale::LINSPACED) {
        if(time_is_real) {
            for(const auto &t : num::LinSpaced(settings::flbit::time_num_steps, time_start.real(), time_final.real())) { time_points.emplace_back(t); }
        } else {
            for(const auto &t : num::LinSpaced(settings::flbit::time_num_steps, time_start.imag(), time_final.imag())) {
                time_points.emplace_back(cplx_t(0.0, t));
            }
        }
    }

    //    tools::log->debug("Created {} time points:\n{}", time_points.size(), time_points);
    // Sanity check
    if(time_points.front().real() != settings::flbit::time_start_real) throw except::logic_error("Time start real mismatch");
    if(time_points.front().imag() != settings::flbit::time_start_imag) throw except::logic_error("Time start imag mismatch");
    if(time_points.back().real() != settings::flbit::time_final_real) throw except::logic_error("Time final real mismatch");
    if(time_points.back().imag() != settings::flbit::time_final_imag) throw except::logic_error("Time final imag mismatch");
    if(time_points.size() != settings::flbit::time_num_steps)
        throw except::logic_error("Got time_points.size():[{}] != settings::flbit::time_num_steps:[{}]", time_points.size(), settings::flbit::time_num_steps);
}

void flbit::create_hamiltonian_gates() {
    // Create the hamiltonian gates with n-site terms
    //
    // Note for 2-body terms
    // We want MPO's with only 2-body terms of the Hamiltonian, for distances |i-j| <= J2_span sites.
    // Let's assume we have L==8 and J2_span == 3, then we apply time-evolution
    // L     :  0,1,2,3,4,5,6,7
    // mpo[0]: [0,1,2,3]
    // mpo[1]:   [1,2,3,4]
    // mpo[2]:     [2,3,4,5]
    // mpo[3]:       [3,4,5,6]
    // mpo[4]:         [4,5,6,7]
    // Note that the interaction on sites {1,2} gets applied twice, and {3,4} thrice.
    // To compensate for this, pass a "0" to nbody, which tells the mpo-generator to divide J[i,j] by
    // the number of times it is applied. This ensures that each interaction is time-evolved the right number of times.

    bool has_swap_gates = not ham_swap_gates_1body.empty() or not ham_swap_gates_2body.empty() or not ham_swap_gates_3body.empty();
    bool has_slow_gates = not ham_gates_1body.empty() or not ham_gates_2body.empty() or not ham_gates_3body.empty();

    if(has_swap_gates) tools::log->warn("Hamiltonian swap gates have already been constructed: Normally this is done just once.");
    if(has_slow_gates) tools::log->warn("Hamiltonian gates have already been constructed: Normally this is done just once.");

    ham_gates_1body.clear();
    ham_gates_2body.clear();
    ham_gates_3body.clear();
    ham_gates_Lbody.clear();
    ham_swap_gates_1body.clear();
    ham_swap_gates_2body.clear();
    ham_swap_gates_3body.clear();
    ham_swap_gates_Lbody.clear();

    auto L          = settings::model::model_size;
    auto list_1body = num::range<size_t>(0, L - 0, 1);
    auto list_2body = num::range<size_t>(0, L - 1, 1);
    auto list_3body = num::range<size_t>(0, L - 2, 1);
    if(settings::flbit::use_swap_gates) {
        tools::log->info("Creating Hamiltonian swap gates");
        auto t_swaphamgates = tid::tic_scope("create_ham_gates_swap");
        for(auto pos : list_1body) {
            auto range = 0ul;                                      // Number of site indices
            auto sites = num::range<size_t>(pos, pos + range + 1); // A list of site indices. +1 to include last site (in this case the only site)
            auto nbody = std::vector<size_t>{1};                   // A list of included nbody interaction terms (1: on-site terms, 2: pairwise, and so on)
            auto spins = tensors.state->get_spin_dims(sites);      // A list of spin dimensions for each site (should all be 2 for two-level systems)
            tools::log->debug("Generating {}-body hamiltonian on sites {}", nbody, sites);
            ham_swap_gates_1body.emplace_back(tensors.model->get_multisite_ham_t({pos}, nbody), sites, spins);
        }
        auto J2_ctof = std::min(settings::model::lbit::J2_span, L - 1); // Max distance |i-j| to the furthest interacting site L-1
        for(auto posL : list_2body) {
            auto maxR = std::min<size_t>(posL + J2_ctof, L - 1);
            if(maxR == posL) continue;
            for(auto posR : num::range<size_t>(posL, maxR + 1)) { // maxR+1 to include the last site in range
                if(posL == posR) continue;
                if(posL >= L) throw except::logic_error("posL {} >= L {}", posL, L);
                if(posR >= L) throw except::logic_error("posR {} >= L {}", posR, L);
                auto sites = std::vector<size_t>{posL, posR}; // +1 to include last site
                auto nbody = std::vector<size_t>{2};
                auto spins = tensors.state->get_spin_dims(sites);
                tools::log->debug("Generating {}-body hamiltonian on sites {}", nbody, sites);
                // Accept all swap gates even if all elements are near zero on gates for remote sites,
                // since at large times t these can become relevant again by exp(-itH)
                ham_swap_gates_2body.emplace_back(tensors.model->get_multisite_ham_t(sites, nbody), sites, spins);
            }
            tensors.model->clear_cache();
        }
        // Ignore Hamiltonians with entries smaller than J2_zero: the timescale is too small to resolve them.

        for(auto pos : list_3body) {
            auto range = 2ul;                                                                  // Distance to next-nearest neighbor when 3 sites interact
            auto sites = num::range<size_t>(pos, std::clamp<size_t>(pos + range + 1, 0ul, L)); // +1 to include last site
            auto nbody = std::vector<size_t>{3};
            auto spins = tensors.state->get_spin_dims(sites);
            tools::log->debug("Generating {}-body hamiltonian on sites {}", nbody, sites);
            ham_swap_gates_3body.emplace_back(tensors.model->get_multisite_ham_t(sites, nbody), sites, spins);
        }
        tensors.model->clear_cache();

        if(L <= 6) { // Used for test/debug on small systems
            auto list_Lbody = num::range<size_t>(0, L - 0, 1);
            auto nbody      = std::vector<size_t>{1, 2, 3};
            auto spins      = tensors.state->get_spin_dims(list_Lbody);
            ham_swap_gates_Lbody.emplace_back(tensors.model->get_multisite_ham_t(list_Lbody, nbody), list_Lbody, spins);
        }

        for(const auto &[idx, ham] : iter::enumerate(ham_swap_gates_1body)) {
            if(tenx::isZero(ham.op_t)) tools::log->warn("ham1 swap [{}] is all zeros", idx);
        }
        for(const auto &[idx, ham] : iter::enumerate(ham_swap_gates_2body)) {
            if(tenx::isZero(ham.op_t)) tools::log->warn("ham2 swap [{}] for sites {} is a zero-matrix", idx, ham.pos);
        }
        for(const auto &[idx, ham] : iter::enumerate(ham_swap_gates_3body))
            if(tenx::isZero(ham.op_t)) tools::log->warn("ham3 swap [{}] is all zeros", idx);
    } else {
        tools::log->info("Creating Hamiltonian gates");
        auto t_hamgates = tid::tic_scope("create_ham_gates");
        for(auto pos : list_1body) {
            auto sites = std::vector<size_t>{pos};            // A list of site indices
            auto nbody = std::vector<size_t>{1};              // A list of included nbody interaction terms (1: on-site terms, 2: pairwise, and so on)
            auto spins = tensors.state->get_spin_dims(sites); // A list of spin dimensions for each site (should all be 2 for two-level systems)
            tools::log->debug("Generating {}-body hamiltonian on sites {}", nbody, sites);
            ham_gates_1body.emplace_back(tensors.model->get_multisite_ham_t({pos}, nbody), sites, spins);
        }
        tensors.model->clear_cache();

        auto J2_ctof = std::min(settings::model::lbit::J2_span, L - 1); // Max distance |i-j| to the furthest interacting site L-1
        for(auto posL : list_2body) {
            auto posR = posL + J2_ctof;
            if(J2_ctof == 0) break;
            if(posR >= L) break;
            auto sites = num::range<size_t>(posL, posR + 1); // +1 to include last site, which is otherwise not included in range
            auto nbody = std::vector<size_t>{0, 2};          // zero is a flag to enable compensation for double-counting
            auto spins = tensors.state->get_spin_dims(sites);
            tools::log->debug("Generating {}-body hamiltonian on sites {}", nbody, sites);
            ham_gates_2body.emplace_back(tensors.model->get_multisite_ham_t(sites, nbody), sites, spins);
        }
        tensors.model->clear_cache();

        for(auto posL : list_3body) {
            auto range = 2ul;                                                                    // Distance to next-nearest neighbor when 3 sites interact
            auto sites = num::range<size_t>(posL, std::clamp<size_t>(posL + range + 1, 0ul, L)); // +1 to include last site
            auto nbody = std::vector<size_t>{3};
            auto spins = tensors.state->get_spin_dims(sites);
            tools::log->debug("Generating {}-body hamiltonian on sites {}", nbody, sites);
            ham_gates_3body.emplace_back(tensors.model->get_multisite_ham_t(sites, nbody), sites, spins);
        }
        tensors.model->clear_cache();

        if(L <= 6) { // Used for test/debug on small systems
            auto list_Lbody = num::range<size_t>(0, L - 0, 1);
            auto nbody      = std::vector<size_t>{1, 2, 3};
            auto spins      = tensors.state->get_spin_dims(list_Lbody);
            ham_gates_Lbody.emplace_back(tensors.model->get_multisite_ham_t(list_Lbody, nbody), list_Lbody, spins);
        }
        tensors.model->clear_cache();
        for(const auto &[idx, ham] : iter::enumerate(ham_gates_1body))
            if(tenx::isZero(ham.op_t)) tools::log->warn("ham1[{}] is all zeros", idx);

        for(const auto &[idx, ham] : iter::enumerate(ham_gates_2body))
            if(tenx::isZero(ham.op_t)) tools::log->warn("ham2[{}] is all zeros", idx);

        for(const auto &[idx, ham] : iter::enumerate(ham_gates_3body))
            if(tenx::isZero(ham.op_t)) tools::log->warn("ham3[{}] is all zeros", idx);
    }
    tensors.model->clear_cache();
}

void flbit::update_time_evolution_gates() {
    // Create the time evolution operators
    if(time_points.empty()) create_time_points();
    update_time_step();
    bool has_swap_gates = not ham_swap_gates_1body.empty() or not ham_swap_gates_2body.empty() or not ham_swap_gates_3body.empty();
    bool has_slow_gates = not ham_gates_1body.empty() or not ham_gates_2body.empty() or not ham_gates_3body.empty();

    if(not has_swap_gates and not has_slow_gates) throw except::logic_error("Hamiltonian gates have not been constructed");
    if(has_swap_gates and has_slow_gates) tools::log->warn("Both swap/non-swap gates have been constructed: Normally only one type should be used");
    auto delta_t = status.delta_t.to_floating_point<cplx_t>();
    if(has_swap_gates) {
        auto t_upd = tid::tic_scope("upd_time_evo_swap_gates");
        tools::log->debug("Updating time evolution swap gates to iter {} | Δt = ({:.2e}, {:.2e})", status.iter, f128_t(std::real(delta_t)),
                          f128_t(std::imag(delta_t)));
        time_swap_gates_1body = qm::lbit::get_time_evolution_swap_gates(delta_t, ham_swap_gates_1body, settings::flbit::time_gate_id_threshold);
        time_swap_gates_2body = qm::lbit::get_time_evolution_swap_gates(delta_t, ham_swap_gates_2body, settings::flbit::time_gate_id_threshold);
        time_swap_gates_3body = qm::lbit::get_time_evolution_swap_gates(delta_t, ham_swap_gates_3body, settings::flbit::time_gate_id_threshold);
        if(settings::model::model_size <= 6)
            time_swap_gates_Lbody = qm::lbit::get_time_evolution_swap_gates(delta_t, ham_swap_gates_Lbody, settings::flbit::time_gate_id_threshold);
    }
    if(has_slow_gates) {
        auto t_upd = tid::tic_scope("upd_time_evo_gates");
        tools::log->debug("Updating time evolution gates to iter {} | Δt = ({:.2e}, {:.2e})", status.iter, f128_t(std::real(delta_t)),
                          f128_t(std::imag(delta_t)));
        time_gates_1body = qm::lbit::get_time_evolution_gates(delta_t, ham_gates_1body, settings::flbit::time_gate_id_threshold);
        time_gates_2body = qm::lbit::get_time_evolution_gates(delta_t, ham_gates_2body, settings::flbit::time_gate_id_threshold);
        time_gates_3body = qm::lbit::get_time_evolution_gates(delta_t, ham_gates_3body, settings::flbit::time_gate_id_threshold);
        if(settings::model::model_size <= 6)
            time_gates_Lbody = qm::lbit::get_time_evolution_gates(delta_t, ham_gates_Lbody, settings::flbit::time_gate_id_threshold);
    }
}

void flbit::create_unitary_circuit_gates() {
    if(unitary_gates_2site_layers.size() == settings::model::lbit::u_depth) return;
    std::vector<double> fields;
    for(const auto &field : tensors.model->get_parameter("J1_rand")) fields.emplace_back(static_cast<double>(std::any_cast<real_t>(field)));
    unitary_gates_2site_layers.resize(settings::model::lbit::u_depth);
    uprop              = qm::lbit::UnitaryGateProperties(fields);
    uprop.keep_circuit = settings::storage::table::random_unitary_circuit::policy != StoragePolicy::NONE;
    tools::log->info("Creating unitary circuit of 2-site gates {}", uprop.string());
    for(auto &ulayer : unitary_gates_2site_layers) ulayer = qm::lbit::create_unitary_2site_gate_layer(uprop);
    if(settings::flbit::use_mpo_circuit) {
        unitary_gates_mpo_layers.clear();
        for(const auto &ulayer : unitary_gates_2site_layers) unitary_gates_mpo_layers.emplace_back(qm::lbit::get_unitary_mpo_layer(ulayer));
        //        unitary_gates_mpo_layer_full = qm::lbit::merge_unitary_mpo_layers(mpo_layers);
        ledge.resize(1);
        redge.resize(1);
        ledge.setConstant(cplx(1.0, 0.0));
        redge.setConstant(cplx(1.0, 0.0));
    }
}

void flbit::time_evolve_lbit_state() {
    auto t_evo          = tid::tic_scope("time_evo");
    auto svd_cfg        = svd::config(status.bond_lim, status.trnc_lim);
    svd_cfg.svd_lib     = svd::lib::lapacke;
    svd_cfg.svd_rtn     = svd::rtn::geauto;
    bool has_swap_gates = not time_swap_gates_1body.empty() or not time_swap_gates_2body.empty() or not time_swap_gates_3body.empty();
    bool has_slow_gates = not time_gates_1body.empty() or not time_gates_2body.empty() or not time_gates_3body.empty();
    if(has_swap_gates and has_slow_gates) throw except::logic_error("Both swap and non-swap time evolution gates found");
    if(not has_swap_gates and not has_slow_gates) throw except::logic_error("None of swap or non-swap time evolution gates found");
    auto delta_t = status.delta_t.to_floating_point<cplx_t>();
    if(has_swap_gates) {
        tools::log->debug("Applying time evolution swap gates Δt = ({:.2e}, {:.2e})", f128_t(std::real(delta_t)), f128_t(std::imag(delta_t)));
        tools::finite::mps::apply_swap_gates(*state_lbit, time_swap_gates_1body, CircuitOp::NONE, GateMove::AUTO, svd_cfg);
        tools::finite::mps::apply_swap_gates(*state_lbit, time_swap_gates_2body, CircuitOp::NONE, GateMove::AUTO, svd_cfg);
        tools::finite::mps::apply_swap_gates(*state_lbit, time_swap_gates_3body, CircuitOp::NONE, GateMove::AUTO, svd_cfg);
    }
    if(has_slow_gates) {
        tools::log->debug("Applying time evolution gates Δt = ({:.2e}, {:.2e})", f128_t(std::real(delta_t)), f128_t(std::imag(delta_t)));
        tools::finite::mps::apply_gates(*state_lbit, time_gates_1body, CircuitOp::NONE, true, GateMove::AUTO, svd_cfg);
        tools::finite::mps::apply_gates(*state_lbit, time_gates_2body, CircuitOp::NONE, true, GateMove::AUTO, svd_cfg);
        tools::finite::mps::apply_gates(*state_lbit, time_gates_3body, CircuitOp::NONE, true, GateMove::AUTO, svd_cfg);
    }
    tools::finite::mps::normalize_state(*state_lbit, std::nullopt, NormPolicy::IFNEEDED);

    t_evo.toc();

    if constexpr(settings::debug) {
        // Check that we would get back the original state if we time evolved backwards
        auto state_lbit_debug = *state_lbit;
        if(has_swap_gates) {
            tools::log->debug("Applying time evolution swap gates backward Δt = ({:.2e}, {:.2e})", f128_t(std::real(delta_t)), f128_t(std::imag(delta_t)));
            for(auto &g : time_swap_gates_1body) g.unmark_as_used();
            for(auto &g : time_swap_gates_2body) g.unmark_as_used();
            for(auto &g : time_swap_gates_3body) g.unmark_as_used();
            tools::finite::mps::apply_swap_gates(state_lbit_debug, time_swap_gates_3body, CircuitOp::ADJ, GateMove::AUTO, svd_cfg);
            tools::finite::mps::apply_swap_gates(state_lbit_debug, time_swap_gates_2body, CircuitOp::ADJ, GateMove::AUTO, svd_cfg);
            tools::finite::mps::apply_swap_gates(state_lbit_debug, time_swap_gates_1body, CircuitOp::ADJ, GateMove::AUTO, svd_cfg);
        }
        if(has_slow_gates) {
            tools::log->debug("Applying time evolution gates backward Δt = ({:.2e}, {:.2e})", f128_t(std::real(delta_t)), f128_t(std::imag(delta_t)));
            for(auto &g : time_gates_1body) g.unmark_as_used();
            for(auto &g : time_gates_2body) g.unmark_as_used();
            for(auto &g : time_gates_3body) g.unmark_as_used();
            tools::finite::mps::apply_gates(state_lbit_debug, time_gates_3body, CircuitOp::ADJ, true, GateMove::AUTO, svd_cfg);
            tools::finite::mps::apply_gates(state_lbit_debug, time_gates_2body, CircuitOp::ADJ, true, GateMove::AUTO, svd_cfg);
            tools::finite::mps::apply_gates(state_lbit_debug, time_gates_1body, CircuitOp::ADJ, true, GateMove::AUTO, svd_cfg);
        }
        tools::finite::mps::normalize_state(state_lbit_debug, std::nullopt, NormPolicy::IFNEEDED);
        auto overlap = tools::finite::ops::overlap(*state_lbit_init, state_lbit_debug);
        tools::log->info("Debug overlap after time evolution: {:.16f}", overlap);
        if(std::abs(overlap - 1.0) > 10 * status.trnc_lim)
            throw except::runtime_error("State overlap after backwards time evolution is not 1: Got {:.16f}", overlap);
    }
}

void flbit::transform_to_real_basis() {
    if(unitary_gates_2site_layers.size() != settings::model::lbit::u_depth) create_unitary_circuit_gates();
    auto t_map    = tid::tic_scope("l2r");
    tensors.state = std::make_unique<StateFinite>(*state_lbit);
    tensors.state->set_name("state_real");
    tensors.state->clear_cache();
    tensors.state->clear_measurements();
    tensors.clear_measurements();
    tensors.clear_cache();

    auto svd_cfg    = svd::config(status.bond_lim, status.trnc_lim);
    svd_cfg.svd_lib = svd::lib::lapacke;
    svd_cfg.svd_rtn = svd::rtn::geauto;
    if(settings::flbit::use_mpo_circuit) {
        svd_cfg.rank_max = static_cast<long>(static_cast<double>(tensors.state->find_largest_bond()) * 4);
        tools::log->debug("Transforming {} to {} using {} unitary mpo layers", state_lbit->get_name(), tensors.state->get_name(),
                          unitary_gates_mpo_layers.size());
        for(const auto &[idx_layer, mpo_layer] : iter::enumerate(unitary_gates_mpo_layers)) {
            tools::finite::ops::apply_mpos(*tensors.state, mpo_layer, ledge, redge, true);
            if((idx_layer + 1) % 1 == 0) {
                tools::finite::mps::normalize_state(*tensors.state, svd_cfg, NormPolicy::ALWAYS);
                svd_cfg.rank_max = static_cast<long>(static_cast<double>(tensors.state->find_largest_bond()) * 4);
            }
        }
    } else {
        tools::log->debug("Transforming {} to {} using {} unitary layers", state_lbit->get_name(), tensors.state->get_name(),
                          unitary_gates_2site_layers.size());
        tools::finite::mps::apply_circuit(*tensors.state, unitary_gates_2site_layers, CircuitOp::ADJ, false, true, GateMove::ON, svd_cfg);
    }

    tools::finite::mps::normalize_state(*tensors.state, svd_cfg, NormPolicy::IFNEEDED);
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
        // Check that the transform backwards is equal to the original state
        if(settings::flbit::use_mpo_circuit) {
            auto state_lbit_debug = *tensors.state;
            for(const auto &[idx_layer, mpo_layer] : iter::enumerate(unitary_gates_mpo_layers)) {
                tools::finite::ops::apply_mpos(state_lbit_debug, mpo_layer, ledge, redge, false);
                if((idx_layer + 1) % 1 == 0) {
                    tools::finite::mps::normalize_state(state_lbit_debug, svd_cfg, NormPolicy::ALWAYS);
                    svd_cfg.rank_max = static_cast<long>(static_cast<double>(tensors.state->find_largest_bond()) * 2);
                }
            }
            tools::finite::mps::normalize_state(state_lbit_debug, std::nullopt, NormPolicy::IFNEEDED);
            auto overlap = tools::finite::ops::overlap(*state_lbit, state_lbit_debug);
            tools::log->info("Debug overlap after unitary circuit: {:.16f}", overlap);
            if(std::abs(overlap - 1.0) > 10 * status.trnc_lim)
                throw except::runtime_error("State overlap after transform back from real is not 1: Got {:.16f}", overlap);
        }
        {
            auto state_lbit_debug = *tensors.state;
            tools::finite::mps::apply_circuit(state_lbit_debug, unitary_gates_2site_layers, CircuitOp::NONE, false, true, GateMove::ON, svd_cfg);
            auto overlap = tools::finite::ops::overlap(*state_lbit, state_lbit_debug);
            tools::log->info("Debug overlap after unitary circuit: {:.16f}", overlap);
            if(std::abs(overlap - 1.0) > 10 * status.trnc_lim)
                throw except::runtime_error("State overlap after transform back from real is not 1: Got {:.16f}", overlap);
        }
    }
}

void flbit::transform_to_lbit_basis() {
    if(unitary_gates_2site_layers.size() != settings::model::lbit::u_depth) create_unitary_circuit_gates();
    auto t_map = tid::tic_scope("r2l");
    state_lbit = std::make_unique<StateFinite>(*tensors.state);
    state_lbit->set_name("state_lbit");
    state_lbit->clear_cache();
    state_lbit->clear_measurements();
    auto svd_cfg    = svd::config(status.bond_lim, status.trnc_lim);
    svd_cfg.svd_lib = svd::lib::lapacke;
    svd_cfg.svd_rtn = svd::rtn::geauto;
    if(settings::flbit::use_mpo_circuit) {
        svd_cfg.rank_max = static_cast<long>(static_cast<double>(state_lbit->find_largest_bond()) * 4);
        tools::log->info("Transforming {} to {} using {} unitary mpo layers", tensors.state->get_name(), state_lbit->get_name(),
                         unitary_gates_mpo_layers.size());
        for(const auto &[idx_layer, mpo_layer] : iter::enumerate_reverse(unitary_gates_mpo_layers)) {
            tools::finite::ops::apply_mpos(*state_lbit, mpo_layer, ledge, redge, false);
            if((idx_layer) % 1 == 0) {
                tools::log->info("Normalizing with rank_max {} | max bond {}", svd_cfg.rank_max.value(), state_lbit->find_largest_bond());
                tools::finite::mps::normalize_state(*state_lbit, svd_cfg, NormPolicy::ALWAYS);
                svd_cfg.rank_max = static_cast<long>(static_cast<double>(state_lbit->find_largest_bond()) * 4);
            }
        }
    } else {
        tools::log->info("Transforming {} to {} using {} unitary layers", tensors.state->get_name(), state_lbit->get_name(), unitary_gates_2site_layers.size());
        tools::finite::mps::apply_circuit(*state_lbit, unitary_gates_2site_layers, CircuitOp::NONE, false, true, GateMove::ON, svd_cfg);
    }

    //    auto svd_cfg = svd::config(status.bond_lim, status.trnc_lim);
    tools::finite::mps::normalize_state(*state_lbit, std::nullopt, NormPolicy::IFNEEDED);
    status.position  = tensors.get_position<long>();
    status.direction = tensors.state->get_direction();
    tools::log->debug("time r2l: {:.3e} s", t_map->get_last_interval());
    if constexpr(settings::debug) {
        auto t_dbg = tid::tic_scope("debug");
        if(tools::log->level() <= spdlog::level::debug) {
            tools::log->debug("{} bond dimensions: {}", state_lbit->get_name(), tools::finite::measure::bond_dimensions(*state_lbit));
            tools::log->debug("{} bond dimensions: {}", tensors.state->get_name(), tools::finite::measure::bond_dimensions(*tensors.state));
        }
        // Check normalization
        for(const auto &mps : state_lbit->mps_sites) mps->assert_normalized();

        // Double-check the that transform operation backwards is equal to the original state
        if(settings::flbit::use_mpo_circuit) {
            auto state_real_debug = *state_lbit;
            for(const auto &[idx_layer, mpo_layer] : iter::enumerate(unitary_gates_mpo_layers)) {
                tools::finite::ops::apply_mpos(state_real_debug, mpo_layer, ledge, redge, true);
                if((idx_layer + 1) % 1 == 0) { tools::finite::mps::normalize_state(state_real_debug, svd_cfg, NormPolicy::ALWAYS); }
            }
            tools::finite::mps::normalize_state(state_real_debug, std::nullopt, NormPolicy::IFNEEDED);
            auto overlap = tools::finite::ops::overlap(*tensors.state, state_real_debug);
            tools::log->info("Debug overlap: {:.16f}", overlap);
            if(std::abs(overlap - 1.0) > 10 * status.trnc_lim)
                throw except::runtime_error("State overlap after transform back from lbit is not 1: Got {:.16f}", overlap);
        } else {
            auto state_real_debug = *state_lbit;
            tools::finite::mps::apply_circuit(state_real_debug, unitary_gates_2site_layers, CircuitOp::ADJ, false, true, GateMove::ON, svd_cfg);
            auto overlap = tools::finite::ops::overlap(*tensors.state, state_real_debug);
            tools::log->info("Debug overlap: {:.16f}", overlap);
            if(std::abs(overlap - 1.0) > 10 * status.trnc_lim)
                throw except::runtime_error("State overlap after transform back from lbit is not 1: Got {:.16f}", overlap);
        }
    }
}

void flbit::write_to_file(StorageEvent storage_event, CopyPolicy copy_policy) {
    AlgorithmFinite::write_to_file(*tensors.state, *tensors.model, *tensors.edges, storage_event, copy_policy);
    if(not state_lbit) transform_to_lbit_basis();
    AlgorithmFinite::write_to_file(*state_lbit, *tensors.model, *tensors.edges, storage_event, copy_policy);

    if(h5file and storage_event == StorageEvent::INIT) {
        using namespace settings::flbit;
        auto state_prefix = fmt::format("/fLBIT/{}", tensors.state->get_name());
        h5file->createGroup(state_prefix);
        h5file->writeAttribute(enum2sv(time_scale), state_prefix, "time_scale");
        h5file->writeAttribute(time_num_steps, state_prefix, "time_num_steps");
        h5file->writeAttribute(std::complex(time_start_real, time_start_imag), state_prefix, "time_start");
        h5file->writeAttribute(std::complex(time_final_real, time_final_imag), state_prefix, "time_final");
    }

    // Save the unitaries once
    if(h5file and storage_event == StorageEvent::MODEL) {
        auto t_h5 = tid::tic_scope("h5");
        qm::lbit::write_unitary_circuit_parameters(*h5file, "/fLBIT/model/unitary_circuit", uprop.circuit); // Writes a table with gate parameters for resuming
        h5file->writeAttribute(uprop.depth, "/fLBIT/model/unitary_circuit", "u_depth");
        h5file->writeAttribute(uprop.fmix, "/fLBIT/model/unitary_circuit", "u_fmix");
        h5file->writeAttribute(uprop.tstd, "/fLBIT/model/unitary_circuit", "u_tstd");
        h5file->writeAttribute(uprop.cstd, "/fLBIT/model/unitary_circuit", "u_cstd");
        h5file->writeAttribute(uprop.g8w8, "/fLBIT/model/unitary_circuit", "u_g8w8");
        h5file->writeAttribute(uprop.type, "/fLBIT/model/unitary_circuit", "u_type");
    }
    //    if(h5file and storage_event == StorageEvent::MODEL) {
    //        auto        t_h5      = tid::tic_scope("h5");
    //        auto        t_event   = tid::tic_scope(enum2sv(storage_event), tid::highest);
    //        auto        t_gates   = tid::tic_scope("gates");
    //        std::string grouppath = "/fLBIT/model/unitary_gates";
    //        if(h5file->linkExists(grouppath)) return;
    //        for(const auto &[idx_layer, layer] : iter::enumerate(unitary_gates_2site_layers)) {
    //            std::string layerpath = fmt::format("{}/layer_{}", grouppath, idx_layer);
    //            for(const auto &[idx_gate, u] : iter::enumerate(layer)) {
    //                std::string gatepath = fmt::format("{}/u_{}", layerpath, idx_gate);
    //                h5file->writeDataset(u.op, gatepath);
    //                h5file->writeAttribute(u.dim, gatepath, "dim");
    //                h5file->writeAttribute(u.pos, gatepath, "pos");
    //            }
    //            h5file->writeAttribute(static_cast<size_t>(idx_layer), layerpath, "idx_layer");
    //            h5file->writeAttribute(layer.size(), layerpath, "num_gates");
    //        }
    //        h5file->writeAttribute(unitary_gates_2site_layers.size(), grouppath, "num_layers");
    //    }

    // Save the lbit analysis once
    if(h5file and storage_event == StorageEvent::MODEL) {
        auto t_h5    = tid::tic_scope("h5");
        auto t_event = tid::tic_scope(enum2sv(storage_event), tid::highest);
        if(h5file and h5file->linkExists("/fLBIT/model/lbits")) return;
        auto nsamps = settings::flbit::cls::num_rnd_circuits;
        if(nsamps > 0) {
            auto                usites = std::vector<size_t>{settings::model::model_size};
            auto                ug8w8s = std::vector<UnitaryGateWeight>{settings::model::lbit::u_g8w8};
            auto                utypes = std::vector<UnitaryGateType>{settings::model::lbit::u_type};
            auto                udpths = std::vector<size_t>{settings::model::lbit::u_depth};
            auto                ufmixs = std::vector<double>{settings::model::lbit::u_fmix};
            auto                utstds = std::vector<double>{settings::model::lbit::u_tstd};
            auto                ucstds = std::vector<double>{settings::model::lbit::u_cstd};
            auto                randhf = settings::flbit::cls::randomize_hfields;
            std::vector<double> fields;
            for(const auto &field : tensors.model->get_parameter("J1_rand")) fields.emplace_back(static_cast<double>(std::any_cast<real_t>(field)));
            auto uprop_default    = qm::lbit::UnitaryGateProperties(fields);
            uprop_default.ulayers = unitary_gates_2site_layers;
            auto lbitSA           = qm::lbit::get_lbit_support_analysis(uprop_default, udpths, ufmixs, utstds, ucstds, ug8w8s);
            if(h5file and settings::storage::dataset::lbit_analysis::level != StorageLevel::NONE) {
                // Put the sample dimension first so that we can collect many simulations in dmrg-meld along the 0'th dim
                auto label_dist = std::vector<std::string>{"sample", "|i-j|"};
                auto shape_avgs = std::vector<long>{1, lbitSA.corravg.size()};
                auto shape_data = std::vector<long>{static_cast<long>(nsamps), static_cast<long>(settings::model::model_size),
                                                    static_cast<long>(settings::model::model_size)};

                auto label_data = std::vector<std::string>{"sample", "i", "j"};
                if(settings::storage::dataset::lbit_analysis::level >= StorageLevel::LIGHT) {
                    h5file->writeDataset(lbitSA.corrmat, "/fLBIT/model/lbits/corrmat", H5D_CHUNKED, shape_data);
                    h5file->writeAttribute(label_data, "/fLBIT/model/lbits/corrmat", "dimensions");
                    h5file->writeAttribute("The operator support matrix O(i,j) = (1/2^L) Tr(tau_i^z sigma_j^z)", "/fLBIT/model/lbits/corrmat", "description");
                }
                if(settings::storage::dataset::lbit_analysis::level == StorageLevel::FULL) {
                    h5file->writeDataset(lbitSA.cls_avg_fit, "/fLBIT/model/lbits/cls_avg_fit");
                    h5file->writeDataset(lbitSA.cls_avg_rms, "/fLBIT/model/lbits/cls_avg_rms");
                    h5file->writeDataset(lbitSA.cls_avg_rsq, "/fLBIT/model/lbits/cls_avg_rsq");
                    h5file->writeDataset(lbitSA.cls_typ_fit, "/fLBIT/model/lbits/cls_typ_fit");
                    h5file->writeDataset(lbitSA.cls_typ_rms, "/fLBIT/model/lbits/cls_typ_rms");
                    h5file->writeDataset(lbitSA.cls_typ_rsq, "/fLBIT/model/lbits/cls_typ_rsq");

                    h5file->writeDataset(lbitSA.corravg, "/fLBIT/model/lbits/corravg", H5D_CHUNKED, shape_avgs);
                    h5file->writeAttribute(label_dist, "/fLBIT/model/lbits/corravg", "dimensions");
                    h5file->writeAttribute("Site arithmetic average <<O(|i-j|)>>", "/fLBIT/model/lbits/corravg", "description");

                    h5file->writeDataset(lbitSA.corrtyp, "/fLBIT/model/lbits/corrtyp", H5D_CHUNKED, shape_avgs);
                    h5file->writeAttribute(label_dist, "/fLBIT/model/lbits/corrtyp", "dimensions");
                    h5file->writeAttribute("Site geometric average <<O(|i-j|)>>_typ", "/fLBIT/model/lbits/corrtyp", "description");

                    h5file->writeDataset(lbitSA.correrr, "/fLBIT/model/lbits/correrr", H5D_CHUNKED, shape_avgs);
                    h5file->writeAttribute(label_dist, "/fLBIT/model/lbits/correrr", "dimensions");
                    h5file->writeAttribute("Standard error of <<O(|i-j|)>>", "/fLBIT/model/lbits/correrr", "description");

                    h5file->writeDataset(lbitSA.corroff, "/fLBIT/model/lbits/corroff", H5D_CHUNKED, shape_data);
                    h5file->writeAttribute(label_data, "/fLBIT/model/lbits/corroff", "dimensions");
                    h5file->writeAttribute("The operator support matrix with shifted columns O(i,j) --> O(i,|i-j|)", "/fLBIT/model/lbits/corroff",
                                           "description");
                }
                h5file->writeAttribute(udpths, "/fLBIT/model/lbits", "u_depth");
                h5file->writeAttribute(ufmixs, "/fLBIT/model/lbits", "u_fmix");
                h5file->writeAttribute(utstds, "/fLBIT/model/lbits", "u_tstd");
                h5file->writeAttribute(ucstds, "/fLBIT/model/lbits", "u_cstd");
                h5file->writeAttribute(enum2sv(ug8w8s), "/fLBIT/model/lbits", "u_g8w8");
                h5file->writeAttribute(enum2sv(utypes), "/fLBIT/model/lbits", "u_type");
                h5file->writeAttribute(nsamps, "/fLBIT/model/lbits", "samples");
                h5file->writeAttribute(randhf, "/fLBIT/model/lbits", "randomize_hfields");
                h5file->writeAttribute(settings::flbit::cls::mpo_circuit_svd_bondlim, "/fLBIT/model/lbits", "u_bond");
                h5file->writeAttribute(settings::flbit::cls::mpo_circuit_svd_trnclim, "/fLBIT/model/lbits", "u_trnc");
                h5file->writeAttribute(settings::flbit::cls::mpo_circuit_switchdepth, "/fLBIT/model/lbits", "mpo_switchdepth");
            }
        }
        if(settings::flbit::cls::exit_when_done) exit(0);
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
    report += fmt::format("Sₑ(L/2):{:<18.16f} ", tools::finite::measure::entanglement_entropy_midchain(*tensors.state));
    if(tensors.state->measurements.number_entropy_midchain) // This one is expensive
        report += fmt::format("Sₙ(L/2):{:<18.16f} ", tensors.state->measurements.number_entropy_midchain.value());

    if(state_lbit->measurements.number_entropy_midchain) // This one is expensive
        report += fmt::format("Sₙ(L/2) lbit:{:<10.8f} ", state_lbit->measurements.number_entropy_midchain.value());

    report += fmt::format("χ:{:<3}|{:<3}|{:<3} ", settings::get_bond_max(status.algo_type), status.bond_lim,
                          tools::finite::measure::bond_dimension_midchain(*tensors.state));
    if(settings::flbit::time_scale == TimeScale::LOGSPACED)
        report += fmt::format("ptime:{:<} ", fmt::format("{:>.2e}s", status.phys_time.to_floating_point<real>()));
    if(settings::flbit::time_scale == TimeScale::LINSPACED)
        report += fmt::format("ptime:{:<} ", fmt::format("{:>.6f}s", status.phys_time.to_floating_point<real>()));
    report += fmt::format("wtime:{:<}(+{:<}) ", fmt::format("{:>.1f}s", tid::get_unscoped("t_tot").get_time()),
                          fmt::format("{:>.1f}s", tid::get_unscoped("fLBIT")["run"].get_lap()));
    report += fmt::format("mem[rss {:<.1f}|peak {:<.1f}|vm {:<.1f}]MB ", debug::mem_rss_in_mb(), debug::mem_hwm_in_mb(), debug::mem_vm_in_mb());
    tools::log->info(report);
}