#include "xdmrg.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "fdmrg.h"
#include "general/iter.h"
#include "io/fmt_custom.h"
#include "math/eig.h"
#include "math/linalg/matrix.h"
#include "math/num.h"
#include "math/rnd.h"
#include "math/svd.h"
#include "qm/time.h"
#include "tensors/edges/EdgesFinite.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/site/mpo/MpoSite.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include "tools/common/prof.h"
#include "tools/finite/h5.h"
#include "tools/finite/measure.h"
#include "tools/finite/mpo.h"
#include "tools/finite/mps.h"
#include "tools/finite/multisite.h"
#include "tools/finite/opt.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"
#include "tools/finite/print.h"
#include <h5pp/details/h5ppFile.h>

xdmrg::xdmrg(std::shared_ptr<h5pp::File> h5ppFile_) : AlgorithmFinite(std::move(h5ppFile_), settings::xdmrg::ritz, AlgorithmType::xDMRG) {
    tools::log->trace("Constructing class_xdmrg");
    tensors.state->set_name("state_emid");
}

void xdmrg::resume() {
    // Resume can imply many things
    // 1) Resume a simulation which terminated prematurely
    // 2) Resume a previously successful simulation. This may be desireable if the config
    //    wants something that is not present in the file.
    //      a) A certain number of states
    //      b) A state inside a particular energy window
    //      c) The ground or "roof" states
    // To guide the behavior, we check the setting ResumePolicy.

    auto states_that_may_resume =
        tools::common::h5::resume::find_states_that_may_resume(*h5file, settings::storage::resume_policy, status.algo_type, "state_emid");
    if(states_that_may_resume.empty()) throw except::state_error("no resumable states were found");
    for(const auto &[state_prefix, algo_stop] : states_that_may_resume) {
        tools::log->info("Resuming [{}] | previous stop reason: {} | resume policy: {} ", state_prefix, enum2sv(algo_stop),
                         enum2sv(settings::storage::resume_policy));
        try {
            tools::finite::h5::load::simulation(*h5file, state_prefix, tensors, status, status.algo_type);
        } catch(const except::load_error &le) {
            tools::log->error("{}", le.what());
            continue;
        }

        // Our first task is to decide on a state name for the newly loaded state
        // The simplest is to infer it from the state prefix itself
        auto name = tools::common::h5::resume::extract_state_name(state_prefix);
        tensors.state->set_name(name);

        // Set a possibly new energy target
        // status.energy_tgt = settings::xdmrg::energy_spectrum_shift;

        // Reload the bond and truncation error limits (could be different in the config compared to the status we just loaded)
        double long_max                   = static_cast<double>(std::numeric_limits<long>::max());
        double bond_max                   = std::min(long_max, std::pow(2.0, settings::model::model_size / 2));
        status.bond_max                   = std::min({status.bond_max, safe_cast<long>(bond_max), settings::get_bond_max(status.algo_type)});
        status.bond_min                   = std::max(status.bond_min, settings::get_bond_min(status.algo_type));
        status.bond_lim                   = std::min(status.bond_lim, status.bond_max);
        status.bond_limit_has_reached_max = status.bond_lim == status.bond_max;

        status.trnc_min                   = settings::precision::svd_truncation_min;
        status.trnc_max                   = settings::precision::svd_truncation_max;
        status.trnc_limit_has_reached_min = status.trnc_lim == status.trnc_min;

        // Apply shifts and compress the model
        tensors.move_center_point_to_inward_edge();
        set_parity_shift_mpo();
        set_parity_shift_mpo_squared();
        set_energy_shift_mpo();
        rebuild_tensors(); // Rebuilds and compresses mpos, then rebuilds the environments
        update_precision_limit();
        update_dmrg_blocksize();
        // Initialize a custom task list
        std::deque<xdmrg_task> task_list;

        if(status.algorithm_has_succeeded)
            task_list = {xdmrg_task::POST_PRINT_RESULT};
        else
            task_list = {xdmrg_task::INIT_CLEAR_CONVERGENCE, xdmrg_task::FIND_EXCITED_STATE,
                         xdmrg_task::POST_DEFAULT}; // Probably a savepoint. Simply "continue" the algorithm until convergence
        run_task_list(task_list);
    }
}

void xdmrg::run_default_task_list() {
    std::deque<xdmrg_task> default_task_list = {
        xdmrg_task::INIT_DEFAULT,
        xdmrg_task::FIND_EXCITED_STATE,
        xdmrg_task::POST_DEFAULT,
    };

    run_task_list(default_task_list);
}

void xdmrg::run_task_list(std::deque<xdmrg_task> &task_list) {
    while(not task_list.empty()) {
        auto task = task_list.front();
        switch(task) {
            case xdmrg_task::INIT_RANDOMIZE_MODEL: initialize_model(); break;
            case xdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE: initialize_state(ResetReason::INIT, StateInit::RANDOM_PRODUCT_STATE); break;
            case xdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE: initialize_state(ResetReason::INIT, StateInit::RANDOM_ENTANGLED_STATE); break;
            case xdmrg_task::INIT_RANDOMIZE_FROM_CURRENT_STATE: initialize_state(ResetReason::INIT, StateInit::RANDOMIZE_PREVIOUS_STATE); break;
            case xdmrg_task::INIT_BOND_LIMITS: init_bond_dimension_limits(); break;
            case xdmrg_task::INIT_TRNC_LIMITS: init_truncation_error_limits(); break;
            case xdmrg_task::INIT_ENERGY_TARGET: init_energy_target(); break;
            case xdmrg_task::INIT_WRITE_MODEL: write_to_file(StorageEvent::MODEL); break;
            case xdmrg_task::INIT_CLEAR_STATUS: status.clear(); break;
            case xdmrg_task::INIT_CLEAR_CONVERGENCE: clear_convergence_status(); break;
            case xdmrg_task::INIT_DEFAULT: run_preprocessing(); break;

            case xdmrg_task::FIND_ENERGY_RANGE: find_energy_range(); break;
            case xdmrg_task::FIND_EXCITED_STATE:
                tensors.state->set_name("state_emid");
                run_algorithm();
                break;
            case xdmrg_task::POST_WRITE_RESULT: write_to_file(StorageEvent::FINISHED, CopyPolicy::FORCE); break;
            case xdmrg_task::POST_PRINT_RESULT: print_status_full(); break;
            case xdmrg_task::POST_PRINT_TIMERS: tools::common::timer::print_timers(); break;
            case xdmrg_task::POST_RBDS_ANALYSIS: run_rbds_analysis(); break;
            case xdmrg_task::POST_RTES_ANALYSIS: run_rtes_analysis(); break;
            case xdmrg_task::POST_DEFAULT: run_postprocessing(); break;
            case xdmrg_task::TIMER_RESET: tid::reset("xDMRG"); break;
        }
        task_list.pop_front();
    }
    if(not task_list.empty()) {
        for(auto &task : task_list) tools::log->critical("Unfinished task: {}", enum2sv(task));
        throw except::runtime_error("Simulation ended with unfinished tasks");
    }
}

void xdmrg::init_energy_target(std::optional<double> energy_density_target) {
    switch(status.opt_ritz) {
        case OptRitz::NONE: throw std::logic_error("status.opt_ritz == OptRitz::NONE is invalid under xdmrg");
        case OptRitz::SR: {
            tools::log->warn("status.opt_ritz == OptRitz::SR should be handled with fdmrg instead of xdmrg");
            status.energy_tgt = 0.0;
            break;
            // throw std::logic_error("status.opt_ritz == OptRitz::SR should be handled with fdmrg instead of xdmrg");
        }
        case OptRitz::LR: {
            tools::log->warn("status.opt_ritz == OptRitz::LR should be handled with fdmrg instead of xdmrg");
            status.energy_tgt = 0.0;
            break;
            // throw std::logic_error("status.opt_ritz == OptRitz::LR should be handled with fdmrg instead of xdmrg");
        }
        case OptRitz::SM: {
            // When the Hamiltonian is traceless, the energy level nearest zero is closest to the infinite-temperature limit.
            // Therefore we expect the energy target to be == 0. However, in some cases we get a symmetric energy spectrum, with
            // every energy level having a counterpart with opposite sign (e.g. Ising-Majorana with g == 0).
            // In this case we can break the degeneracy by setting a tiny shift ~1e-10 to bias xDMRG towards one of the states closest to 0.
            status.energy_tgt = settings::xdmrg::energy_spectrum_shift;
            break;
        }
        case OptRitz::IS: {
            status.energy_tgt = tools::finite::measure::energy(tensors); // Should take the energy from the initial state
            break;
        }
        case OptRitz::TE: {
            if(not energy_density_target) energy_density_target = settings::xdmrg::energy_density_target;
            if(energy_density_target.value() < 0.0 or energy_density_target.value() > 1.0)
                throw except::runtime_error(fmt::format(
                    "xdmrg::init_energy_target: with OptRitz::TE: invalid energy_density_target: Expected value in range [0.0 - 1.0], got: [{:.8f}]",
                    energy_density_target));
            // Set energy boundaries. This function is supposed to run after find_energy_range!
            if(status.energy_max == status.energy_min)
                throw except::runtime_error("xdmrg::init_energy_target: with OptRitz::TE Failed because energy_max == {} and energy_min == {}\n"
                                            "Try running find_energy_range() first",
                                            status.energy_max, status.energy_min);

            status.energy_dens_target = energy_density_target.value();
            status.energy_tgt         = status.energy_min + status.energy_dens_target * (status.energy_max - status.energy_min);
            tools::log->info("Energy minimum     = {:.8f}", status.energy_min);
            tools::log->info("Energy maximum     = {:.8f}", status.energy_max);
            tools::log->info("Energy target      = {:.8f}", status.energy_tgt);
            break;
        }
    }
}

void xdmrg::run_preprocessing() {
    tools::log->info("Running {} preprocessing", status.algo_type_sv());
    auto t_pre = tid::tic_scope("pre");
    status.clear();
    init_bond_dimension_limits();
    init_truncation_error_limits();
    initialize_model(); // First use of random!

    initialize_state(ResetReason::INIT, settings::strategy::initial_state); // Second use of random!
    find_energy_range();
    init_energy_target();
    set_parity_shift_mpo();
    set_parity_shift_mpo_squared();
    set_energy_shift_mpo();
    rebuild_tensors(); // Rebuilds and compresses mpos, then rebuilds the environments
    update_precision_limit();
    write_to_file(StorageEvent::MODEL);

    // auto imodel = tools::finite::mpo::get_inverted_mpos(tensors.model->get_all_mpo_tensors(MposWithEdges::ON));

    if(tensors.get_length<long>() <= 4) {
        // Print the spectrum if small
        // tensors.clear_cache();
        auto svd_solver = svd::solver();
        auto L          = tensors.get_length<long>();
        auto sites      = num::range<size_t>(0, L);
        auto ham1       = tensors.model->get_multisite_ham(sites);
        auto ham2       = tensors.model->get_multisite_ham_squared(sites);
        auto norm_est   = tensors.model->get_energy_upper_bound();
        // auto        ham1i      = svd_solver.pseudo_inverse(ham1_);
        eig::solver solver1, solver2;
        solver1.eig<eig::Form::SYMM>(ham1.data(), ham1.dimension(0), eig::Vecs::OFF);
        solver2.eig<eig::Form::SYMM>(ham2.data(), ham2.dimension(0), eig::Vecs::OFF);
        // solver1i.eig<eig::Form::SYMM>(ham1i.data(), ham1i.dimension(0));

        auto            evals1 = eig::view::get_eigvals<real>(solver1.result);
        auto            evals2 = eig::view::get_eigvals<real>(solver2.result);
        Eigen::VectorXd diffs1 = Eigen::VectorXd::Zero(evals1.size());
        Eigen::VectorXd diffs2 = Eigen::VectorXd::Zero(evals1.size());
        auto            N1     = evals1.size() - 1;
        auto            N2     = evals2.size() - 1;
        diffs1.topRows(N1)     = (evals1.bottomRows(N1) - evals1.topRows(N1));
        diffs2.topRows(N2)     = (evals2.bottomRows(N2) - evals2.topRows(N2));

        // auto evals1i = eig::view::get_eigvals<real>(solver1i.result);
        fmt::print("{:^8} {:<20}\n", " ", "H¹");
        for(long idx = 0; idx < evals1.size(); ++idx) {
            // if(std::abs(evals1[idx]) > 1.1) continue;
            fmt::print("idx {:2}: {:20.16f} {:>10.3e}\n", idx, evals1[idx], diffs1[idx]);
        }
        fmt::print("{:^8} {:<20} {:<20} {:<20}\n", " ", "H²", "diff", "sqrt(H²)");
        for(long idx = 0; idx < evals2.size(); ++idx) {
            // if(std::abs(evals2[idx]) > 1.1) continue;
            fmt::print("idx {:2}: {:20.16f} {:>10.3e} {:20.16f}\n", idx, evals2[idx], diffs2[idx], std::sqrt(evals2[idx]));
        }
        auto h5file = h5pp::File("../../output/spectrum.h5", h5pp::FileAccess::RENAME);
        h5file.writeDataset(evals1, "H_evals");
        h5file.writeDataset(diffs1, "H_diffs");
        h5file.writeDataset(evals2, "H2_evals");
        h5file.writeDataset(diffs2, "H2_diffs");
        // fmt::print("{:^8} {:<20} {:<20}\n", " ", "H¹", "H⁻¹");
        // for(long idx = 0; idx < std::min(evals1_.size(), evals1i.size()); ++idx) {
        //     fmt::print("idx {:2}: {:20.16f} {:20.16f}\n", idx, evals1_[idx], evals1i[idx]);
        // }
        fmt::print("\n");
        fmt::print("Hamiltonian norm estimate: {:.16f}\n", norm_est);
        // Try the iterative scheme
        // auto impos = tools::finite::mpo::get_inverted_mpos(tensors.model->get_compressed_mpos_squared(MposWithEdges::ON));
        // auto impos = tools::finite::mpo::get_inverted_mpos(tensors.model->get_all_mpo_tensors(MposWithEdges::ON));

        // auto imodel = *tensors.model;
        // for(auto &&[pos, impo] : iter::enumerate(imodel.MPO)) impo->set_mpo(impos[pos]);
        // auto ham3i = imodel.get_multisite_ham(sites);
        // // auto ham3i = svd::solver::pseudo_inverse(ham3i);
        // eig::solver solver3i;
        // solver3i.eig<eig::Form::SYMM>(ham3i.data(), ham3i.dimension(0));
        // auto evals3i = eig::view::get_eigvals<real>(solver3i.result);

        // fmt::print("{:^8} {:<20} {:<20} {:<20} {:<20}\n", " ", "H¹", "H⁻¹", "H⁻¹(iter)", "diff");
        // for(long idx = 0; idx < std::min(evals1i.size(), evals3i.size()); ++idx) {
        //     fmt::print("idx {:2}: {:20.16f} {:20.16f} {:20.16f} {:.4e}\n", idx, evals1_[idx], evals1i[idx], evals3i[idx],
        //                std::abs(evals1i[idx] - evals3i[idx]));
        // }
        // fmt::print("\n");

        // tools::log->info("Iterative inverse  H⁻¹");
        // for(long idx = 0; idx < evals3i.size(); ++idx) { fmt::print("idx {:2}: {:20.16f}\n", idx, evals3i[idx]); }
        // exit(0);
    }
    tools::log->info("Finished {} preprocessing", status.algo_type_sv());
}

void xdmrg::run_algorithm() {
    if(tensors.state->get_name().empty()) tensors.state->set_name("state_emid");
    tools::log->info("Starting {} simulation of model [{}] for state [{}] with ritz [{}]", status.algo_type_sv(), enum2sv(settings::model::model_type),
                     tensors.state->get_name(), enum2sv(status.opt_ritz));
    auto t_run       = tid::tic_scope("run");
    status.algo_stop = AlgorithmStop::NONE;

    while(true) {
        tools::log->trace("Starting step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);
        update_state();
        print_status();
        check_convergence();
        write_to_file();

        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);

        // It's important not to perform the last move, so we break now: that last state would not get optimized
        if(status.algo_stop != AlgorithmStop::NONE) break;

        // Prepare for next step

        // Updating bond dimension must go first since it decides based on truncation error, but a projection+normalize resets truncation.
        update_bond_dimension_limit();   // Updates the bond dimension if the state precision is being limited by bond dimension
        update_truncation_error_limit(); // Updates the truncation error limit if the state is being truncated
        update_eigs_tolerance();         // Updates the tolerance on the iterative eigensolver
        update_dmrg_blocksize();         // Updates the number sites used in dmrg steps using the information typical scale
        try_retargeting();
        try_projection(); // Tries to project the state to the nearest global spin parity sector along settings::strategy::target_axis
        // try_moving_sites(); // Tries to overcome an entanglement barrier by moving sites around the lattice, to optimize non-nearest neighbors

        // expand_environment(EnvExpandMode::ENE, EnvExpandSide::BACKWARD);
        move_center_point(); // Moves the center point AC to the next site and increments status.iter and status.step
        status.wall_time = tid::get_unscoped("t_tot").get_time();
        status.algo_time = t_run->get_time();
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", status.algo_type_sv(), tensors.state->get_name(), status.algo_stop_sv());
    status.algorithm_has_finished = true;
    //    tools::finite::measure::parity_components(*tensors.state, qm::spin::half::sz);
}

void xdmrg::update_state() {
    using namespace tools::finite;
    using namespace tools::finite::opt;
    auto t_step   = tid::tic_scope("step");
    auto opt_meta = get_opt_meta();
    tools::log->debug("Starting {} iter {} | step {} | pos {} | dir {} | ritz {} | type {}", status.algo_type_sv(), status.iter, status.step, status.position,
                      status.direction, enum2sv(settings::get_ritz(status.algo_type)), enum2sv(opt_meta.optType));
    // Try activating the sites asked for;
    tensors.activate_sites(opt_meta.chosen_sites);
    if(tensors.active_sites.empty()) {
        tools::log->warn("Failed to activate sites");
        return;
    }
    tensors.rebuild_edges();

    tools::log->debug("Updating state: {}", opt_meta.string()); // Announce the current configuration for optimization
    auto bond_dims_old = tensors.state->get_mps_dims_active();
    // Expand the environment to grow the bond dimension in 1-site dmrg
    if(opt_meta.expand_mode != EnvExpandMode::NONE) {
        expand_environment(opt_meta.expand_mode, opt_meta.expand_side);
        opt_meta.problem_dims = tools::finite::multisite::get_dimensions(*tensors.state);
        opt_meta.problem_size = tools::finite::multisite::get_problem_size(*tensors.state);
        opt_meta.optSolver    = opt_meta.problem_size <= settings::precision::eig_max_size ? OptSolver::EIG : OptSolver::EIGS;
    }
    auto variance_after_exp = tools::finite::measure::energy_variance(tensors);
    auto bond_dims_exp      = tensors.state->get_mps_dims_active();

    // Run the optimization
    auto initial_state = opt::get_opt_initial_mps(tensors, opt_meta);
    auto opt_state     = opt::find_excited_state(tensors, initial_state, status, opt_meta);

    // Determine the quality of the optimized state.
    opt_state.set_relchange(opt_state.get_variance() / var_latest);
    opt_state.set_bond_limit(opt_meta.svd_cfg->rank_max.value());
    opt_state.set_trnc_limit(opt_meta.svd_cfg->truncation_limit.value());
    /* clang-format off */
    opt_meta.optExit = OptExit::SUCCESS;
    if(opt_state.get_grad_max()       > 1.000                         ) opt_meta.optExit |= OptExit::FAIL_GRADIENT;
    if(opt_state.get_rnorm()          > settings::precision::eigs_tol_max) opt_meta.optExit |= OptExit::FAIL_RESIDUAL;
    if(opt_state.get_eigs_nev()       == 0 and
       opt_meta.optSolver              == OptSolver::EIGS             ) opt_meta.optExit |= OptExit::FAIL_RESIDUAL; // No convergence
    if(opt_state.get_overlap()        < 0.010                         ) opt_meta.optExit |= OptExit::FAIL_OVERLAP;
    if(opt_state.get_relchange()      > 1.001                         ) opt_meta.optExit |= OptExit::FAIL_WORSENED;
    else if(opt_state.get_relchange() > 0.999                         ) opt_meta.optExit |= OptExit::FAIL_NOCHANGE;

    opt_state.set_optexit(opt_meta.optExit);
    /* clang-format on */

    // Estimate the convergence rate (requires nev > 1)
    double convr = std::numeric_limits<double>::quiet_NaN();
    if(opt_state.get_optcost() == OptCost::VARIANCE and settings::xdmrg::try_shifting_when_degen > 0) {
        // auto ev0   = opt_state.get_eigs_eigval();
        // auto ev1   = opt_state.next_evs.front().eigs_eigval;
        // auto evn   = std::pow(tensors.model->get_energy_upper_bound(), 2.0);
        // auto kappa = (evn - ev0) / (ev1 - ev0);
        // auto gamma = (std::sqrt(kappa) - 1) / (std::sqrt(kappa) + 1);

        // Compare the progress in the objective function
        auto old_eigvalue = initial_state.get_variance() + std::pow(initial_state.get_energy(), 2.0); // Plain H²
        auto new_eigvalue = opt_state.get_eigs_eigval();                                              // Should also be plain H²
        auto num_matvecs  = static_cast<double>(opt_state.get_mv());
        auto time_eigs    = opt_state.get_time();
        if(new_eigvalue < old_eigvalue and num_matvecs > 0 and time_eigs > 0.0) {
            convr = (1.0 - new_eigvalue / old_eigvalue) / num_matvecs / time_eigs;
            tools::log->info("H² {:.3e} -> {:.3e} in {:<10} matvecs: convr {:.3e}", old_eigvalue, new_eigvalue, num_matvecs, convr);
        }

        // convr      = std::log(0.5) / std::log(gamma);
    }
    for(auto site : opt_state.get_sites()) dmrg_degeneracy_score.at(site) = convr;

    tools::log->trace("Optimization [{}|{}]: {}. Variance change {:8.2e} --> {:8.2e} ({:.3f} %)", enum2sv(opt_meta.optCost), enum2sv(opt_meta.optSolver),
                      flag2str(opt_meta.optExit), var_latest, opt_state.get_variance(), opt_state.get_relchange() * 100);
    if(opt_state.get_relchange() > 1000) tools::log->warn("Variance increase by x {:.2e}", opt_state.get_relchange());

    if(tools::log->level() <= spdlog::level::debug) {
        tools::log->debug("Optimization result: {:<24} | E {:<20.16f}| σ²H {:<8.2e} | rnorm {:8.2e} | overlap {:.16f} | "
                          "sites {} | conv {:.2f} | {:20} | {} | time {:.2e} s",
                          opt_state.get_name(), opt_state.get_energy(), opt_state.get_variance(), opt_state.get_rnorm(), opt_state.get_overlap(),
                          opt_state.get_sites(), convr, fmt::format("[{}][{}]", enum2sv(opt_state.get_optcost()), enum2sv(opt_state.get_optsolver())),
                          flag2str(opt_state.get_optexit()), opt_state.get_time());
    }

    tensors.state->tag_active_sites_normalized(false);

    // Do the truncation with SVD
    // TODO: We may need to detect here whether the truncation error limit needs lowering due to a variance increase in the svd merge
    auto logPolicy = LogPolicy::SILENT;
    if constexpr(settings::debug) logPolicy = LogPolicy::VERBOSE;
    tensors.merge_multisite_mps(opt_state.get_tensor(), opt_meta.svd_cfg, logPolicy);
    tensors.rebuild_edges(); // This will only do work if edges were modified, which is the case in 1-site dmrg.
    if constexpr(settings::debug) {
        if(tools::log->level() <= spdlog::level::trace) tools::log->trace("Truncation errors: {::8.3e}", tensors.state->get_truncation_errors_active());
    }

    if constexpr(settings::debug) {
        auto variance_after_svd = tools::finite::measure::energy_variance(tensors);
        tools::log->debug("Before update: variance {:8.2e} | mps dims {}", var_latest, bond_dims_old);
        tools::log->debug("After  expns.: variance {:8.2e} | mps dims {}", variance_after_exp, bond_dims_exp);
        tools::log->debug("After  merge : variance {:8.2e} | mps dims {}", variance_after_svd, tensors.state->get_mps_dims_active());
        tools::log->debug("Variance change from  SVD: {:.16f}%", 100 * variance_after_svd / opt_state.get_variance());
    }

    // Update current energy density ε
    if(status.opt_ritz == OptRitz::TE)
        status.energy_dens = (tools::finite::measure::energy(tensors) - status.energy_min) / (status.energy_max - status.energy_min);

    tools::log->trace("Updating variance record holder");
    auto ene                      = tools::finite::measure::energy(tensors);
    auto var                      = tools::finite::measure::energy_variance(tensors);
    status.energy_variance_lowest = std::min(var, status.energy_variance_lowest);
    var_delta                     = var - var_latest;
    ene_delta                     = ene - ene_latest;
    var_change                    = var / var_latest;
    var_latest                    = var;
    ene_latest                    = ene;

    last_optsolver = opt_state.get_optsolver();
    last_optalgo   = opt_state.get_optalgo();
    last_optcost   = opt_state.get_optcost();
    if constexpr(settings::debug) tensors.assert_validity();
}

void xdmrg::find_energy_range() {
    // We only need to find an energy range if we are targeting a particular energy density window or target
    if(status.opt_ritz != OptRitz::TE) return; // We only need the extremal for OptRitz::TED

    tools::log->trace("Finding energy range");
    auto t_init = tid::tic_scope("init");
    // Here we define a set of tasks for fdmrg in order to produce the lowest and highest energy eigenstates,
    // We don't want it to randomize its own model, so we implant our current model before running the tasks.

    std::deque<fdmrg_task> gs_tasks = {fdmrg_task::INIT_CLEAR_STATUS, fdmrg_task::INIT_BOND_LIMITS, fdmrg_task::INIT_TRNC_LIMITS,
                                       fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE, fdmrg_task::FIND_GROUND_STATE};

    std::deque<fdmrg_task> hs_tasks = {fdmrg_task::INIT_CLEAR_STATUS, fdmrg_task::INIT_BOND_LIMITS, fdmrg_task::INIT_TRNC_LIMITS,
                                       fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE, fdmrg_task::FIND_HIGHEST_STATE};
    // Find the lowest energy state
    {
        auto  t_gs = tid::tic_scope("fDMRG");
        fdmrg fdmrg_gs{};
        *fdmrg_gs.tensors.model = *tensors.model; // Copy the model
        fdmrg_gs.tensors.state->set_name("state_emin");
        tools::log = tools::Logger::setLogger(fmt::format("{}-gs", status.algo_type_sv()), settings::console::loglevel, settings::console::timestamp);
        fdmrg_gs.run_task_list(gs_tasks);
        status.energy_min = tools::finite::measure::energy(fdmrg_gs.tensors);
        fdmrg_gs.h5file   = h5file;
        write_to_file(*fdmrg_gs.tensors.state, *fdmrg_gs.tensors.model, *fdmrg_gs.tensors.edges, StorageEvent::EMIN, CopyPolicy::OFF);
    }

    // Find the highest energy state
    {
        auto  t_hs = tid::tic_scope("fDMRG");
        fdmrg fdmrg_hs{};
        *fdmrg_hs.tensors.model = *tensors.model; // Copy the model
        fdmrg_hs.tensors.state->set_name("state_emax");
        tools::log = tools::Logger::setLogger(fmt::format("{}-hs", status.algo_type_sv()), settings::console::loglevel, settings::console::timestamp);
        fdmrg_hs.run_task_list(hs_tasks);
        status.energy_max = tools::finite::measure::energy(fdmrg_hs.tensors);
        fdmrg_hs.h5file   = h5file;
        write_to_file(*fdmrg_hs.tensors.state, *fdmrg_hs.tensors.model, *fdmrg_hs.tensors.edges, StorageEvent::EMAX, CopyPolicy::OFF);
    }

    // Reset our logger
    tools::log = tools::Logger::getLogger(fmt::format("{}", status.algo_type_sv()));
}

void xdmrg::set_energy_shift_mpo() {
    // In xdmrg we find an excited energy eigenstate by optimizing the energy variance of some state close to a target energy.
    // We can target a particular energy by setting an energy shift (equal to the target energy), which then becomes the energy minimum
    // once we fold the spectrum by squaring the Hamiltonian (i.e. we optimize (H-E_tgt)²):
    //      Var H = <(H-E_tgt)²> - <H-E_tgt>²     = <H²> - 2<H>E_tgt + E_tgt² - (<H> - E_tgt)²
    //                                            =  H²  - 2*E*E_tgt + E_tgt² - E² + 2*E*E_tgt - E_tgt²
    //                                            =  H²  - E²
    // The first term <(H-E_tgt)²> is computed using a double-layer of mpos with energy shifted by E_tgt.
    // If we didn't shift mpo's, the last line, H²-E² is subtraction of two large numbers --> catastrophic cancellation --> loss of precision.
    // However, by minimizing the variance of shifted mpos:
    //              Var H = <(H-E_tgt)²> - <H-E_tgt>² = <(H-E_tgt)²> - (E-E_tgt)²
    // we get the subtraction of two very small terms since E-E_shf should be small.

    if(not tensors.position_is_inward_edge()) return;
    tensors.set_energy_shift_mpo(status.energy_tgt);
    tensors.rebuild_mpo();          // The shift clears our squared mpo's. So we have to rebuild them.
    tensors.rebuild_mpo_squared();  // The shift clears our squared mpo's. So we have to rebuild them.
    tensors.compress_mpo_squared(); // Compress the mpo's if compression is enabled
    tensors.rebuild_edges();        // The shift modified all our mpo's. So we have to rebuild all the edges.
    if constexpr(settings::debug) tensors.assert_validity();
    if(std::abs(tensors.model->get_energy_shift_mpo() - status.energy_tgt) > 1e-15)
        throw except::runtime_error("Energy shift mismatch: {:.16f} != {:.16f}", tensors.model->get_energy_shift_mpo(), status.energy_tgt);
}

void xdmrg::try_retargeting() {
    if(not tensors.position_is_inward_edge()) return;
    if(status.algorithm_has_stuck_for == 0) return;
    if(settings::xdmrg::try_shifting_when_degen <= 0) return;
    auto convmax = num::nanmax(dmrg_degeneracy_score);
    if(std::isnan(convmax)) return;
    tools::log->info("convergence difficulty score: {:.1f} (limit {:.1f}) | {::.1f}", convmax, settings::xdmrg::try_shifting_when_degen, dmrg_degeneracy_score);
    if(convmax > settings::xdmrg::try_shifting_when_degen) {
        auto mu                                = 0.0;
        auto sig                               = tensors.model->get_energy_upper_bound() / 100.0;
        settings::xdmrg::energy_spectrum_shift = rnd::normal(mu, sig);
        status.energy_tgt                      = settings::xdmrg::energy_spectrum_shift;
        dmrg_degeneracy_score                  = std::vector(tensors.get_length<size_t>(), std::numeric_limits<double>::quiet_NaN());
        tools::finite::mps::normalize_state(*tensors.state, svd::config(settings::xdmrg::bond_min, status.trnc_max), NormPolicy::ALWAYS);
        clear_convergence_status();
        set_energy_shift_mpo();
    }
}

void xdmrg::update_time_step() {
    tools::log->trace("Updating time step");
    status.delta_t = std::complex<double>(1e-6, 0);
}


