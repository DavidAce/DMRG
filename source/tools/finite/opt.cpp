//
// Created by david on 2019-03-18.
//
#include <algorithms/class_algorithm_status.h>
#include <glog/logging.h>
#include <string>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt-internal/ceres_base.h>
#include <tools/finite/opt-internal/ceres_direct_functor.h>
#include <tools/finite/opt-internal/ceres_subspace_functor.h>
#include <tools/finite/opt-internal/opt-internal.h>
#include <tools/finite/opt-internal/report.h>
#include <tools/finite/opt_mps.h>

tools::finite::opt::opt_mps tools::finite::opt::find_excited_state(const class_tensors_finite &tensors, const class_algorithm_status &status, OptMode optMode,
                                                                   OptSpace optSpace, OptType optType) {
    double  energy_reduced = tools::finite::measure::energy_reduced(tensors);
    opt_mps initial_mps("current state", tensors.get_multisite_mps(), tensors.active_sites,
                        tools::finite::measure::energy(tensors) - energy_reduced, // Eigval
                        energy_reduced,                                           // Energy reduced for full system
                        tools::finite::measure::energy_variance(tensors),
                        1.0, // Overlap
                        tensors.get_length());

    return find_excited_state(tensors, initial_mps, status, optMode, optSpace, optType);
}

tools::finite::opt::opt_mps tools::finite::opt::find_excited_state(const class_tensors_finite &tensors, const opt_mps &initial_mps,
                                                                   const class_algorithm_status &status, OptMode optMode, OptSpace optSpace, OptType optType) {
    auto t_opt = tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt"]->tic_token();
    tools::log->trace("Starting optimization: mode [{}] | space [{}] | type [{}] | position [{}] | sites {} | shape {} = {}", enum2str(optMode),
                      enum2str(optSpace), enum2str(optType), status.position, tensors.active_sites, tensors.active_problem_dims(),
                      tensors.active_problem_size());

    using namespace opt::internal;
    static bool googleLogginghasInitialized = false;
    if(not googleLogginghasInitialized) {
        googleLogginghasInitialized = true;
#if defined(CERCES_INTERNAL_MINIGLOG_GLOG_LOGGING_H_)
        google::InitGoogleLogging(const_cast<char *>(tools::log->name().c_str()));
#else
        google::InitGoogleLogging(tools::log->name().c_str());
        google::SetStderrLogging(3);
//#pragma message "Using verbose glog"
//        google::SetStderrLogging(0);
#endif
    }

    /* clang-format off */

    /* Alglib has a good read comparing NCG and LBFGS https://www.alglib.net/optimization/lbfgsandcg.php
     * Some noteworthy points:
     *      "L-BFGS algorithm is quite easy to tune.
     *       The only parameter to tune is number of correction pairs M - number of function/gradient pairs
     *       used to build quadratic model (i.e. rank). This parameter is passed to the function which is used
     *       to create optimizer object. Increase of M decreases number of function evaluations and iterations
     *       needed to converge, but increases computational overhead associated with iteration. On well-conditioned
     *       problems M can be as small as 3-10. If function changes rapidly in some directions and slowly in
     *       other ones, then you can try increasing M in order to increase convergence."
     *
     *  Testing here shows that lbfgs is much better than NCG.
     *  Testing also shows that increasing lbfgs rank can give better results at some cost in some cases,
     *  in particular when a problem is ill-conditioned.
     *
     */

    /*
     * When a problem is ill-conditioned, for instance when the difference in eigenvalues is small,
     * we need to increase the lbfgs rank.
     * Read more here:
     *   kth.diva-portal.org/smash/get/diva2:1438308/FULLTEXT01.pdf (page 22)
     *   https://www.deeplearningbook.org/contents/optimization.html (page 307)
     *
     */


    ceres_default_options.line_search_type                           = ceres::LineSearchType::WOLFE;
    ceres_default_options.line_search_interpolation_type             = ceres::LineSearchInterpolationType::CUBIC;
    ceres_default_options.line_search_direction_type                 = ceres::LineSearchDirectionType::LBFGS;
    ceres_default_options.max_num_iterations                         = 20000;
    ceres_default_options.max_lbfgs_rank                             = 128; // Tested: around 8-32 seems to be a good compromise, anything larger incurs a large overhead. The overhead means 2x computation time at ~64
    ceres_default_options.use_approximate_eigenvalue_bfgs_scaling    = true;  // Tested: True makes a huge difference, takes longer steps at each iteration and generally converges faster/to better variance
    ceres_default_options.min_line_search_step_size                  = 1e-20;//std::numeric_limits<double>::epsilon();
    ceres_default_options.max_line_search_step_contraction           = 1e-3; // 1e-3
    ceres_default_options.min_line_search_step_contraction           = 0.6; // 0.6
    ceres_default_options.max_line_search_step_expansion             = 10; // 10
    ceres_default_options.max_num_line_search_step_size_iterations   = 200;
    ceres_default_options.max_num_line_search_direction_restarts     = 80; //5
    ceres_default_options.line_search_sufficient_function_decrease   = 1e-4; //Tested, doesn't seem to matter between [1e-1 to 1e-4]. Default is fine: 1e-4
    ceres_default_options.line_search_sufficient_curvature_decrease  = 0.9;//0.9 // This one should be above 0.5. Below, it makes retries at every step and starts taking twice as long for no added benefit. Tested 0.9 to be sweetspot
    ceres_default_options.max_solver_time_in_seconds                 = 60*60;//60*2;
    ceres_default_options.function_tolerance                         = 1e-14; // Tested, 1e-6 seems to be a sweetspot
    ceres_default_options.gradient_tolerance                         = 1e-14;
    ceres_default_options.parameter_tolerance                        = 1e-20;
    ceres_default_options.minimizer_progress_to_stdout               = false; //tools::log->level() <= spdlog::level::trace;
    ceres_default_options.logging_type                               = ceres::LoggingType::PER_MINIMIZER_ITERATION;

    if(status.algorithm_has_stuck_for > 0){
        ceres_default_options.max_lbfgs_rank                            = 128; // Tested: around 8-32 seems to be a good compromise,but larger is more precise sometimes. Overhead goes from 1.2x to 2x computation time at in 8 -> 64
    }



    //    Progress log definitions:
    //    f is the value of the objective function.
    //    d is the change in the value of the objective function if the step computed in this iteration is accepted.
    //    g is the inf norm of the gradient (i.e. the largest element of |grad f| )
    //    h is the change in the parameter vector.
    //    s is the optimal step length computed by the line search.
    //    it is the time take by the current iteration.
    //    tt is the total time taken by the minimizer.

    /* clang-format on */
    if(initial_mps.get_sites() != tensors.active_sites)
        throw std::runtime_error(fmt::format("mismatch in active sites: initial_mps {} | active {}", initial_mps.get_sites(), tensors.active_sites));

    initial_mps.validate_candidate();
    opt_mps result;
    switch(optSpace) {
        /* clang-format off */
        case OptSpace::SUBSPACE_ONLY:       result = internal::ceres_subspace_optimization(tensors,initial_mps,status, optType, optMode,optSpace); break;
        case OptSpace::SUBSPACE_AND_DIRECT: result = internal::ceres_subspace_optimization(tensors,initial_mps,status, optType, optMode,optSpace); break;
        case OptSpace::DIRECT:              result = internal::ceres_direct_optimization(tensors,initial_mps,status, optType,optMode,optSpace); break;
        case OptSpace::KRYLOV_ENERGY:        result = internal::krylov_energy_optimization(tensors,initial_mps,status, optType,optMode,optSpace); break;
        case OptSpace::KRYLOV_VARIANCE:      result = internal::krylov_variance_optimization(tensors,initial_mps,status, optType,optMode,optSpace); break;
            /* clang-format on */
    }
    // Finish up and print reports
    reports::print_bfgs_report();
    reports::print_time_report();
    result.set_optspace(optSpace);
    result.set_optmode(optMode);
    result.validate_result();
    return result;
}

tools::finite::opt::opt_mps tools::finite::opt::find_ground_state(const class_tensors_finite &tensors, const class_algorithm_status & status,  StateRitz ritz) {
    return internal::ground_state_optimization(tensors,status, ritz);
}

double tools::finite::opt::internal::windowed_func_abs(double x, double window) {
    if(std::abs(x) >= window)
        return std::abs(x) - window;
    else
        return 0;
}
double tools::finite::opt::internal::windowed_grad_abs(double x, double window) {
    if(std::abs(x) >= window)
        return sgn(x);
    else
        return 0.0;
}

double tools::finite::opt::internal::windowed_func_pow(double x, double window) {
    if(std::abs(x) >= window)
        return x * x - window * window;
    else
        return 0.0;
}
double tools::finite::opt::internal::windowed_grad_pow(double x, double window) {
    if(std::abs(x) >= window)
        return 2.0 * x;
    else
        return 0.0;
}

std::pair<double, double> tools::finite::opt::internal::windowed_func_grad(double x, double window) {
    double func = 0;
    double grad = 0;
    if(std::abs(x) >= window) {
        if(x > 0) {
            func = (x - window) * (x - window);
            grad = 2 * (x - window);
        } else {
            func = (x + window) * (x + window);
            grad = 2 * (x + window);
        }
    }

    return std::make_pair(func, grad);
}

long tools::finite::opt::internal::get_ops(long d, long chiL, long chiR, long m) {
    if(chiR > chiL) return get_ops(d, chiR, chiL, m);
    if(d > m) {
        // d first
        long step1 = chiL * chiL * chiR * m * m * d;
        long step2 = chiL * chiR * d * m * m * m * (d * m + 1);
        long step3 = chiL * chiR * d * m * (chiR * m * m * m + m * m + 1);
        return step1 + step2 + step3;
    } else {
        // m first
        long step1 = chiL * chiL * chiR * m * d;
        long step2 = chiL * chiR * d * m * m * (d * m + 1);
        long step3 = chiL * chiR * chiR * d * m * (m + 1);
        return step1 + step2 + step3;
    }
}

template<typename FunctorType>
tools::finite::opt::internal::CustomLogCallback<FunctorType>::CustomLogCallback(const FunctorType &functor_) : functor(functor_) {
    if(not log) log = tools::Logger::setLogger("xDMRG");
    log->set_level(tools::log->level());
    if(log->level() == spdlog::level::debug) {
        freq_log_iter = 10;
        freq_log_time = 10;
        init_log_time = 1;
    }
}

template<typename FunctorType>
ceres::CallbackReturnType tools::finite::opt::internal::CustomLogCallback<FunctorType>::operator()(const ceres::IterationSummary &summary) {
    //    if(summary.iteration > 100 and summary.cost_change < 1e-5 and summary.gradient_max_norm < 10.0 and summary.step_norm < 1e-6) return
    //    ceres::SOLVER_ABORT;
    functor.set_delta_f(std::abs(summary.cost_change/summary.cost));
    if(not log) return ceres::SOLVER_CONTINUE;
    if(log->level() >= spdlog::level::debug) return ceres::SOLVER_CONTINUE;
    if(summary.iteration - last_log_iter < freq_log_iter and summary.iteration > init_log_iter) return ceres::SOLVER_CONTINUE;
    if(summary.cumulative_time_in_seconds - last_log_time < freq_log_time and summary.cumulative_time_in_seconds > init_log_time) return ceres::SOLVER_CONTINUE;
    last_log_time = summary.cumulative_time_in_seconds;
    last_log_iter = summary.iteration;
    /* clang-format off */
    log->trace("LBFGS: iter {:>5} f {:>8.5f} |Δf| {:>3.2e} |∇f|∞ {:>3.2e} "
               "|ΔΨ| {:3.2e} |Ψ|-1 {:3.2e} ls {:3.2e} evals {:>4}/{:<4} "
               "t_step {:<} t_iter {:<} t_tot {:<} GOp/s {:<4.2f} | energy {:<18.15f} log₁₀var {:<6.6f}",
               summary.iteration,
               summary.cost,
               summary.cost_change,
               summary.gradient_max_norm,
               summary.step_norm, // By lbfgs
               functor.get_norm_offset(),
               std::abs(summary.step_size), // By line search
               functor.get_count() - last_count,
               functor.get_count(),
               fmt::format("{:>6.0f} ms",summary.step_solver_time_in_seconds * 1000),
               fmt::format("{:>6.0f} ms",summary.iteration_time_in_seconds * 1000),
               fmt::format("{:>5.3f} s",summary.cumulative_time_in_seconds),
               static_cast<double>(functor.get_ops()) / functor.t_vH2->get_last_interval()/1e9,
               functor.get_energy_per_site(),
               std::log10(functor.get_variance())
    );
    last_count = functor.get_count();
    /* clang-format on */
    return ceres::SOLVER_CONTINUE;
}

template class tools::finite::opt::internal::CustomLogCallback<tools::finite::opt::internal::ceres_subspace_functor<double>>;
template class tools::finite::opt::internal::CustomLogCallback<tools::finite::opt::internal::ceres_direct_functor<double>>;
template class tools::finite::opt::internal::CustomLogCallback<tools::finite::opt::internal::ceres_subspace_functor<std::complex<double>>>;
template class tools::finite::opt::internal::CustomLogCallback<tools::finite::opt::internal::ceres_direct_functor<std::complex<double>>>;