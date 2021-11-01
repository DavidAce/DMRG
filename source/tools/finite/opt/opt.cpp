#include "../opt_meta.h"
#include "../opt_mps.h"
#include "ceres_base.h"
#include "ceres_direct_functor.h"
#include "ceres_subspace_functor.h"
#include "opt-internal.h"
#include "report.h"
#include <algorithms/AlgorithmStatus.h>
#include <config/debug.h>
#include <glog/logging.h>
#include <math/num.h>
#include <string>
#include <tensors/model/ModelFinite.h>
#include <tensors/TensorsFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>

namespace settings {
    constexpr static bool debug_functor = false;
}

tools::finite::opt::opt_mps tools::finite::opt::find_excited_state(const TensorsFinite &tensors, const AlgorithmStatus &status, OptMeta &meta) {
    auto    t_opt = tid::tic_scope("opt");
    opt_mps initial_mps("current state", tensors.get_multisite_mps(), tensors.active_sites,
                        //                        tools::finite::measure::energy(tensors) - energy_reduced, // Eigval
                        tools::finite::measure::energy_minus_energy_reduced(tensors), // Eigval
                        tools::finite::measure::energy_reduced(tensors),              // Energy reduced for full system
                        tools::finite::measure::energy_variance(tensors),
                        1.0, // Overlap
                        tensors.get_length());
    initial_mps.set_max_grad(tools::finite::measure::max_gradient(initial_mps.get_tensor(), tensors));
    t_opt.toc(); // The next call to find_excited state opens the "opt" scope again.
    return find_excited_state(tensors, initial_mps, status, meta);
}

tools::finite::opt::opt_mps tools::finite::opt::find_excited_state(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status,
                                                                   OptMeta &meta) {
    auto t_opt = tid::tic_scope("opt");
    tools::log->trace("Starting optimization: mode [{}] | space [{}] | type [{}] | position [{}] | sites {} | shape {} = {}", enum2sv(meta.optMode),
                      enum2sv(meta.optSpace), enum2sv(meta.optType), status.position, tensors.active_sites, tensors.active_problem_dims(),
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
    ceres_default_options.max_lbfgs_rank                             = 16; // Tested: around 8-32 seems to be a good compromise, anything larger incurs a large overhead. The overhead means 2x computation time at ~64
    ceres_default_options.use_approximate_eigenvalue_bfgs_scaling    = true;  // Tested: True makes a huge difference, takes longer steps at each iteration and generally converges faster/to better variance
    ceres_default_options.min_line_search_step_size                  = 1e-16;//std::numeric_limits<double>::epsilon();
    ceres_default_options.max_line_search_step_contraction           = 1e-3; // 1e-3
    ceres_default_options.min_line_search_step_contraction           = 0.6; // 0.6
    ceres_default_options.max_line_search_step_expansion             = 10; // 10
    ceres_default_options.max_num_line_search_step_size_iterations   = 200;
    ceres_default_options.max_num_line_search_direction_restarts     = 10; //5
    ceres_default_options.line_search_sufficient_function_decrease   = 1e-4; //1e-4; Tested, doesn't seem to matter between [1e-1 to 1e-4]. Default is fine: 1e-4
    ceres_default_options.line_search_sufficient_curvature_decrease  = 0.8;//0.9 // This one should be above 0.5. Below, it makes retries at every step and starts taking twice as long for no added benefit. Tested 0.9 to be sweetspot
    ceres_default_options.max_solver_time_in_seconds                 = 60*60;//60*2;
    ceres_default_options.function_tolerance                         = 1e-10; // Tested, 1e-6 seems to be a sweetspot
    ceres_default_options.gradient_tolerance                         = 1e-1;
    ceres_default_options.parameter_tolerance                        = 1e-16;
    ceres_default_options.minimizer_progress_to_stdout               = false; //tools::log->level() <= spdlog::level::trace;
    ceres_default_options.update_state_every_iteration               = false;
    ceres_default_options.logging_type                               = ceres::LoggingType::PER_MINIMIZER_ITERATION;

    if(status.algorithm_has_stuck_for > 0){
        // Eigenvalue bfgs scaling performs badly when the problem is ill-conditioned (sensitive to some parameters).
        // Empirically, we observe that the gradient is still large when lbfgs has finished,
        // and often the result is a tiny bit worse than what we started with.
        // When this happens, it's not worth trying to get LBFGS to converge: instead, try an eigensolver on the reduced operator H²
        ceres_default_options.max_lbfgs_rank                         = 32; // Tested: around 8-32 seems to be a good compromise,but larger is more precise sometimes. Overhead goes from 1.2x to 2x computation time at in 8 -> 64
    }

    /* clang-format off */
    // Apply overrides if there are any
    if(meta.max_grad_tolerance)       ceres_default_options.gradient_tolerance                      = meta.max_grad_tolerance.value();
    if(meta.ceres_max_num_iterations) ceres_default_options.max_num_iterations                      = meta.ceres_max_num_iterations.value();
    if(meta.ceres_max_lbfgs_rank)     ceres_default_options.max_lbfgs_rank                          = meta.ceres_max_lbfgs_rank.value();
    if(meta.ceres_function_tolerance) ceres_default_options.function_tolerance                      = meta.ceres_function_tolerance.value();
    if(meta.ceres_eigenvalue_scaling) ceres_default_options.use_approximate_eigenvalue_bfgs_scaling = meta.ceres_eigenvalue_scaling.value();
    /* clang-format on */

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

    opt_mps result;
    switch(meta.optSpace) {
        /* clang-format off */
        case OptSpace::SUBSPACE:   result = internal::ceres_subspace_optimization(tensors, initial_mps, status, meta); break;
        case OptSpace::DIRECT:     result = internal::ceres_direct_optimization(tensors, initial_mps, status, meta); break;
        case OptSpace::KRYLOV:     result = internal::krylov_optimization(tensors,initial_mps, status, meta); break;
            /* clang-format on */
    }

    tools::finite::opt::reports::print_krylov_report();
    tools::finite::opt::reports::print_bfgs_report();
    tools::finite::opt::reports::print_time_report();

    // Finish up
    result.set_optspace(meta.optSpace);
    result.set_optmode(meta.optMode);
    result.set_optexit(meta.optExit);
    result.validate_result();
    return result;
}

tools::finite::opt::opt_mps tools::finite::opt::find_ground_state(const TensorsFinite &tensors, const AlgorithmStatus &status, OptMeta &conf) {
    return internal::ground_state_optimization(tensors, status, conf);
}

double tools::finite::opt::internal::windowed_func_abs(double x, double window) {
    if(std::abs(x) >= window)
        return std::abs(x) - window;
    else
        return 0;
}
double tools::finite::opt::internal::windowed_grad_abs(double x, double window) {
    if(std::abs(x) >= window)
        return num::sign(x);
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
    if constexpr(settings::debug) {
        freq_log_iter = 1;
        freq_log_time = 0;
        init_log_time = 0;
    }

    if constexpr(settings::debug_functor) {
        log->set_level(spdlog::level::debug);
        freq_log_iter = 1;
        freq_log_time = 0;
        init_log_time = 0;
    }
}

template<typename FunctorType>
ceres::CallbackReturnType tools::finite::opt::internal::CustomLogCallback<FunctorType>::operator()(const ceres::IterationSummary &summary) {
    functor.set_delta_f(std::abs(summary.cost_change / summary.cost));
    if(not log) return ceres::SOLVER_CONTINUE;
    if(log->level() >= spdlog::level::info) return ceres::SOLVER_CONTINUE;
    bool log_iter = summary.iteration - last_log_iter >= freq_log_iter and                  // log every freq_log_iter
                    summary.iteration >= init_log_iter;                                     // when at least init_log_iter have passed
    bool log_time = summary.cumulative_time_in_seconds - last_log_time >= freq_log_time and // log every last_log_time
                    summary.cumulative_time_in_seconds >= init_log_time;                    // when at least init_log_time have passed
    if(not log_iter or not log_time) return ceres::SOLVER_CONTINUE;
    last_log_time = summary.cumulative_time_in_seconds;
    last_log_iter = summary.iteration;
    /* clang-format off */
    log->debug(FMT_STRING("LBFGS: "
               "it {:>5} f {:>8.5f} |Δf| {:>8.2e} "
               "|∇f| {:>8.2e} |∇f|ᵐᵃˣ {:>8.2e} "
               "|ΔΨ| {:8.2e} |Ψ|-1 {:8.2e} "
               "ops {:>4}/{:<4} t {:>8.2e} s op {:>8.2e} s Gop/s {:>5.1f} "
               "ls:[ss {:8.2e} fe {:>4} ge {:>4} it {:>4}] "
//               "| energy {:<18.15f} lg var {:<6.6f}"
               ),
               summary.iteration,
               summary.cost,
               summary.cost_change,
               summary.gradient_norm,
               summary.gradient_max_norm,
               summary.step_norm, // By lbfgs
               functor.get_norm_offset(),
               functor.get_count() - last_count,
               functor.get_count(),
               summary.cumulative_time_in_seconds,
               summary.iteration_time_in_seconds,
               static_cast<double>(functor.get_ops()) / functor.t_H2n->get_last_interval()/1e9,
               std::abs(summary.step_size), // By line search
               summary.line_search_function_evaluations,
               summary.line_search_gradient_evaluations,
               summary.line_search_iterations
//             ,functor.get_energy_per_site(),
//               std::log10(functor.get_variance()
               );
    last_count = functor.get_count();
    /* clang-format on */
    return ceres::SOLVER_CONTINUE;
}

namespace tools::finite::opt::internal {
    template class CustomLogCallback<ceres_subspace_functor<double>>;
    template class CustomLogCallback<ceres_direct_functor<double, LagrangeNorm::ON>>;
    template class CustomLogCallback<ceres_direct_functor<double, LagrangeNorm::OFF>>;
    template class CustomLogCallback<ceres_subspace_functor<std::complex<double>>>;
    template class CustomLogCallback<ceres_direct_functor<std::complex<double>, LagrangeNorm::ON>>;
    template class CustomLogCallback<ceres_direct_functor<std::complex<double>, LagrangeNorm::OFF>>;
}
