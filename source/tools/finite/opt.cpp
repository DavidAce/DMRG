//
// Created by david on 2019-03-18.
//
#include <algorithms/class_algorithm_status.h>
#include <general/nmspc_tensor_extra.h>
#include <glog/logging.h>
#include <string>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt.h>
#include <tools/finite/opt_state.h>

tools::finite::opt::opt_state tools::finite::opt::find_excited_state(const class_tensors_finite &  tensors,
                                                                                       const class_algorithm_status &status, OptMode optMode, OptSpace optSpace, OptType optType) {
    std::vector<size_t> sites = tensors.active_sites;
    opt_state           initial_tensor("current state", tensors.state->get_multisite_mps(), sites,
                                       tools::finite::measure::energy(tensors) - tensors.model->get_energy_reduced(), // Eigval
                                       tensors.model->get_energy_reduced(),                                           // Energy reduced for full system
                                       tools::finite::measure::energy_variance(tensors),
                                       1.0, // Overlap
                                       tensors.get_length());

    return find_excited_state(tensors, initial_tensor, status, optMode, optSpace, optType);
}




tools::finite::opt::opt_state tools::finite::opt::find_excited_state(const class_tensors_finite &tensors, const opt_state & initial_tensor, const class_algorithm_status &status,
                                                                      OptMode optMode, OptSpace optSpace, OptType optType) {
    tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt"]->tic();
    tools::log->debug("Starting optimization: mode [{}] | space [{}] | type [{}] | position [{}] | shape {} = {}", enum2str(optMode), enum2str(optSpace), enum2str(optType),
                      tensors.get_position(), tensors.state->active_dimensions(), tensors.state->active_problem_size());

    using namespace opt::internal;
    static bool googleLogginghasInitialized = false;
    if(not googleLogginghasInitialized) {
        google::InitGoogleLogging(tools::log->name().c_str());
        googleLogginghasInitialized = true;
        google::SetStderrLogging(3);
    }

    /* clang-format off */
    ceres_default_options.line_search_type                           = ceres::LineSearchType::WOLFE;
    ceres_default_options.line_search_interpolation_type             = ceres::LineSearchInterpolationType::CUBIC;
    ceres_default_options.line_search_direction_type                 = ceres::LineSearchDirectionType::LBFGS;
    ceres_default_options.nonlinear_conjugate_gradient_type          = ceres::NonlinearConjugateGradientType::FLETCHER_REEVES;
    ceres_default_options.max_num_iterations                         = 4000;
    ceres_default_options.max_lbfgs_rank                             = 16; // Tested: around 8-32 seems to be a good compromise, anything larger incurs a large overhead. The overhead means 2x computation time at ~256
    ceres_default_options.use_approximate_eigenvalue_bfgs_scaling    = true;  // Tested: True makes a huge difference, takes longer steps at each iteration and generally converges faster/to better variance
    ceres_default_options.min_line_search_step_size                  = 1e-64;//  std::numeric_limits<double>::epsilon();
    ceres_default_options.max_line_search_step_contraction           = 1e-3;
    ceres_default_options.min_line_search_step_contraction           = 0.6;
    ceres_default_options.max_line_search_step_expansion             = 10;
    ceres_default_options.max_num_line_search_step_size_iterations   = 20;//20;
    ceres_default_options.max_num_line_search_direction_restarts     = 5;//2;
    ceres_default_options.line_search_sufficient_function_decrease   = 1e-4; //Tested, doesn't seem to matter between [1e-1 to 1e-4]. Default is fine: 1e-4
    ceres_default_options.line_search_sufficient_curvature_decrease  = 0.9; // This one should be above 0.5. Below, it makes retries at every step and starts taking twice as long for no added benefit. Tested 0.9 to be sweetspot
    ceres_default_options.max_solver_time_in_seconds                 = 60*10;//60*2;
    ceres_default_options.function_tolerance                         = 1e-5; // Tested, 1e-6 seems to be a sweetspot
    ceres_default_options.gradient_tolerance                         = 1e-4; // Not tested yet
    ceres_default_options.parameter_tolerance                        = 1e-10;
    ceres_default_options.minimizer_progress_to_stdout               = false; //tools::log->level() <= spdlog::level::trace;
    ceres_default_options.logging_type                               = ceres::LoggingType::PER_MINIMIZER_ITERATION;
    if(status.algorithm_has_got_stuck){
        ceres_default_options.max_num_iterations                        = 8000;
        ceres_default_options.function_tolerance                        = 1e-6;
        ceres_default_options.gradient_tolerance                        = 1e-6;
        ceres_default_options.parameter_tolerance                       = 1e-12;
        ceres_default_options.max_solver_time_in_seconds                = 60*20;//60*2;
        ceres_default_options.use_approximate_eigenvalue_bfgs_scaling   = true;  // True makes a huge difference, takes longer steps at each iteration!!
        ceres_default_options.minimizer_progress_to_stdout              = false;// tools::log->level() <= spdlog::level::debug;
    }
    if(status.algorithm_has_stuck_for > 1){
        ceres_default_options.function_tolerance                        = 1e-7;
        ceres_default_options.gradient_tolerance                        = 1e-6;
        ceres_default_options.parameter_tolerance                       = 1e-16;
        ceres_default_options.use_approximate_eigenvalue_bfgs_scaling   = true;  // True makes a huge difference, takes longer steps at each iteration. May not always be optimal according to library.
    }
    /* clang-format on */

    //    Progress log definitions:
    //    f is the value of the objective function.
    //    d is the change in the value of the objective function if the step computed in this iteration is accepted.
    //    g is the inf norm of the gradient (i.e. the largest element of |grad f| )
    //    h is the change in the parameter vector.
    //    s is the optimal step length computed by the line search.
    //    it is the time take by the current iteration.
    //    tt is the total time taken by the minimizer.
    opt_state result;
    switch(optSpace) {
            /* clang-format off */
        case OptSpace::SUBSPACE_ONLY:       result = internal::ceres_subspace_optimization(tensors,initial_tensor,status, optType, optMode,optSpace); break;
        case OptSpace::SUBSPACE_AND_DIRECT: result = internal::ceres_subspace_optimization(tensors,initial_tensor,status, optType, optMode,optSpace); break;
        case OptSpace::DIRECT:              result = internal::ceres_direct_optimization(tensors,initial_tensor,status, optType,optMode,optSpace); break;
            /* clang-format on */
    }
    tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt"]->toc();
    // Finish up and print reports
    reports::print_bfgs_report();
    reports::print_time_report();
    return result;
}

Eigen::Tensor<std::complex<double>, 3> tools::finite::opt::find_ground_state(const class_tensors_finite &tensors, StateRitz ritz) {
    return internal::ground_state_optimization(tensors, ritz);
}

double tools::finite::opt::internal::windowed_func_abs(double x, double window) {
    if(std::abs(x) >= window) return std::abs(x) - window;
    else
        return 0;
}
double tools::finite::opt::internal::windowed_grad_abs(double x, double window) {
    if(std::abs(x) >= window) return sgn(x);
    else
        return 0.0;
}

double tools::finite::opt::internal::windowed_func_pow(double x, double window) {
    if(std::abs(x) >= window) return x * x - window * window;
    else
        return 0.0;
}
double tools::finite::opt::internal::windowed_grad_pow(double x, double window) {
    if(std::abs(x) >= window) return 2.0 * x;
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
        //        func = x*x - window*window;
        //        grad = 2*x;
    }
    return std::make_pair(func, grad);
}


long tools::finite::opt::internal::get_ops_v1(long d, long chiL, long chiR, long m) {
    // d first
    // This is the one chosen by eigen whenever d > m
    long step1   = chiL * chiL * chiR * m * m * d;
    long step2_d = chiL * chiR * d * m * m * m * (d * m + 1);
    long step3   = chiL * chiR * d * m * m * m * (chiR * m + 1);
    long step4_d = chiL * chiR * d * m * (d * m * m * m + m * m + 1);
    return step1 + step2_d + step3 + step4_d;
}

long tools::finite::opt::internal::get_ops_v2(long d, long chiL, long chiR, long m) {
    // d first
    // This is the one chosen by eigen whenever d > m
    // Same as v1, just swap chiL and chiR
    return tools::finite::opt::internal::get_ops_v1(d,chiR,chiL,m);
}

long tools::finite::opt::internal::get_ops_v3(long d, long chiL, long chiR,long m) {
    // d first
    // This is the one chosen by eigen whenever d > m
    long step1   = chiL * chiL * chiR * m * m * d;
    long step2_d = chiL * chiR * d * m * m * m * (d * m + 1);
    long step3_d = chiL * chiR * d * m * m * m * (d * m + 1);
    long step4   = chiL * chiR * d * m * (chiR * m * m * m + m * m + 1);
    return step1 + step2_d + step3_d + step4;
}