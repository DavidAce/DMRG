//
// Created by david on 2019-03-18.
//

#include <glog/logging.h>
#include <tensors/model/class_mpo_base.h>
#include <algorithms/class_algorithm_status.h>
#include <tensors/state/class_state_finite.h>
#include <tensors/class_tensors_finite.h>
#include <string>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt.h>

Eigen::Tensor<class_state_finite::Scalar,3>
tools::finite::opt::find_excited_state(const class_tensors_finite &tensors, const class_algorithm_status &status, OptMode optMode, OptSpace optSpace, OptType optType){
    tools::log->trace("Optimizing excited state");
    using namespace opt::internal;
    using namespace Textra;
    static bool googleLogginghasInitialized = false;
    if(not googleLogginghasInitialized){
        google::InitGoogleLogging(tools::log->name().c_str());
        googleLogginghasInitialized = true;
        google::SetStderrLogging(3);
    }

    std::stringstream problem_report;
    auto dims = tensors.state->active_dimensions();
    problem_report << fmt::format("Starting optimization: ")
                   << fmt::format("mode [ {} ] "             , optMode)
                   << fmt::format("space [ {} ] "            , optSpace)
                   << fmt::format("type [ {} ] "             , optType)
                   << fmt::format("position [ {} ] "         , tensors.get_position())
                   << fmt::format("shape [ {} {} {} ] = {} " , dims[0],dims[1],dims[2], tensors.state->active_problem_size());
    tools::log->debug(problem_report.str());


    ceres_default_options.line_search_type = ceres::LineSearchType::WOLFE;
    ceres_default_options.line_search_interpolation_type = ceres::LineSearchInterpolationType::CUBIC;
    ceres_default_options.line_search_direction_type = ceres::LineSearchDirectionType::LBFGS;
    ceres_default_options.nonlinear_conjugate_gradient_type = ceres::NonlinearConjugateGradientType::FLETCHER_REEVES;
    ceres_default_options.max_num_iterations = 2000;
    ceres_default_options.max_lbfgs_rank     = 16; // Tested: around 8-32 seems to be a good compromise, anything larger incurs a large overhead. The overhead means 2x computation time at ~256
    ceres_default_options.use_approximate_eigenvalue_bfgs_scaling = true;  // Tested: True makes a huge difference, takes longer steps at each iteration and generally converges faster/to better variance
    ceres_default_options.min_line_search_step_size = 1e-64;//  std::numeric_limits<double>::epsilon();
    ceres_default_options.max_line_search_step_contraction = 1e-3;
    ceres_default_options.min_line_search_step_contraction = 0.6;
    ceres_default_options.max_line_search_step_expansion = 10;
    ceres_default_options.max_num_line_search_step_size_iterations  = 40;//20;
    ceres_default_options.max_num_line_search_direction_restarts    = 50;//2;
    ceres_default_options.line_search_sufficient_function_decrease  = 1e-3; //Tested, doesn't seem to matter between [1e-1 to 1e-4]. Default is fine: 1e-4
    ceres_default_options.line_search_sufficient_curvature_decrease = 0.8; //This one should be above 0.5. Below, it makes retries at every step and starts taking twice as long for no added benefit.
    ceres_default_options.max_solver_time_in_seconds = 60*10;//60*2;
    ceres_default_options.function_tolerance = 1e-6; // Tested, 1e-6 seems to be a sweetspot
    ceres_default_options.gradient_tolerance = 1e-4; // Not tested yet
    ceres_default_options.parameter_tolerance = 1e-6;
    ceres_default_options.minimizer_progress_to_stdout = tools::log->level() <= spdlog::level::trace;
    ceres_default_options.logging_type = ceres::LoggingType::PER_MINIMIZER_ITERATION;


    if(status.algorithm_has_got_stuck){
        ceres_default_options.max_num_iterations  = 4000;
        ceres_default_options.function_tolerance  = 1e-8;
        ceres_default_options.gradient_tolerance  = 1e-6;
        ceres_default_options.parameter_tolerance = 1e-8;
        ceres_default_options.max_solver_time_in_seconds = 60*20;//60*2;
        ceres_default_options.use_approximate_eigenvalue_bfgs_scaling = true;  // True makes a huge difference, takes longer steps at each iteration!!
        ceres_default_options.minimizer_progress_to_stdout = tools::log->level() <= spdlog::level::debug;
    }
    if(status.algorithm_has_stuck_for > 1){
        ceres_default_options.use_approximate_eigenvalue_bfgs_scaling = true;  // True makes a huge difference, takes longer steps at each iteration. May not always be optimal according to library.
    }


//    Progress log definitions:
//    f is the value of the objective function.
//    d is the change in the value of the objective function if the step computed in this iteration is accepted.
//    g is the inf norm of the gradient (i.e. the largest element of |grad f| )
//    h is the change in the parameter vector.
//    s is the optimal step length computed by the line search.
//    it is the time take by the current iteration.
//    tt is the total time taken by the minimizer.


    switch (optSpace){
        case OptSpace::SUBSPACE_ONLY:    return internal::ceres_subspace_optimization(tensors,status, optType, optMode,optSpace);
        case OptSpace::SUBSPACE_AND_DIRECT:  return internal::ceres_subspace_optimization(tensors,status, optType, optMode,optSpace);
        case OptSpace::DIRECT:      return internal::ceres_direct_optimization(tensors,status, optType,optMode,optSpace);
    }
    throw std::logic_error("No valid optimization type given");
}

Eigen::Tensor<std::complex<double>,4> tools::finite::opt::find_ground_state(
        const class_tensors_finite &tensors,  StateRitz ritz){
    return internal::ground_state_optimization(tensors,ritz);
}





void tools::finite::opt::internal::reset_timers(){
//    tools::common::profile::t_opt-> reset();
//    tools::common::profile::t_eig-> reset();
//    tools::common::profile::t_tot-> reset();
    tools::common::profile::t_ham-> reset();
    tools::common::profile::t_vH2v->reset();
    tools::common::profile::t_vHv ->reset();
    tools::common::profile::t_vH2 ->reset();
    tools::common::profile::t_vH  ->reset();
    tools::common::profile::t_op  ->reset();
}


double tools::finite::opt::internal::windowed_func_abs(double x,double window){
    if (std::abs(x) >= window){
        return std::abs(x)-window;
    }else{
        return 0;
    }
}
double tools::finite::opt::internal::windowed_grad_abs(double x,double window){
    if (std::abs(x) >= window){
        return sgn(x);
    }else{
        return 0.0;
    }
}



double tools::finite::opt::internal::windowed_func_pow(double x,double window){
    if (std::abs(x) >= window){
        return x*x - window*window;
    }else{
        return 0.0;
    }
}
double tools::finite::opt::internal::windowed_grad_pow(double x,double window){
    if (std::abs(x) >= window){
        return 2.0*x;
    }else{
        return 0.0;
    }
}



std::pair<double,double> tools::finite::opt::internal::windowed_func_grad(double x,double window){
    double func = 0;
    double grad = 0;
    if (std::abs(x) >= window){
        if(x > 0){
            func = (x-window)*(x-window);
            grad = 2*(x-window);
        }else{
            func = (x+window)*(x+window);
            grad = 2*(x+window);
        }
//        func = x*x - window*window;
//        grad = 2*x;
    }
    return std::make_pair(func,grad);
}

