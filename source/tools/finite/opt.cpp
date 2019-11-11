//
// Created by david on 2019-03-18.
//
#include <string>
#include <sstream>
#include <tools/finite/opt.h>
#include <state/class_state_finite.h>
#include <simulation/class_simulation_status.h>
#include <model/class_model_base.h>
#include <simulation/nmspc_settings.h>
#include <glog/logging.h>

Eigen::Tensor<class_state_finite::Scalar,3>
tools::finite::opt::find_excited_state(const class_state_finite &state, const class_simulation_status &sim_status, OptMode optMode, OptSpace optSpace, OptType optType){
    tools::log->trace("Optimizing excited state");
    using namespace opt::internal;
    using namespace Textra;
    static bool googleLogginghasInitialized = false;
    if(not googleLogginghasInitialized){
        google::InitGoogleLogging(tools::log->name().c_str());
        googleLogginghasInitialized = true;
    }

    std::stringstream problem_report;
    auto dims = state.active_dimensions();
    problem_report << fmt::format("Starting optimization: ")
                   << fmt::format("mode [ {} ] "             , optMode.str())
                   << fmt::format("space [ {} ] "            , optSpace.str())
                   << fmt::format("type [ {} ] "             , optType.str())
                   << fmt::format("position [ {} ] "         , state.get_position())
                   << fmt::format("shape [ {} {} {} ] = {} " , dims[0],dims[1],dims[2], state.active_problem_size());
    tools::log->debug(problem_report.str());


//    ceres::GradientProblemSolver::Options options;
    ceres_default_options.line_search_type = ceres::LineSearchType::WOLFE;
    ceres_default_options.line_search_interpolation_type = ceres::LineSearchInterpolationType::CUBIC;
    ceres_default_options.line_search_direction_type = ceres::LineSearchDirectionType::LBFGS;
    ceres_default_options.nonlinear_conjugate_gradient_type = ceres::NonlinearConjugateGradientType::POLAK_RIBIERE;
    ceres_default_options.max_num_iterations = 1000;
    ceres_default_options.max_lbfgs_rank     = 250;
    ceres_default_options.use_approximate_eigenvalue_bfgs_scaling = true;  // True makes a huge difference, takes longer steps at each iteration!!
    ceres_default_options.max_line_search_step_expansion = 1e4;// 100.0;
    ceres_default_options.min_line_search_step_size = 1e-128;//  std::numeric_limits<double>::epsilon();
    ceres_default_options.max_line_search_step_contraction = 1e-3;
    ceres_default_options.min_line_search_step_contraction = 0.6;
    ceres_default_options.max_num_line_search_step_size_iterations  = 40;//20;
    ceres_default_options.max_num_line_search_direction_restarts    = 5;//2;
    ceres_default_options.line_search_sufficient_function_decrease  = 1e-2;
    ceres_default_options.line_search_sufficient_curvature_decrease = 0.4; //0.5;
    ceres_default_options.max_solver_time_in_seconds = 60*5;//60*2;
    ceres_default_options.function_tolerance = 1e-5;
    ceres_default_options.gradient_tolerance = 1e-2;
    ceres_default_options.parameter_tolerance = 1e-128;//std::numeric_limits<double>::epsilon();//1e-12;
    ceres_default_options.minimizer_progress_to_stdout = tools::log->level() <= spdlog::level::trace;
    ceres_default_options.logging_type = ceres::LoggingType::SILENT;
    if(sim_status.simulation_has_got_stuck){
//        options.min_line_search_step_size = std::numeric_limits<double>::epsilon();
        ceres_default_options.function_tolerance = 1e-7; //Operations are cheap in subspace, so you can afford low tolerance
        ceres_default_options.max_num_iterations = 1000;
        ceres_default_options.gradient_tolerance = 1e-4;
        ceres_default_options.max_solver_time_in_seconds = 60*10;//60*2;
    }


//    Progress log definitions:
//    f is the value of the objective function.
//    d is the change in the value of the objective function if the step computed in this iteration is accepted.
//    g is the max norm of the gradient (i.e. the largest element of |grad f|)
//    h is the change in the parameter vector.
//    s is the optimal step length computed by the line search.
//    it is the time take by the current iteration.
//    tt is the total time taken by the minimizer.


    switch (optSpace.option){
        case opt::SPACE::SUBSPACE:    return internal::ceres_subspace_optimization(state,sim_status, optType, optMode);
        case opt::SPACE::DIRECT:      return internal::ceres_direct_optimization(state,sim_status, optType);
    }
}

Eigen::Tensor<std::complex<double>,4> tools::finite::opt::find_ground_state(
        const class_state_finite &state, std::string ritz){
    return internal::ground_state_optimization(state,ritz);

}





void tools::finite::opt::internal::reset_timers(){
    t_opt-> reset();
    t_eig-> reset();
    t_ham-> reset();
    t_tot-> reset();
    t_vH2v->reset();
    t_vHv ->reset();
    t_vH2 ->reset();
    t_vH  ->reset();
    t_op  ->reset();
}



//
//template<typename Scalar>
//tools::finite::opt::internal::MultiComponents<Scalar>::MultiComponents(const class_state_finite & state){
//    tools::log->trace("Generating multi components");
//    if constexpr (std::is_same<Scalar,double>::value){
//        mpo                          = state.get_multimpo().real();
//        auto & envL_cplx  = state.get_ENVL(state.active_sites.front());
//        auto & envR_cplx  = state.get_ENVR(state.active_sites.back());
//        auto & env2L_cplx = state.get_ENV2L(state.active_sites.front());
//        auto & env2R_cplx = state.get_ENV2R(state.active_sites.back());
//
//        envL  = envL_cplx.block.real();         envR  = envR_cplx.block.real();
//        env2L = env2L_cplx.block.real();        env2R = env2R_cplx.block.real();
//    }
//
//    if constexpr (std::is_same<Scalar,std::complex<double>>::value){
//        mpo                          = state.get_multimpo();
//        auto & envL_cplx  = state.get_ENVL(state.active_sites.front());
//        auto & envR_cplx  = state.get_ENVR(state.active_sites.back());
//        auto & env2L_cplx = state.get_ENV2L(state.active_sites.front());
//        auto & env2R_cplx = state.get_ENV2R(state.active_sites.back());
//        envL  = envL_cplx.block;         envR  = envR_cplx.block;
//        env2L = env2L_cplx.block;        env2R = env2R_cplx.block;
//    }
//
//
//
//    dsizes        = state.active_dimensions();
//    tools::log->trace("Finished building multicomponents");
//}
//
//template struct tools::finite::opt::internal::MultiComponents<double>;
//template struct tools::finite::opt::internal::MultiComponents<std::complex<double>>;


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
        func = x*x - window*window;
        grad = 2*x;
    }
    return std::make_pair(func,grad);
}

