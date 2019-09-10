//
// Created by david on 2019-07-09.
//

#include "ceres_direct_functor.h"
#include <spdlog/spdlog.h>
#include <general/class_tic_toc.h>
#include <state/class_finite_state.h>
#include <state/class_environment.h>
#include <simulation/class_simulation_status.h>
#include <general/nmspc_random_numbers.h>
#include <simulation/nmspc_settings.h>
#include <ceres/ceres.h>


Eigen::Tensor<std::complex<double>,3>
tools::finite::opt::internals::ceres_direct_optimization(const class_finite_state &state,
                                                         const class_simulation_status &sim_status, OptType optType){
//    opt::internals::old_direct_optimization(state,sim_status,optType);

    tools::log->trace("Optimizing in CERES mode");
    tools::common::profile::t_opt.tic();
    using Scalar = std::complex<double>;
    t_opt->tic();
    auto theta = state.get_multitheta();
    Eigen::VectorXcd theta_start  = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size());
    std::vector<reports::direct_opt_tuple> opt_log;
    t_opt->toc();
    if (tools::log->level() <= spdlog::level::debug){
        double energy_0   = tools::finite::measure::multisite::energy_per_site(state,theta);
        double variance_0 = tools::finite::measure::multisite::energy_variance_per_site(state,theta);
        double variance_r = tools::finite::measure::reduced::energy_variance_per_site(state,theta);
        t_opt->toc();
        opt_log.emplace_back("Initial"          ,theta_start.size(), energy_0, std::log10(variance_0), 1.0, theta_start.norm(), 0 ,0,t_opt->get_last_time_interval());
        opt_log.emplace_back("Initial (reduced)",theta_start.size(), energy_0, std::log10(variance_r), 1.0, theta_start.norm(), 0 ,0,t_opt->get_last_time_interval());
    }else{
        t_opt->toc();
    }
    double energy_new,variance_new,overlap_new;

//    ceres::GradientProblemSolver::Options options;
//    options.line_search_type = ceres::LineSearchType::WOLFE;
//    options.line_search_interpolation_type = ceres::LineSearchInterpolationType::CUBIC;
//    options.line_search_direction_type = ceres::LineSearchDirectionType::LBFGS;
//    options.nonlinear_conjugate_gradient_type = ceres::NonlinearConjugateGradientType::POLAK_RIBIERE;
//    options.minimizer_progress_to_stdout = Logger::getLogLevel(tools::log) < 2 ? true : false;
//    options.max_num_iterations = 150;
//    options.max_lbfgs_rank     = 250;
//    options.use_approximate_eigenvalue_bfgs_scaling = true;
//    options.max_line_search_step_expansion = 100.0;
//    options.min_line_search_step_size = 1e-8;
//    options.max_line_search_step_contraction = 1e-3;
//    options.min_line_search_step_contraction = 0.6;
//    options.max_num_line_search_step_size_iterations  = 20;
//    options.max_num_line_search_direction_restarts    = 2;
//    options.line_search_sufficient_function_decrease  = 1e-2;
//    options.line_search_sufficient_curvature_decrease = 0.5;
//    options.max_solver_time_in_seconds = 60*2;
//    options.function_tolerance = 1e-4;
//    options.gradient_tolerance = 1e-8;
//    options.parameter_tolerance = 1e-12;



    ceres::GradientProblemSolver::Options options;
    options.line_search_type = ceres::LineSearchType::WOLFE;
    options.line_search_interpolation_type = ceres::LineSearchInterpolationType::CUBIC;
    options.line_search_direction_type = ceres::LineSearchDirectionType::LBFGS;
    options.nonlinear_conjugate_gradient_type = ceres::NonlinearConjugateGradientType::POLAK_RIBIERE;
    options.max_num_iterations = 300;
    options.max_lbfgs_rank     = 250;
    options.use_approximate_eigenvalue_bfgs_scaling = true;
    options.max_line_search_step_expansion = 100;// 100.0;
    options.min_line_search_step_size = 1e-9;
    options.max_line_search_step_contraction = 1e-3;
    options.min_line_search_step_contraction = 0.6;
    options.max_num_line_search_step_size_iterations  = 30;//20;
    options.max_num_line_search_direction_restarts    = 5;//2;
    options.line_search_sufficient_function_decrease  = 1e-2;// 1e-2;
    options.line_search_sufficient_curvature_decrease = 0.5; //0.5;
    options.max_solver_time_in_seconds = 60*5;//60*2;
    options.function_tolerance = 1e-4;
    options.gradient_tolerance = 1e-8;
    options.parameter_tolerance = 1e-12;//1e-12;

    options.minimizer_progress_to_stdout = tools::log->level() <= spdlog::level::debug;

    ceres::GradientProblemSolver::Summary summary;
    int counter,iter;
    t_opt->tic();
    switch (optType){
        case OptType::CPLX:{
            Eigen::VectorXd  theta_start_cast = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double*> (theta_start.data()), 2*theta_start.size());
            auto * functor = new ceres_direct_functor<std::complex<double>>(state, sim_status);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running L-BFGS");
            ceres::Solve(options, problem, theta_start_cast.data(), &summary);

            iter         = (int)summary.iterations.size();
            counter      = functor->get_count();
            energy_new   = functor->get_energy() ;
            variance_new = functor->get_variance();
            theta_start  = Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<Scalar*> (theta_start_cast.data()), theta_start_cast.size()/2).normalized();
//            delete functor;
            break;
        }
        case OptType::REAL:{
            Eigen::VectorXd  theta_start_cast = theta_start.real();
            auto * functor = new ceres_direct_functor<double>(state, sim_status);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running L-BFGS");
            ceres::Solve(options, problem, theta_start_cast.data(), &summary);
            iter        = (int)summary.iterations.size();
            counter      = functor->get_count();
            energy_new   = functor->get_energy();
            variance_new = functor->get_variance();
            theta_start  = theta_start_cast.normalized().cast<Scalar>();
//            delete functor;
            break;
        }
    }
    t_opt->toc();

    if (tools::log->level() <= spdlog::level::debug){

        auto theta_old = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size());
        overlap_new  = std::abs(theta_old.dot(theta_start));
        opt_log.emplace_back("Ceres L-BFGS",theta.size(), energy_new, std::log10(variance_new), overlap_new, theta_start.norm(), iter,counter, t_opt->get_last_time_interval());

        // Sanity check
        t_opt->tic();
        auto theta_san      = Textra::Matrix_to_Tensor(theta_start, state.active_dimensions());
        double energy_san   = tools::finite::measure::multisite::energy_per_site(state,theta_san);
        double variance_san = tools::finite::measure::multisite::energy_variance_per_site(state,theta_san);
        t_opt->toc();
        opt_log.emplace_back("Sanity check",theta_san.size(), energy_san, std::log10(variance_san), overlap_new, theta_start.norm(), 0,0, t_opt->get_last_time_interval());

        double variance_acc = tools::finite::measure::reduced::energy_variance_per_site(state,theta_san);
        opt_log.emplace_back("Sanity check (reduced)",theta_san.size(), energy_san, std::log10(variance_acc), overlap_new, theta_start.norm(), 0,0, t_opt->get_last_time_interval());


    }

    // Finish up and print reports
    tools::log->trace("Finished Ceres. Exit status: {}. Message: {}", ceres::TerminationTypeToString(summary.termination_type) , summary.message.c_str());
//    std::cout << summary.FullReport() << "\n";
    reports::print_report(opt_log);
    reports::print_report(std::make_tuple(
            tools::finite::opt::internals::t_vH2v->get_measured_time(),
            tools::finite::opt::internals::t_vHv->get_measured_time(),
            tools::finite::opt::internals::t_vH2->get_measured_time(),
            tools::finite::opt::internals::t_vH->get_measured_time(),
            tools::finite::opt::internals::t_op->get_measured_time()
    ));

    // Return something strictly better
    tools::common::profile::t_opt.toc();
    if (variance_new < 1e-11 or variance_new < tools::finite::measure::multisite::energy_variance_per_site(state,theta)){
        state.unset_measurements();
        tools::log->debug("Returning new theta");
        state.tag_active_sites_have_been_updated(true);
        return  Textra::Matrix_to_Tensor(theta_start, state.active_dimensions());

    }else{
        state.tag_active_sites_have_been_updated(false);
        tools::log->debug("Returning old theta");
        return  theta;
    }

}


