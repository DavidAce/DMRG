//
// Created by david on 2019-07-09.
//

#include "ceres_direct_functor.h"
#include <general/class_tic_toc.h>
#include <state/class_state_finite.h>
#include <simulation/class_simulation_status.h>
#include <simulation/nmspc_settings.h>
#include <ceres/ceres.h>



Eigen::Tensor<class_state_finite::Scalar,3>
tools::finite::opt::internal::ceres_direct_optimization(const class_state_finite &state,
                                                         const class_simulation_status &sim_status,
                                                         OptType optType){
    return ceres_direct_optimization(state,state.get_multitheta(),sim_status,optType);
}


Eigen::Tensor<class_state_finite::Scalar,3>
tools::finite::opt::internal::ceres_direct_optimization(const class_state_finite &state,
                                                         const Eigen::Tensor<class_state_finite::Scalar,3> & theta_initial,
                                                         const class_simulation_status &sim_status, OptType optType){
    tools::log->trace("Optimizing in DIRECT mode");
    tools::common::profile::t_opt.tic();
    using Scalar = std::complex<double>;
    auto & theta_old          = state.get_multitheta();
    auto   theta_old_vec      = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta_old.data(),theta_old.size());
    auto   theta_initial_vec  = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta_initial.data(), theta_initial.size());
    std::vector<reports::direct_opt_tuple> opt_log;

    if (tools::log->level() <= spdlog::level::debug){
        t_opt->tic();
        double energy_old   = tools::finite::measure::energy_per_site(state);
        double variance_old = tools::finite::measure::energy_variance_per_site(state);
        t_opt->toc();
        opt_log.emplace_back("Current state" ,theta_old.size(), energy_old, std::log10(variance_old), 1.0, theta_old_vec.norm(), 0 ,0,t_opt->get_last_time_interval());

        t_opt->tic();
        double energy_initial   = tools::finite::measure::multisite::energy_per_site(state,theta_initial);
        double variance_initial = tools::finite::measure::multisite::energy_variance_per_site(state,theta_initial);
        t_opt->toc();
        opt_log.emplace_back("Initial guess" , theta_initial.size(), energy_initial, std::log10(variance_initial), 1.0, theta_initial_vec.norm(), 0 , 0, t_opt->get_last_time_interval());
    }
    double energy_new = 0,variance_new = 0,overlap_new = 0;
    auto options = ceres_default_options;

    ceres::GradientProblemSolver::Summary summary;
    int counter,iter;
    t_opt->tic();
    Eigen::VectorXcd theta_new;
    switch (optType.option){
        case opt::TYPE::CPLX:{
            Eigen::VectorXd  theta_start_cast = Eigen::Map<const Eigen::VectorXd>(reinterpret_cast<const double*> (theta_initial_vec.data()), 2 * theta_initial_vec.size());
            auto * functor = new ceres_direct_functor<std::complex<double>>(state, sim_status);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running L-BFGS");
            ceres::Solve(options, problem, theta_start_cast.data(), &summary);

            iter         = (int)summary.iterations.size();
            counter      = functor->get_count();
            energy_new   = functor->get_energy() ;
            variance_new = functor->get_variance();
            theta_new    = Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<Scalar*> (theta_start_cast.data()), theta_start_cast.size() / 2).normalized();
            break;
        }
        case opt::TYPE::REAL:{
            Eigen::VectorXd  theta_start_cast = theta_initial_vec.real();
            auto * functor = new ceres_direct_functor<double>(state, sim_status);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running L-BFGS");
            ceres::Solve(options, problem, theta_start_cast.data(), &summary);
            iter        = (int)summary.iterations.size();
            counter      = functor->get_count();
            energy_new   = functor->get_energy();
            variance_new = functor->get_variance();
            theta_new    = theta_start_cast.normalized().cast<Scalar>();
            break;
        }
    }
    t_opt->toc();

    if (tools::log->level() <= spdlog::level::debug){

//        auto theta_old = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size());
        overlap_new  = std::abs(theta_old_vec.dot(theta_new));
        opt_log.emplace_back("Ceres L-BFGS", theta_new.size(), energy_new, std::log10(variance_new), overlap_new, theta_new.norm(), iter, counter, t_opt->get_last_time_interval());

        // Sanity check
        t_opt->tic();
        auto theta_san      = Textra::MatrixTensorMap(theta_new, state.active_dimensions());
        double energy_san   = tools::finite::measure::multisite::energy_per_site(state,theta_san);
        double variance_san = tools::finite::measure::multisite::energy_variance_per_site(state,theta_san);
        t_opt->toc();
        opt_log.emplace_back("Sanity check", theta_san.size(), energy_san, std::log10(variance_san), overlap_new, theta_new.norm(), 0, 0, t_opt->get_last_time_interval());

        //double variance_acc = tools::finite::measure::reduced::energy_variance_per_site(state,theta_san);
        //opt_log.emplace_back("Sanity check (reduced)",theta_san.size(), energy_san, std::log10(variance_acc), overlap_new, theta_initial_vec.norm(), 0,0, t_opt->get_last_time_interval());


    }

    // Finish up and print reports
    tools::log->debug("Finished LBFGS after {} seconds ({} iters). Exit status: {}. Message: {}",summary.total_time_in_seconds, summary.iterations.size(), ceres::TerminationTypeToString(summary.termination_type) , summary.message.c_str());
//    tools::log->trace("Finished Ceres. Exit status: {}. Message: {}", ceres::TerminationTypeToString(summary.termination_type) , summary.message.c_str());
//    std::cout << summary.FullReport() << "\n";
    reports::print_report(opt_log);
    reports::print_report(std::make_tuple(
            tools::finite::opt::internal::t_vH2v->get_measured_time(),
            tools::finite::opt::internal::t_vHv->get_measured_time(),
            tools::finite::opt::internal::t_vH2->get_measured_time(),
            tools::finite::opt::internal::t_vH->get_measured_time(),
            tools::finite::opt::internal::t_op->get_measured_time()
    ));

    tools::common::profile::t_opt.toc();

    tools::log->debug("Returning new theta from DIRECT optimization");
    return  Textra::MatrixTensorMap(theta_new, state.active_dimensions());



//    if (variance_new < 1.0 * tools::finite::measure::energy_variance_per_site(state)){
//        // Only an improvement of 1% is considered to be an actual improvement
//        tools::log->debug("Returning new (better) theta");
//        state.tag_active_sites_have_been_updated(true);
//        return  Textra::Matrix_to_Tensor(theta_new, state.active_dimensions());
//
//    }
//    else if (variance_new < 100.0 * tools::finite::measure::energy_variance_per_site(state)) {
//        // Allow for variance to increase a bit to come out of local minima
//        tools::log->debug("Returning new (worse) theta");
//        state.tag_active_sites_have_been_updated(false);
//        return  Textra::Matrix_to_Tensor(theta_new, state.active_dimensions());
//    }
//    else{
//        tools::log->debug("Direct optimization didn't improve variance.");
//        tools::log->debug("Returning old theta");
//        state.tag_active_sites_have_been_updated(variance_new <= settings::precision::variance_convergence_threshold);
//        return  theta_old;
//    }

}


