//
// Created by david on 2019-07-09.
//

#include "ceres_direct_functor.h"
#include <general/class_tic_toc.h>
#include <state/class_state_finite.h>
#include <simulation/class_simulation_status.h>
#include <simulation/nmspc_settings.h>
#include <tools/finite/measure.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

Eigen::Tensor<class_state_finite::Scalar,3>
tools::finite::opt::internal::ceres_direct_optimization(const class_state_finite &state,
                                                         const class_simulation_status &sim_status,
                                                         OptType optType, OptMode optMode, OptSpace optSpace){
    return ceres_direct_optimization(state,state.get_multitheta(),sim_status,optType,optMode,optSpace);
}


Eigen::Tensor<class_state_finite::Scalar,3>
tools::finite::opt::internal::ceres_direct_optimization(const class_state_finite &state,
                                                         const Eigen::Tensor<class_state_finite::Scalar,3> & theta_initial,
                                                         const class_simulation_status &sim_status, OptType optType, OptMode optMode, OptSpace optSpace){
    tools::log->trace("Optimizing in DIRECT mode");
    tools::common::profile::t_opt->tic();
    using Scalar = std::complex<double>;
    auto & theta_old          = state.get_multitheta();
    auto   theta_old_vec      = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta_old.data(),theta_old.size());
    auto   theta_initial_vec  = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta_initial.data(), theta_initial.size());
    std::vector<reports::direct_opt_tuple> opt_log;
    class_tic_toc t_local(true,5,"");
    if (tools::log->level() <= spdlog::level::debug){
        t_local.tic();
        double energy_old   = tools::finite::measure::energy_per_site(state,theta_old);
        double variance_old = tools::finite::measure::energy_variance_per_site(state,theta_old);
        t_local.toc();
        opt_log.emplace_back("Current state" ,theta_old.size(), energy_old, std::log10(variance_old), 1.0, theta_old_vec.norm(), 0 ,0, t_local.get_last_time_interval());

        t_local.tic();
        double energy_initial   = tools::finite::measure::multisite::energy_per_site(state,theta_initial);
        double variance_initial = tools::finite::measure::multisite::energy_variance_per_site(state,theta_initial);
        t_local.toc();
        opt_log.emplace_back("Initial guess" , theta_initial.size(), energy_initial, std::log10(variance_initial), 1.0, theta_initial_vec.norm(), 0 , 0, t_local.get_last_time_interval());
    }

    auto options = ceres_default_options;
    ceres::GradientProblemSolver::Summary summary;
    int counter,iter;
    double variance_ceres, energy_ceres;
    t_local.tic();
    Eigen::VectorXcd theta_new;
    switch (optType){
        case OptType::CPLX:{
            Eigen::VectorXd  theta_start_cast = Eigen::Map<const Eigen::VectorXd>(reinterpret_cast<const double*> (theta_initial_vec.data()), 2 * theta_initial_vec.size());
            auto * functor = new ceres_direct_functor<std::complex<double>>(state, sim_status);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running L-BFGS");
            ceres::Solve(options, problem, theta_start_cast.data(), &summary);
            iter         = (int)summary.iterations.size();
            counter      = functor->get_count();
            variance_ceres = functor->get_variance();
            energy_ceres = functor->get_energy();
            theta_new    = Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<Scalar*> (theta_start_cast.data()), theta_start_cast.size() / 2).normalized();
            break;
        }
        case OptType::REAL:{
            Eigen::VectorXd  theta_start_cast = theta_initial_vec.real();
            auto * functor = new ceres_direct_functor<double>(state, sim_status);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running L-BFGS");
            ceres::Solve(options, problem, theta_start_cast.data(), &summary);
            iter        = (int)summary.iterations.size();
            counter      = functor->get_count();
            variance_ceres = functor->get_variance();
            energy_ceres = functor->get_energy();
            theta_new    = theta_start_cast.normalized().cast<Scalar>();
            break;
        }
    }
    t_local.toc();

    double overlap_new  = std::abs(theta_old_vec.dot(theta_new));

    opt_log.emplace_back("LBFGS direct", theta_new.size(), energy_ceres, std::log10(variance_ceres), overlap_new, theta_new.norm(), iter, counter,t_local.get_last_time_interval());
    if (tools::log->level() <= spdlog::level::debug) {
        t_local.tic();
        auto theta_new_map  = Textra::MatrixTensorMap(theta_new, state.active_dimensions());
        double energy_new   = tools::finite::measure::energy_per_site(state,theta_new_map);
        double variance_new = tools::finite::measure::energy_variance_per_site(state,theta_new_map);
        t_local.toc();
        opt_log.emplace_back("LBFGS check", theta_new.size(), energy_new, std::log10(variance_new), overlap_new, theta_new.norm(), 0,0,
                             t_local.get_last_time_interval());
    }
    // Finish up and print reports
    tools::log->debug("Finished LBFGS after {} seconds ({} iters). Exit status: {}. Message: {}",summary.total_time_in_seconds, summary.iterations.size(), ceres::TerminationTypeToString(summary.termination_type) , summary.message.c_str());
    if(optSpace == OptSpace::DIRECT){
        reports::print_report(opt_log);
        reports::print_report(std::make_tuple(
            tools::common::profile::t_vH2v->get_measured_time(),
            tools::common::profile::t_vHv->get_measured_time(),
            tools::common::profile::t_vH2->get_measured_time(),
            tools::common::profile::t_vH->get_measured_time(),
            tools::common::profile::t_op->get_measured_time()
        ));
    }


    tools::common::profile::t_opt->toc();

    tools::log->trace("Returning theta from optimization mode {} space {}",optMode,optSpace);
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


