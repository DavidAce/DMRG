//
// Created by david on 2019-07-09.
//

#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include "ceres_direct_functor.h"
#include <general/class_tic_toc.h>
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <tensors/state/class_state_finite.h>
#include <tensors/class_tensors_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>

Eigen::Tensor<class_state_finite::Scalar, 3> tools::finite::opt::internal::ceres_direct_optimization(const class_tensors_finite &tensors,
                                                                                                     const class_algorithm_status &status, OptType optType,
                                                                                                     OptMode optMode, OptSpace optSpace) {
    return ceres_direct_optimization(tensors, tensors.state->get_multisite_tensor(), status, optType, optMode, optSpace);
}

Eigen::Tensor<class_state_finite::Scalar, 3>
    tools::finite::opt::internal::ceres_direct_optimization(const class_tensors_finite &tensors, const Eigen::Tensor<class_state_finite::Scalar, 3> &theta_initial,
                                                            const class_algorithm_status &status, OptType optType, OptMode optMode, OptSpace optSpace) {
    const auto & state = *tensors.state;
    const auto & model = *tensors.model;
    const auto & edges = *tensors.edges;

    tools::log->trace("Optimizing in DIRECT mode");
    tools::common::profile::t_opt->tic();
    using Scalar                    = std::complex<double>;
    auto &        theta_old         = state.get_multisite_tensor();
    auto          theta_old_vec     = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(theta_old.data(), theta_old.size());
    auto          theta_initial_vec = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(theta_initial.data(), theta_initial.size());
    class_tic_toc t_local(true, 5, "");
    if(tools::log->level() <= spdlog::level::debug and optSpace == OptSpace::DIRECT) {
        t_local.tic();
        double energy_old   = tools::finite::measure::energy_per_site(tensors);
        double variance_old = tools::finite::measure::energy_variance_per_site(tensors);
        t_local.toc();
        reports::bfgs_log.emplace_back("Current state", theta_old.size(), 0, energy_old, std::log10(variance_old), 1.0, theta_old_vec.norm(), 1, 1,
                                       t_local.get_last_time_interval());
    }
    if(tools::log->level() <= spdlog::level::trace and optSpace == OptSpace::DIRECT and settings::debug) {
        t_local.tic();
        double energy_initial   = tools::finite::measure::energy_per_site(theta_initial,model,edges);
        double variance_initial = tools::finite::measure::energy_variance_per_site(theta_initial,model,edges);
        t_local.toc();
        reports::bfgs_log.emplace_back("Initial guess", theta_initial.size(), 0, energy_initial, std::log10(variance_initial), 1.0, theta_initial_vec.norm(), 1,
                                       1, t_local.get_last_time_interval());
    }

    std::vector<double>                                      bfgs_tols = {0.8};
    std::vector<std::pair<double, Eigen::Tensor<Scalar, 3>>> optimized_results;
    for(auto &bfgs_tol : bfgs_tols) {
        auto options                                      = ceres_default_options;
        options.line_search_sufficient_curvature_decrease = bfgs_tol;
        ceres::GradientProblemSolver::Summary summary;
        size_t                                counter = 0, iter = 0;
        double                                variance_ceres, energy_ceres;
        t_local.tic();
        Eigen::VectorXcd theta_new;
        switch(optType) {
            case OptType::CPLX: {
                Eigen::VectorXd theta_start_cast =
                    Eigen::Map<const Eigen::VectorXd>(reinterpret_cast<const double *>(theta_initial_vec.data()), 2 * theta_initial_vec.size());
                auto *                 functor = new ceres_direct_functor<std::complex<double>>(tensors, status);
                ceres::GradientProblem problem(functor);
                tools::log->trace("Running L-BFGS");
                ceres::Solve(options, problem, theta_start_cast.data(), &summary);
                iter           = summary.iterations.size();
                counter        = functor->get_count();
                variance_ceres = functor->get_variance();
                energy_ceres   = functor->get_energy();
                theta_new      = Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<Scalar *>(theta_start_cast.data()), theta_start_cast.size() / 2).normalized();
                break;
            }
            case OptType::REAL: {
                Eigen::VectorXd        theta_start_cast = theta_initial_vec.real();
                auto *                 functor          = new ceres_direct_functor<double>(tensors, status);
                ceres::GradientProblem problem(functor);
                tools::log->trace("Running L-BFGS");
                ceres::Solve(options, problem, theta_start_cast.data(), &summary);
                iter           = summary.iterations.size();
                counter        = functor->get_count();
                variance_ceres = functor->get_variance();
                energy_ceres   = functor->get_energy();
                theta_new      = theta_start_cast.normalized().cast<Scalar>();
                break;
            }
        }
        t_local.toc();

        double overlap_new = std::abs(theta_old_vec.dot(theta_new));

        reports::bfgs_log.emplace_back("LBFGS direct", theta_new.size(), options.max_lbfgs_rank, energy_ceres, std::log10(variance_ceres), overlap_new,
                                       theta_new.norm(), iter, counter, t_local.get_last_time_interval());
        optimized_results.emplace_back(std::make_pair(variance_ceres, Textra::MatrixTensorMap(theta_new, state.active_dimensions())));
        if(tools::log->level() <= spdlog::level::trace and optSpace == OptSpace::DIRECT and settings::debug) {
            t_local.tic();
            auto   theta_new_map = Textra::MatrixTensorMap(theta_new, state.active_dimensions());
            double energy_new    = tools::finite::measure::energy_per_site(theta_new_map,tensors);
            double variance_new  = tools::finite::measure::energy_variance_per_site(theta_new_map, tensors);
            t_local.toc();
            reports::bfgs_log.emplace_back("LBFGS check", theta_new.size(), 0, energy_new, std::log10(variance_new), overlap_new, theta_new.norm(), 1, 1,
                                           t_local.get_last_time_interval());
        }
        tools::log->debug("Finished LBFGS after {} seconds ({} iters). Exit status: {}. Message: {}", summary.total_time_in_seconds, summary.iterations.size(),
                          ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
        reports::time_log.emplace_back(tools::common::profile::t_vH2v->get_measured_time(), tools::common::profile::t_vHv->get_measured_time(),
                                       tools::common::profile::t_vH2->get_measured_time(), tools::common::profile::t_vH->get_measured_time(),
                                       tools::common::profile::t_op->get_measured_time());
        tools::common::profile::t_vH2v->reset();
        tools::common::profile::t_vHv->reset();
        tools::common::profile::t_vH2->reset();
        tools::common::profile::t_vH->reset();
        tools::common::profile::t_op->reset();
    }

    tools::common::profile::t_opt->toc();

    // Finish up and print reports
    if(optSpace == OptSpace::DIRECT) {
        reports::print_bfgs_report();
        reports::print_time_report();
    }

    tools::log->trace("Returning theta from optimization mode {} space {}", optMode, optSpace);
    // Sort thetas in ascending order in variance
    std::sort(optimized_results.begin(), optimized_results.end(), [](auto &left, auto &right) { return left.first < right.first; });
    // Return the best theta
    return optimized_results.front().second;

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
