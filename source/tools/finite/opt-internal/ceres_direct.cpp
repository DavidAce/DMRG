//
// Created by david on 2019-07-09.
//

#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include "ceres_direct_functor.h"
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <general/class_tic_toc.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt_tensor.h>

tools::finite::opt::opt_tensor tools::finite::opt::internal::ceres_direct_optimization(const class_tensors_finite &  tensors,
                                                                                       const class_algorithm_status &status, OptType optType, OptMode optMode,
                                                                                       OptSpace optSpace) {
    std::vector<size_t> sites (tensors.active_sites.begin(),tensors.active_sites.end());
    opt_tensor initial_tensor("current state", tensors.state->get_multisite_tensor(), sites,
                                       tools::finite::measure::energy(tensors) - tensors.model->get_energy_reduced(), // Eigval
                                       tensors.model->get_energy_reduced(), // Energy reduced for full system
                                       tools::finite::measure::energy_variance(tensors),
                                       1.0, // Overlap
                                       tensors.get_length());


    return ceres_direct_optimization(tensors, initial_tensor, status, optType, optMode, optSpace);
}

tools::finite::opt::opt_tensor tools::finite::opt::internal::ceres_direct_optimization(const class_tensors_finite &tensors, const opt_tensor &initial_tensor,
                                                                                       const class_algorithm_status &status, OptType optType, OptMode optMode,
                                                                                       OptSpace optSpace) {
    tools::log->trace("Optimizing in DIRECT mode");
    tools::common::profile::t_opt_dir->tic();

    reports::bfgs_add_entry("Direct", "init", initial_tensor);
    const auto & current_tensor = tensors.state->get_multisite_tensor();
    const auto   current_vector = Eigen::Map<const Eigen::VectorXcd>(current_tensor.data(),current_tensor.size());
    auto       options = ceres_default_options;
    auto       summary = ceres::GradientProblemSolver::Summary();
    opt_tensor optimized_tensor;
    optimized_tensor.set_name(initial_tensor.get_name());
    optimized_tensor.set_sites(initial_tensor.get_sites());
    optimized_tensor.set_length(initial_tensor.get_length());
    tools::common::profile::t_opt_dir_bfgs->tic();
    switch(optType) {
        case OptType::CPLX: {
            // Copy the initial guess and operate directly on it
            optimized_tensor.set_tensor(initial_tensor.get_tensor());
            auto *                 functor = new ceres_direct_functor<std::complex<double>>(tensors, status);

            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS direct cplx");
            ceres::Solve(options, problem, optimized_tensor.get_vector_cplx_as_2xreal().data(), &summary);
            // Copy the results from the functor
            optimized_tensor.set_counter(functor->get_count());
            *tools::common::profile::t_opt_dir_vH2 += *functor->t_vH2;
            *tools::common::profile::t_opt_dir_vH2v += *functor->t_vH2v;
            *tools::common::profile::t_opt_dir_vH += *functor->t_vH;
            *tools::common::profile::t_opt_dir_vHv += *functor->t_vHv;
            break;
        }
        case OptType::REAL: {
            // Here we make a temporary
            auto                   initial_tensor_real = initial_tensor.get_vector_cplx_as_1xreal();
            auto *                 functor             = new ceres_direct_functor<double>(tensors, status);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS direct real");
            ceres::Solve(options, problem, initial_tensor_real.data(), &summary);
            // Copy the results from the functor
            optimized_tensor.set_counter(functor->get_count());
            optimized_tensor.set_tensor_real(initial_tensor_real.data(), initial_tensor.get_tensor().dimensions());
            *tools::common::profile::t_opt_dir_vH2 += *functor->t_vH2;
            *tools::common::profile::t_opt_dir_vH2v += *functor->t_vH2v;
            *tools::common::profile::t_opt_dir_vH += *functor->t_vH;
            *tools::common::profile::t_opt_dir_vHv += *functor->t_vHv;
            break;
        }
    }
    // Copy and set the rest of the tensor metadata
    optimized_tensor.normalize();
    optimized_tensor.set_iter(summary.iterations.size());
    optimized_tensor.set_time(summary.total_time_in_seconds);
    optimized_tensor.set_energy(tools::finite::measure::energy(optimized_tensor.get_tensor(), tensors));
    optimized_tensor.set_variance(tools::finite::measure::energy_variance(optimized_tensor.get_tensor(), tensors));
    optimized_tensor.set_overlap(std::abs(current_vector.dot(optimized_tensor.get_vector())));

    tools::common::profile::t_opt_dir_bfgs->toc();
    reports::time_add_dir_entry();

    tools::log->debug("Finished LBFGS after {} seconds ({} iters). Exit status: {}. Message: {}", summary.total_time_in_seconds, summary.iterations.size(),
                      ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());

    reports::bfgs_add_entry("Direct", "opt", optimized_tensor);

    tools::log->trace("Returning theta from optimization mode {} space {}", optMode, optSpace);
    tools::common::profile::t_opt_dir->toc();
    return optimized_tensor;
}
