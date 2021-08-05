#include "../opt_meta.h"
#include "../opt_mps.h"
#include "ceres_direct_functor.h"
#include "opt-internal.h"
#include "report.h"
#include <algorithms/AlgorithmStatus.h>
#include <ceres/gradient_problem.h>
#include <config/settings.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/state/StateFinite.h>
#include <tensors/TensorsFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>

tools::finite::opt::opt_mps tools::finite::opt::internal::ceres_direct_optimization(const TensorsFinite &tensors, const AlgorithmStatus &status,
                                                                                    OptMeta &meta) {
    auto t_dir = tid::tic_scope("direct");
    std::vector<size_t> sites(tensors.active_sites.begin(), tensors.active_sites.end());
    opt_mps             initial_state("current state", tensors.state->get_multisite_mps(), sites,
                                      tools::finite::measure::energy(tensors) - tensors.model->get_energy_reduced(), // Eigval
                                      tensors.model->get_energy_reduced(),                                           // Energy reduced for full system
                                      tools::finite::measure::energy_variance(tensors),
                                      1.0, // Overlap
                                      tensors.get_length());
    t_dir.toc(); // The next call to ceres_direct_optimization opens the "direct" scope again.
    return ceres_direct_optimization(tensors, initial_state, status, meta);
}

tools::finite::opt::opt_mps tools::finite::opt::internal::ceres_direct_optimization(const TensorsFinite &tensors, const opt_mps &initial_mps,
                                                                                    const AlgorithmStatus &status, OptMeta &meta) {
    tools::log->trace("Optimizing in DIRECT mode");
    auto t_dir = tid::tic_scope("direct");
    initial_mps.validate_candidate();
    if constexpr(settings::debug)
        if(initial_mps.has_nan()) throw std::runtime_error("initial_mps has nan's");
    reports::bfgs_add_entry("Direct", "init", initial_mps);
    const auto &current_mps = tensors.state->get_multisite_mps();
    const auto  current_map = Eigen::Map<const Eigen::VectorXcd>(current_mps.data(), current_mps.size());

    auto    options = internal::ceres_default_options;
    auto    summary = ceres::GradientProblemSolver::Summary();
    opt_mps optimized_mps;
    optimized_mps.set_name(initial_mps.get_name());
    optimized_mps.set_sites(initial_mps.get_sites());
    optimized_mps.set_length(initial_mps.get_length());
    optimized_mps.set_energy_reduced(initial_mps.get_energy_reduced());

    auto t_lbfgs = tid::tic_scope("lbfgs");
    switch(meta.optType) {
        case OptType::CPLX: {
            // Copy the initial guess and operate directly on it
            optimized_mps.set_tensor(initial_mps.get_tensor());
            auto *functor = new ceres_direct_functor<std::complex<double>>(tensors, status);
            if(settings::precision::use_compressed_mpo_squared_otf)
                // Compress the virtual bond between MPO² and the environments
                functor->compress();

            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS direct cplx");
            ceres::Solve(options, problem, optimized_mps.get_vector_cplx_as_2xreal().data(), &summary);
            // Copy the results from the functor
            optimized_mps.set_counter(functor->get_count());
            optimized_mps.set_delta_f(functor->get_delta_f());
            optimized_mps.set_grad_norm(functor->get_grad_max_norm());
            tid::get("vH2") += *functor->t_H2n;
            tid::get("vH2v") += *functor->t_nH2n;
            tid::get("vH") += *functor->t_Hn;
            tid::get("vHv") += *functor->t_nHn;
            tid::get("step") += *functor->t_step;
            break;
        }
        case OptType::REAL: {
            // Here we make a temporary
            auto  initial_state_real = initial_mps.get_vector_cplx_as_1xreal();
            auto *functor            = new ceres_direct_functor<double>(tensors, status);
            if(settings::precision::use_compressed_mpo_squared_otf)
                // Compress the virtual bond between MPO² and the environments
                functor->compress();

            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS direct real");
            if constexpr(settings::debug)
                if(initial_state_real.hasNaN()) throw std::runtime_error("initial_state_real has nan's");
            ceres::Solve(options, problem, initial_state_real.data(), &summary);
            // Copy the results from the functor
            optimized_mps.set_counter(functor->get_count());
            optimized_mps.set_delta_f(functor->get_delta_f());
            optimized_mps.set_grad_norm(functor->get_grad_max_norm());
            optimized_mps.set_tensor_real(initial_state_real.data(), initial_mps.get_tensor().dimensions());
            tid::get("vH2") += *functor->t_H2n;
            tid::get("vH2v") += *functor->t_nH2n;
            tid::get("vH") += *functor->t_Hn;
            tid::get("vHv") += *functor->t_nHn;
            tid::get("step") += *functor->t_step;
            break;
        }
    }
    reports::time_add_opt_entry();
    t_lbfgs.toc();
    // Copy and set the rest of the tensor metadata
    optimized_mps.normalize();
    optimized_mps.set_energy(tools::finite::measure::energy(optimized_mps.get_tensor(), tensors));
    optimized_mps.set_variance(tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors));
    optimized_mps.set_overlap(std::abs(current_map.dot(optimized_mps.get_vector())));
    optimized_mps.set_iter(summary.iterations.size());
    optimized_mps.set_time(summary.total_time_in_seconds);

    int    hrs = static_cast<int>(summary.total_time_in_seconds / 3600);
    int    min = static_cast<int>(std::fmod(summary.total_time_in_seconds, 3600) / 60);
    double sec = std::fmod(std::fmod(summary.total_time_in_seconds, 3600), 60);
    tools::log->debug("Finished LBFGS in {:0<2}:{:0<2}:{:0<.1f} seconds and {} iters. Exit status: {}. Message: {}", hrs, min, sec, summary.iterations.size(),
                      ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
    reports::bfgs_add_entry("Direct", "opt", optimized_mps);

    return optimized_mps;
}
