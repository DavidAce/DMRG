#include "../opt_meta.h"
#include "../opt_mps.h"
#include "lbfgs_variance_functor.h"
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

tools::finite::opt::opt_mps tools::finite::opt::internal::lbfgs_optimize_variance(const TensorsFinite &tensors, const AlgorithmStatus &status, OptMeta &meta) {
    auto                t_var = tid::tic_scope("variance");
    std::vector<size_t> sites(tensors.active_sites.begin(), tensors.active_sites.end());
    opt_mps             initial_state("current state", tensors.state->get_multisite_mps(), sites,
                                      tools::finite::measure::energy(tensors) - tensors.model->get_energy_shift(), // Eigval
                                      tensors.model->get_energy_shift(),                                           // Shifted energy for full system
                                      tools::finite::measure::energy_variance(tensors),
                                      1.0, // Overlap
                                      tensors.get_length());
    t_var.toc(); // The next call to lbfgs_optimize_variance opens the "direct" scope again.
    return lbfgs_optimize_variance(tensors, initial_state, status, meta);
}

tools::finite::opt::opt_mps tools::finite::opt::internal::lbfgs_optimize_variance(const TensorsFinite &tensors, const opt_mps &initial_mps,
                                                                                  const AlgorithmStatus &status, OptMeta &meta) {
    tools::log->trace("Optimizing variance");
    auto t_var = tid::tic_scope("variance");
    initial_mps.validate_basis_vector();
    if constexpr(settings::debug)
        if(initial_mps.has_nan()) throw std::runtime_error("initial_mps has nan's");
    reports::bfgs_add_entry("Direct", "init", initial_mps);
    const auto &current_mps = tensors.state->get_multisite_mps();
    const auto  current_map = Eigen::Map<const Eigen::VectorXcd>(current_mps.data(), current_mps.size());

    auto    options = internal::lbfgs_default_options;
    auto    summary = ceres::GradientProblemSolver::Summary();
    opt_mps optimized_mps;
    optimized_mps.set_name(initial_mps.get_name());
    optimized_mps.set_sites(initial_mps.get_sites());
    optimized_mps.set_length(initial_mps.get_length());
    optimized_mps.set_energy_shift(initial_mps.get_energy_shift());

    auto t_lbfgs = tid::tic_scope("lbfgs");
    switch(meta.optType) {
        case OptType::CPLX: {
            auto  initial_guess = initial_mps.get_initial_state_with_lagrange_multiplier<OptType::CPLX>();
            auto *functor       = new lbfgs_variance_functor<std::complex<double>, LagrangeNorm::ON>(tensors, status);
            if(settings::precision::use_compressed_mpo_squared_otf)
                // Compress the virtual bond between MPO² and the environments
                functor->compress();

            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS direct cplx");

            ceres::Solve(options, problem, optimized_mps.get_vector_cplx_as_2xreal().data(), &summary);
            // Copy the results from the functor
            optimized_mps.set_mv(functor->get_count());
            optimized_mps.set_delta_f(functor->get_delta_f());
            optimized_mps.set_max_grad(functor->get_max_grad_norm());
            optimized_mps.set_tensor_cplx(initial_guess.data(), initial_mps.get_tensor().dimensions());
            tid::get("vH2") += *functor->t_H2n;
            tid::get("vH2v") += *functor->t_nH2n;
            tid::get("vH") += *functor->t_Hn;
            tid::get("vHv") += *functor->t_nHn;
            tid::get("step") += *functor->t_step;
            break;
        }
        case OptType::REAL: {
            // Here we make a temporary
            auto  initial_guess = initial_mps.get_initial_state_with_lagrange_multiplier<OptType::REAL>();
            auto *functor       = new lbfgs_variance_functor<double, LagrangeNorm::ON>(tensors, status);
            if(settings::precision::use_compressed_mpo_squared_otf)
                // Compress the virtual bond between MPO² and the environments
                functor->compress();
            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS direct real");
            ceres::Solve(options, problem, initial_guess.data(), &summary);
            // Copy the results from the functor
            optimized_mps.set_mv(functor->get_count());
            optimized_mps.set_delta_f(functor->get_delta_f());
            optimized_mps.set_max_grad(functor->get_max_grad_norm());
            optimized_mps.set_tensor_real(initial_guess.data(), initial_mps.get_tensor().dimensions());
            tid::get("vH2") += *functor->t_H2n;
            tid::get("vH2v") += *functor->t_nH2n;
            tid::get("vH") += *functor->t_Hn;
            tid::get("vHv") += *functor->t_nHn;
            tid::get("step") += *functor->t_step;
            break;
        }
    }
    reports::time_add_entry();
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
    reports::bfgs_add_entry("variance", "opt", optimized_mps);
    return optimized_mps;
}
