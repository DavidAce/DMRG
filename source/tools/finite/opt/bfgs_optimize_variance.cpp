#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/state/StateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt/bfgs_callback.h"
#include "tools/finite/opt/bfgs_variance_functor.h"
#include "tools/finite/opt/opt-internal.h"
#include "tools/finite/opt/report.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"
#include <ceres/gradient_problem.h>

tools::finite::opt::opt_mps tools::finite::opt::internal::bfgs_optimize_variance(const TensorsFinite &tensors, const opt_mps &initial_mps,
                                                                                 const AlgorithmStatus &status, OptMeta &meta) {
    tools::log->trace("Optimizing variance");
    auto t_var = tid::tic_scope("variance");
    initial_mps.validate_initial_mps();
    if constexpr(settings::debug)
        if(initial_mps.has_nan()) throw std::runtime_error("initial_mps has nan's");
    reports::bfgs_add_entry("variance", "init", initial_mps);
    const auto &current_mps = tensors.state->get_multisite_mps();
    const auto  current_map = Eigen::Map<const Eigen::VectorXcd>(current_mps.data(), current_mps.size());

    auto    options = internal::bfgs_default_options;
    auto    summary = ceres::GradientProblemSolver::Summary();
    opt_mps optimized_mps;
    optimized_mps.set_sites(initial_mps.get_sites());
    optimized_mps.set_length(initial_mps.get_length());
    optimized_mps.set_energy_shift(initial_mps.get_energy_shift());

    auto t_bfgs = tid::tic_scope("bfgs");
    switch(meta.optType) {
        case OptType::CPLX: {
            auto  initial_guess = initial_mps.get_initial_state_with_lagrange_multiplier<OptType::CPLX>();
            auto *functor       = new bfgs_variance_functor<cplx, LagrangeNorm::ON>(tensors, status);
            //            auto *param         = new NormParametrization(*functor);
            if(settings::precision::use_compressed_mpo_squared_otf)
                // Compress the virtual bond between MPO² and the environments
                functor->compress();

            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running BFGS cplx");

            ceres::Solve(options, problem, initial_guess.data(), &summary);
            // Copy the results from the functor
            optimized_mps.set_name("bfgs-cplx");
            optimized_mps.set_mv(functor->get_count());
            optimized_mps.set_delta_f(functor->get_delta_f());
            optimized_mps.set_grad_max(functor->get_max_grad_norm());
            optimized_mps.set_tensor_cplx(initial_guess.data(), initial_mps.get_tensor().dimensions());
            break;
        }
        case OptType::REAL: {
            // Here we make a temporary
            auto  initial_guess = initial_mps.get_initial_state_with_lagrange_multiplier<OptType::REAL>();
            auto *functor       = new bfgs_variance_functor<real, LagrangeNorm::ON>(tensors, status);
            //            auto *param         = new NormParametrization(*functor);
            if(settings::precision::use_compressed_mpo_squared_otf)
                // Compress the virtual bond between MPO² and the environments
                functor->compress();
            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running BFGS real");
            ceres::Solve(options, problem, initial_guess.data(), &summary);
            // Copy the results from the functor
            optimized_mps.set_name("bfgs-real");
            optimized_mps.set_mv(functor->get_count());
            optimized_mps.set_delta_f(functor->get_delta_f());
            optimized_mps.set_grad_max(functor->get_max_grad_norm());
            optimized_mps.set_tensor_real(initial_guess.data(), initial_mps.get_tensor().dimensions());
            break;
        }
    }
    reports::time_add_entry();
    t_bfgs.toc();
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
    tools::log->debug("Finished BFGS in {:0<2}:{:0<2}:{:0<.1f} seconds and {} iters. Exit status: {}. Message: {}", hrs, min, sec, summary.iterations.size(),
                      ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
    reports::bfgs_add_entry("variance", "opt", optimized_mps);
    return optimized_mps;
}
