#include "bfgs_callback.h"
#include "config/debug.h"
#include "debug/exceptions.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/finite/opt/bfgs_base_functor.h"
#include "tools/finite/opt/bfgs_subspace_functor.h"
#include "tools/finite/opt/bfgs_variance_functor.h"

namespace settings {
    constexpr static bool debug_callback = false;
}

template<typename FunctorType>
tools::finite::opt::internal::CustomLogCallback<FunctorType>::CustomLogCallback(const FunctorType &functor_) : functor(functor_) {
    if(not log) log = tools::Logger::setLogger("LBFGS");
    log->set_level(tools::log->level());
    if(log->level() <= spdlog::level::debug) {
        freq_log_iter = 1000;
        freq_log_time = 10;
        init_log_time = 0;
    }
    if constexpr(settings::debug) {
        freq_log_iter = 1;
        freq_log_time = 0;
        init_log_time = 0;
    }

    if constexpr(settings::debug_callback) {
        log->set_level(spdlog::level::debug);
        freq_log_iter = 1;
        freq_log_time = 0;
        init_log_time = 0;
    }
}

template<typename FunctorType>
ceres::CallbackReturnType tools::finite::opt::internal::CustomLogCallback<FunctorType>::operator()(const ceres::IterationSummary &summary) {
    functor.set_delta_f(std::abs(summary.cost_change / summary.cost));
    if(not log) return ceres::SOLVER_CONTINUE;
    if(log->level() > spdlog::level::info) return ceres::SOLVER_CONTINUE;
    double time_since_last_log = summary.cumulative_time_in_seconds - last_log_time;
    double iter_since_last_log = summary.iteration - last_log_iter;

    bool log_iter = iter_since_last_log >= freq_log_iter and             // log every freq_log_iter
                    summary.iteration >= init_log_iter;                  // when at least init_log_iter have passed
    bool log_time = time_since_last_log >= freq_log_time and             // log every last_log_time
                    summary.cumulative_time_in_seconds >= init_log_time; // when at least init_log_time have passed
    if(not log_iter or not log_time) return ceres::SOLVER_CONTINUE;

    /* clang-format off */
    log->info(FMT_STRING(
               "it {:>5} dims {}={} σ² {:>11.5e} f {:>8.5f} |Δf| {:>8.2e} "
               "|∇f| {:>8.2e} |∇f|ᵐᵃˣ {:>8.2e} |∇Ψ|ᵐᵃˣ {:>8.2e} "
               "|ΔΨ| {:8.2e} |Ψ|-1 {:8.2e} rnorm {:8.2e} "
               "ops {:>4}/{:<4} t {:>8.2e} it/s {:>8.2e}  Gop/s {:>5.1f} "
               "ls:[ss {:8.2e} fe {:>4} ge {:>4} it {:>4}] "
//               "| energy {:<18.15f} lg var {:<6.6f}"
               ),
               summary.iteration,
               functor.get_dims(),
               functor.get_size(),
               functor.get_variance(),
               summary.cost,
               summary.cost_change,
               summary.gradient_norm,
               summary.gradient_max_norm,
               functor.get_max_grad_norm(),
               summary.step_norm, // By bfgs
               functor.get_norm_offset(),
               functor.get_resnorm(),
               functor.get_count() - last_count,
               functor.get_count(),
               summary.cumulative_time_in_seconds,
               iter_since_last_log/time_since_last_log,
               static_cast<double>(functor.get_ops()) / functor.t_H2n->get_last_interval()/1e9,
               std::abs(summary.step_size), // By line search
               summary.line_search_function_evaluations,
               summary.line_search_gradient_evaluations,
               summary.line_search_iterations
               );
    last_log_time = summary.cumulative_time_in_seconds;
    last_log_iter = summary.iteration;
    last_count = functor.get_count();
    /* clang-format on */
    return ceres::SOLVER_CONTINUE;
}

namespace tools::finite::opt::internal {
    template class CustomLogCallback<bfgs_subspace_functor<real>>;
    template class CustomLogCallback<bfgs_variance_functor<real, LagrangeNorm::ON>>;
    template class CustomLogCallback<bfgs_variance_functor<real, LagrangeNorm::OFF>>;
    template class CustomLogCallback<bfgs_subspace_functor<cplx>>;
    template class CustomLogCallback<bfgs_variance_functor<cplx, LagrangeNorm::ON>>;
    template class CustomLogCallback<bfgs_variance_functor<cplx, LagrangeNorm::OFF>>;
}
