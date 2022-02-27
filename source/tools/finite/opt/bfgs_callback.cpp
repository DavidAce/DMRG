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
    if(not log) log = tools::Logger::setLogger("xDMRG");
    log->set_level(tools::log->level());
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
    if(log->level() >= spdlog::level::info) return ceres::SOLVER_CONTINUE;
    bool log_iter = summary.iteration - last_log_iter >= freq_log_iter and                  // log every freq_log_iter
                    summary.iteration >= init_log_iter;                                     // when at least init_log_iter have passed
    bool log_time = summary.cumulative_time_in_seconds - last_log_time >= freq_log_time and // log every last_log_time
                    summary.cumulative_time_in_seconds >= init_log_time;                    // when at least init_log_time have passed
    if(not log_iter or not log_time) return ceres::SOLVER_CONTINUE;
    last_log_time = summary.cumulative_time_in_seconds;
    last_log_iter = summary.iteration;
    /* clang-format off */
    log->debug(FMT_STRING("BFGS: "
               "it {:>5} f {:>8.5f} |Δf| {:>8.2e} "
               "|∇f| {:>8.2e} |∇f|ᵐᵃˣ {:>8.2e} |∇Ψ|ᵐᵃˣ {:>8.2e} "
               "|ΔΨ| {:8.2e} |Ψ|-1 {:8.2e} "
               "ops {:>4}/{:<4} t {:>8.2e} s op {:>8.2e} s Gop/s {:>5.1f} "
               "ls:[ss {:8.2e} fe {:>4} ge {:>4} it {:>4}] "
//               "| energy {:<18.15f} lg var {:<6.6f}"
               ),
               summary.iteration,
               summary.cost,
               summary.cost_change,
               summary.gradient_norm,
               summary.gradient_max_norm,
               functor.get_max_grad_norm(),
               summary.step_norm, // By bfgs
               functor.get_norm_offset(),
               functor.get_count() - last_count,
               functor.get_count(),
               summary.cumulative_time_in_seconds,
               summary.iteration_time_in_seconds,
               static_cast<double>(functor.get_ops()) / functor.t_H2n->get_last_interval()/1e9,
               std::abs(summary.step_size), // By line search
               summary.line_search_function_evaluations,
               summary.line_search_gradient_evaluations,
               summary.line_search_iterations
//             ,functor.get_energy_per_site(),
//               std::log10(functor.get_variance()
               );
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
