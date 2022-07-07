#include "algorithms/AlgorithmStatus.h"
#include "config/debug.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "math/num.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt/bfgs_base_functor.h"
#include "tools/finite/opt/bfgs_subspace_functor.h"
#include "tools/finite/opt/bfgs_variance_functor.h"
#include "tools/finite/opt/opt-internal.h"
#include "tools/finite/opt/report.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"
#include <glog/logging.h>
#include <string>

tools::finite::opt::opt_mps tools::finite::opt::get_opt_initial_mps(const TensorsFinite &tensors) {
    auto    t_init = tid::tic_scope("initial_mps");
    opt_mps initial_mps("current mps", tensors.get_multisite_mps(), tensors.active_sites,
                        tools::finite::measure::energy_minus_energy_shift(tensors), // Eigval
                        tools::finite::measure::energy_shift(tensors),              // Shifted energy for full system
                        tools::finite::measure::energy_variance(tensors),
                        1.0, // Overlap to initial state (self) is obviously 1
                        tensors.get_length());
    return initial_mps;
}

tools::finite::opt::opt_mps tools::finite::opt::find_excited_state(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status,
                                                                   OptMeta &meta) {
    auto t_opt = tid::tic_scope("opt");
    tools::log->trace("Starting optimization: mode [{}] | solver [{}] | type [{}] | position [{}] | sites {} | shape {} = {}", enum2sv(meta.optMode),
                      enum2sv(meta.optSolver), enum2sv(meta.optType), status.position, tensors.active_sites, tensors.active_problem_dims(),
                      tensors.active_problem_size());

    using namespace opt::internal;
    static bool googleLogginghasInitialized = false;
    if(not googleLogginghasInitialized) {
        googleLogginghasInitialized = true;
#if defined(CERCES_INTERNAL_MINIGLOG_GLOG_LOGGING_H_)
        google::InitGoogleLogging(const_cast<char *>(tools::log->name().c_str()));
#else
        google::InitGoogleLogging(tools::log->name().c_str());
        google::SetStderrLogging(3);
//#pragma message "Using verbose glog"
//        google::SetStderrLogging(0);
#endif
    }

    /* clang-format off */

    /* Alglib has a good read comparing NCG and LBFGS https://www.alglib.net/optimization/lbfgsandcg.php
     * Some noteworthy points:
     *      "L-BFGS algorithm is quite easy to tune.
     *       The only parameter to tune is number of correction pairs M - number of function/gradient pairs
     *       used to build quadratic model (i.e. rank). This parameter is passed to the function which is used
     *       to create optimizer object. Increase of M decreases number of function evaluations and iterations
     *       needed to converge, but increases computational overhead associated with iteration. On well-conditioned
     *       problems M can be as small as 3-10. If function changes rapidly in some directions and slowly in
     *       other ones, then you can try increasing M in order to increase convergence."
     *
     *  Testing here shows that lbfgs is much better than NCG.
     *  Testing also shows that increasing lbfgs rank can give better results at some cost in some cases,
     *  in particular when a problem is ill-conditioned.
     *
     */

    /*
     * When a problem is ill-conditioned, for instance when the difference in eigenvalues is small,
     * we need to increase the lbfgs rank.
     * Read more here:
     *   kth.diva-portal.org/smash/get/diva2:1438308/FULLTEXT01.pdf (page 22)
     *   https://www.deeplearningbook.org/contents/optimization.html (page 307)
     *
     */


    bfgs_default_options.line_search_type                           = ceres::LineSearchType::WOLFE;
    bfgs_default_options.line_search_interpolation_type             = ceres::LineSearchInterpolationType::CUBIC;
    bfgs_default_options.line_search_direction_type                 = ceres::LineSearchDirectionType::LBFGS;
    bfgs_default_options.max_num_iterations                         = 2000;
    bfgs_default_options.max_lbfgs_rank                             = 16; // Tested: around 8-32 seems to be a good compromise,but larger is more precise sometimes. Overhead goes from 1.2x to 2x computation time at in 8 -> 64

    // Approximate eigenvalue bfgs scaling performs badly when the problem
    // is ill-conditioned (sensitive to some parameters). Empirically,
    // we observe that the gradient is still large when bfgs has finished,
    // and often the result is a tiny bit worse than what we started with.
    // When this happens, it's not worth trying to get BFGS to converge,
    // try an eigensolver on the energy-shifted operator H² instead.
    bfgs_default_options.use_approximate_eigenvalue_bfgs_scaling    = true;  // Tested: True makes a huge difference, takes longer steps at each iteration and generally converges faster/to better variance
    bfgs_default_options.min_line_search_step_size                  = std::numeric_limits<double>::epsilon();
    bfgs_default_options.max_line_search_step_contraction           = 1e-3; // 1e-3
    bfgs_default_options.min_line_search_step_contraction           = 0.6; // 0.6
    bfgs_default_options.max_line_search_step_expansion             = 10; // 10
    bfgs_default_options.max_num_line_search_step_size_iterations   = 20;
    bfgs_default_options.max_num_line_search_direction_restarts     = 5; //5
    bfgs_default_options.line_search_sufficient_function_decrease   = 1e-4; //1e-4; Tested, doesn't seem to matter between [1e-1 to 1e-4]. Default is fine: 1e-4
    bfgs_default_options.line_search_sufficient_curvature_decrease  = 1-1e-3;//0.9 // This one should be above 0.5. Below, it makes retries at every step and starts taking twice as long for no added benefit. Tested 0.9 to be sweetspot
    bfgs_default_options.max_solver_time_in_seconds                 = 60*60;//60*2;
    bfgs_default_options.function_tolerance                         = 1e-6; // Tested, 1e-6 seems to be a sweetspot
    bfgs_default_options.gradient_tolerance                         = 1e-4; // This is the max gradient on f = log Var H
    bfgs_default_options.parameter_tolerance                        = 1e-14;
    bfgs_default_options.minimizer_progress_to_stdout               = false; //tools::log->level() <= spdlog::level::trace;
    bfgs_default_options.update_state_every_iteration               = false;
    bfgs_default_options.logging_type                               = ceres::LoggingType::PER_MINIMIZER_ITERATION;

    /* clang-format off */
    // Apply overrides if there are any
    if(meta.bfgs_max_iter) bfgs_default_options.max_num_iterations  = meta.bfgs_max_iter.value();
    if(meta.bfgs_max_rank) bfgs_default_options.max_lbfgs_rank      = meta.bfgs_max_rank.value();
    if(meta.bfgs_func_tol) bfgs_default_options.function_tolerance  = meta.bfgs_func_tol.value();
    if(meta.bfgs_grad_tol) bfgs_default_options.gradient_tolerance  = meta.bfgs_grad_tol.value();
    if(meta.bfgs_eigenvalue_scaling) bfgs_default_options.use_approximate_eigenvalue_bfgs_scaling = meta.bfgs_eigenvalue_scaling.value();
    /* clang-format on */

    //    Progress log definitions:
    //    f is the value of the objective function.
    //    d is the change in the value of the objective function if the step computed in this iteration is accepted.
    //    g is the inf norm of the gradient (i.e. the largest element of |grad f| )
    //    h is the change in the parameter vector.
    //    s is the optimal step length computed by the line search.
    //    it is the time take by the current iteration.
    //    tt is the total time taken by the minimizer.

    if(initial_mps.get_sites() != tensors.active_sites)
        throw except::runtime_error("mismatch in active sites: initial_mps {} | active {}", initial_mps.get_sites(), tensors.active_sites);
    /* clang-format off */
    opt_mps result;
    if     (meta.optMode == OptMode::OVERLAP  and meta.optSolver == OptSolver::EIGS)  result = internal::eigs_optimize_overlap(tensors, initial_mps, status, meta);
    else if(meta.optMode == OptMode::ENERGY   and meta.optSolver == OptSolver::EIGS)  result = internal::eigs_optimize_energy(tensors, initial_mps, status, meta); // TODO: Implement energy mode
    else if(meta.optMode == OptMode::SUBSPACE and meta.optSolver == OptSolver::EIGS)  result = internal::eigs_optimize_subspace(tensors, initial_mps, status, meta);
    else if(meta.optMode == OptMode::SIMPS    and meta.optSolver == OptSolver::EIGS)  result = internal::eigs_optimize_variance(tensors, initial_mps, status, meta); // TODO: Implement simps mode
    else if(meta.optMode == OptMode::VARIANCE and meta.optSolver == OptSolver::EIGS)  result = internal::eigs_optimize_variance(tensors, initial_mps, status, meta);
    else if(meta.optMode == OptMode::VARIANCE and meta.optSolver == OptSolver::BFGS)  result = internal::bfgs_optimize_variance(tensors, initial_mps, status, meta);
    else
        throw except::logic_error("Incompatible: OptMode [{}] and OptSolver [{}]", enum2sv(meta.optMode), enum2sv(meta.optSolver));
    /* clang-format on */

    tools::finite::opt::reports::print_eigs_report();
    tools::finite::opt::reports::print_bfgs_report();
    tools::finite::opt::reports::print_time_report();

    // Finish up
    result.set_optsolver(meta.optSolver);
    result.set_optmode(meta.optMode);
    result.set_optexit(meta.optExit);
    result.validate_result();
    return result;
}

tools::finite::opt::opt_mps tools::finite::opt::find_ground_state(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status,
                                                                  OptMeta &meta) {
    return internal::eigs_optimize_energy(tensors, initial_mps, status, meta);
}

double tools::finite::opt::internal::windowed_func_abs(double x, double window) {
    if(std::abs(x) >= window)
        return std::abs(x) - window;
    else
        return 0;
}
double tools::finite::opt::internal::windowed_grad_abs(double x, double window) {
    if(std::abs(x) >= window)
        return num::sign(x);
    else
        return 0.0;
}

double tools::finite::opt::internal::windowed_func_pow(double x, double window) {
    if(std::abs(x) >= window)
        return x * x - window * window;
    else
        return 0.0;
}
double tools::finite::opt::internal::windowed_grad_pow(double x, double window) {
    if(std::abs(x) >= window)
        return 2.0 * x;
    else
        return 0.0;
}

std::pair<double, double> tools::finite::opt::internal::windowed_func_grad(double x, double window) {
    double func = 0;
    double grad = 0;
    if(std::abs(x) >= window) {
        if(x > 0) {
            func = (x - window) * (x - window);
            grad = 2 * (x - window);
        } else {
            func = (x + window) * (x + window);
            grad = 2 * (x + window);
        }
    }

    return std::make_pair(func, grad);
}

long tools::finite::opt::internal::get_ops(long d, long chiL, long chiR, long m) {
    if(chiR > chiL) return get_ops(d, chiR, chiL, m);
    if(d > m) {
        // d first
        long step1 = chiL * chiL * chiR * m * m * d;
        long step2 = chiL * chiR * d * m * m * m * (d * m + 1);
        long step3 = chiL * chiR * d * m * (chiR * m * m * m + m * m + 1);
        return step1 + step2 + step3;
    } else {
        // m first
        long step1 = chiL * chiL * chiR * m * d;
        long step2 = chiL * chiR * d * m * m * (d * m + 1);
        long step3 = chiL * chiR * chiR * d * m * (m + 1);
        return step1 + step2 + step3;
    }
}

template<typename FunctorType>
tools::finite::opt::internal::NormParametrization<FunctorType>::NormParametrization(const FunctorType &functor_) : functor(functor_) {}

template<typename FunctorType>
bool tools::finite::opt::internal::NormParametrization<FunctorType>::Plus(const double *v_double_double, const double *d_double_double,
                                                                          double *v_plus_d) const {
    auto t_plus      = tid::tic_scope("plus");
    using VectorType = typename FunctorType::VectorType;
    using ScalarType = typename VectorType::Scalar;

    auto v = Eigen::Map<const VectorType>(reinterpret_cast<const ScalarType *>(v_double_double), functor.size);
    auto d = Eigen::Map<const VectorType>(reinterpret_cast<const ScalarType *>(d_double_double), functor.size);
    auto p = Eigen::Map<VectorType>(reinterpret_cast<ScalarType *>(v_plus_d), functor.size);
    p      = (v + d).normalized();
    return true;
}

template<typename FunctorType>
bool tools::finite::opt::internal::NormParametrization<FunctorType>::PlusJacobian(const double *v_double_double, double *jac_double_double) const {
    auto t_jac       = tid::tic_scope("jac");
    using VectorType = typename FunctorType::VectorType;
    using ScalarType = typename VectorType::Scalar;
    using MatrixType = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    auto   v           = Eigen::Map<const VectorType>(reinterpret_cast<const ScalarType *>(v_double_double), functor.size);
    auto   jac         = Eigen::Map<MatrixType>(reinterpret_cast<ScalarType *>(jac_double_double), functor.size, functor.size);
    double n           = v.norm();
    double one_over_n  = 1.0 / n;
    double one_over_nn = one_over_n * one_over_n;
    //    jac      = (MatrixType::Identity(v.size(), v.size()) - v * v.transpose() / (n * n)) / n;
#pragma omp parallel for schedule(dynamic)
    for(long i = 0; i < jac.rows(); i++) {
        for(long j = i; j < jac.cols(); j++) {
            auto id   = static_cast<double>(i == j);
            auto res  = (id - v[i] * v[j] * one_over_nn) * one_over_n;
            jac(i, j) = res;
            jac(j, i) = res;
        }
    }
    return true;
}

template<typename FunctorType>
int tools::finite::opt::internal::NormParametrization<FunctorType>::AmbientSize() const {
    return functor.num_parameters;
}
template<typename FunctorType>
int tools::finite::opt::internal::NormParametrization<FunctorType>::TangentSize() const {
    return functor.num_parameters;
}

namespace tools::finite::opt::internal {
    template class NormParametrization<bfgs_variance_functor<real, LagrangeNorm::OFF>>;
    template class NormParametrization<bfgs_variance_functor<cplx, LagrangeNorm::OFF>>;
    template class NormParametrization<bfgs_subspace_functor<real>>;
    template class NormParametrization<bfgs_subspace_functor<cplx>>;
}

namespace tools::finite::opt::internal {
    bool comparator::energy(const opt_mps &lhs, const opt_mps &rhs) {
        // The eigenvalue solver on H gives results sorted in energy
        if(lhs.get_eigs_idx() != rhs.get_eigs_idx()) return lhs.get_eigs_idx() < rhs.get_eigs_idx();
        return lhs.get_energy() < rhs.get_energy();
    }

    bool comparator::energy_distance(const opt_mps &lhs, const opt_mps &rhs, double target) {
        if(std::isnan(target)) throw except::logic_error("Energy target for comparison is NAN");
        auto diff_energy_lhs = std::abs(lhs.get_energy() - target);
        auto diff_energy_rhs = std::abs(rhs.get_energy() - target);
        return diff_energy_lhs < diff_energy_rhs;
    }

    bool comparator::variance(const opt_mps &lhs, const opt_mps &rhs) {
        // The eigenvalue solver on (H-E)² gives results sorted in variance
        if(lhs.get_eigs_idx() != rhs.get_eigs_idx()) return lhs.get_eigs_idx() < rhs.get_eigs_idx();
        return lhs.get_variance() < rhs.get_variance();
    }

    bool comparator::gradient(const opt_mps &lhs, const opt_mps &rhs) {
        if(lhs.get_eigs_idx() != rhs.get_eigs_idx()) return lhs.get_eigs_idx() < rhs.get_eigs_idx();
        return lhs.get_grad_max() < rhs.get_grad_max();
    }
    bool comparator::eigval(const opt_mps &lhs, const opt_mps &rhs) {
        if(lhs.get_eigs_idx() != rhs.get_eigs_idx()) return lhs.get_eigs_idx() < rhs.get_eigs_idx();
        return lhs.get_eigs_eigval() < rhs.get_eigs_eigval();
    }

    bool comparator::overlap(const opt_mps &lhs, const opt_mps &rhs) { return lhs.get_overlap() > rhs.get_overlap(); }

    bool comparator::eigval_and_overlap(const opt_mps &lhs, const opt_mps &rhs) {
        double ratio = std::max(lhs.get_eigs_eigval(), rhs.get_eigs_eigval()) / std::min(lhs.get_eigs_eigval(), rhs.get_eigs_eigval());
        if(ratio < 10 and lhs.get_overlap() >= std::sqrt(0.5)) return comparator::overlap(lhs, rhs);
        return comparator::eigval(lhs, rhs);
    }

    Comparator::Comparator(const OptMeta &meta_, double target_energy_) : meta(&meta_), target_energy(target_energy_) {}
    bool Comparator::operator()(const opt_mps &lhs, const opt_mps &rhs) {
        if(not meta) throw except::logic_error("No opt_meta given to comparator");
        // The variance case first
        if(meta->optMode == OptMode::VARIANCE) {
            return comparator::eigval_and_overlap(lhs, rhs);
        } else {
            auto diff = std::abs(lhs.get_eigval() - rhs.get_eigval());
            if(diff < settings::precision::eigs_tolerance) return lhs.get_overlap() > rhs.get_overlap();
            switch(meta->optRitz) {
                case OptRitz::SR: return comparator::energy(lhs, rhs);
                case OptRitz::LR: return comparator::energy(rhs, lhs);
                case OptRitz::SM: return comparator::energy_distance(lhs, rhs, target_energy);
            }
        }
        return comparator::eigval(lhs, rhs);
    }
}

namespace tools::finite::opt::internal {}