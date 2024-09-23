
#include "../opt_meta.h"
#include "../opt_mps.h"
#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "general/iter.h"
#include "general/sfinae.h"
#include "math/eig.h"
#include "math/eig/matvec/matvec_mpo.h"
#include "math/eig/matvec/matvec_mpos.h"
#include "math/eig/matvec/matvec_zero.h"
#include "tensors/edges/EdgesFinite.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/state/StateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt/opt-internal.h"
#include "tools/finite/opt/report.h"
#include <math/num.h>
#include <primme/primme.h>

namespace tools::finite::opt {
    namespace folded_spectrum {
        template<typename Scalar>
        void preconditioner_jacobi(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr) {
            if(x == nullptr) return;
            if(y == nullptr) return;
            if(primme == nullptr) return;
            const auto H_ptr      = static_cast<MatVecMPOS<Scalar> *>(primme->matrix);
            H_ptr->preconditioner = eig::Preconditioner::JACOBI;
            H_ptr->MultPc(x, ldx, y, ldy, blockSize, primme, ierr);
        }

        template<typename Scalar>
        void convTestFun([[maybe_unused]] double *eval, [[maybe_unused]] void *evec, double *rNorm, int *isconv, struct primme_params *primme, int *ierr) {
            if(rNorm == nullptr) return;
            if(primme == nullptr) return;

            double problemNorm;
            if(not primme->massMatrixMatvec) {
                problemNorm = primme->aNorm > 0.0 ? primme->aNorm : primme->stats.estimateLargestSVal;
            } else {
                problemNorm = primme->aNorm > 0.0 && primme->invBNorm > 0.0 ? primme->aNorm * primme->invBNorm : primme->stats.estimateLargestSVal;
            }
            double problemTol = problemNorm * primme->eps;
            double diff_rnorm = 0;
            double diff_eval  = 0;
            if(primme->monitor != nullptr and *rNorm != 0) {
                auto &solver = *static_cast<eig::solver *>(primme->monitor);
                auto &config = solver.config;
                auto &result = solver.result;
                result.meta.recent_evals.push_back(*eval);
                result.meta.recent_rnorms.push_back(*rNorm);
                while(result.meta.recent_evals.size() > 500) result.meta.recent_evals.pop_front();
                while(result.meta.recent_rnorms.size() > 500) result.meta.recent_rnorms.pop_front();
                double diff_evals_sum  = -1.0;
                double diff_rnorms_sum = -1.0;
                if(result.meta.recent_evals.size() >= 250) {
                    auto diff      = num::diff(result.meta.recent_evals);
                    diff_evals_sum = num::sum(diff, 1);
                    // tools::log->info("recent evals: {::.3e}", result.meta.recent_evals);
                    // tools::log->info("diff   evals: {::.3e}", diff);
                }
                if(result.meta.recent_rnorms.size() >= 250) {
                    auto diff       = num::diff(result.meta.recent_rnorms);
                    diff_rnorms_sum = num::sum(diff, 1);
                    // tools::log->info("recent rnorm: {::.3e}", result.meta.recent_rnorms);
                    // tools::log->info("diff   rnorm: {::.3e}", diff);
                }

                bool evals_saturated   = diff_evals_sum > -primme->eps;
                bool rnorms_saturated  = diff_rnorms_sum > -primme->eps;
                result.meta.last_eval  = *eval;
                result.meta.last_rnorm = *rNorm;

                *isconv = *rNorm < problemTol and (evals_saturated or rnorms_saturated) or *rNorm < primme->eps;
                if(*isconv == 1) {
                    tools::log->info("eval {:38.32f} ({:+8.3e}) | rnorm {:38.32f} ({:+8.3e}) | problemNorm {:.3e}", *eval, diff_evals_sum, *rNorm,
                                     diff_rnorms_sum, problemNorm);
                }
            } else {
                *isconv = 0;
            }

            *ierr = 0;
        }

        template<typename Scalar>
        struct opt_mps_init_t {
            Eigen::Tensor<Scalar, 3> mps = {};
            long                     idx = 0;
        };
        template<typename Scalar>
        struct opt_bond_init_t {
            Eigen::Tensor<Scalar, 2> bond = {};
            long                     idx  = 0;
        };
        template<typename Scalar>
        Eigen::Tensor<Scalar, 3> get_initial_guess(const opt_mps &initial_mps, const std::vector<opt_mps> &results) {
            if(results.empty()) {
                if constexpr(std::is_same_v<Scalar, real>)
                    return initial_mps.get_tensor().real();
                else
                    return initial_mps.get_tensor();
            } else {
                // Return whichever of initial_mps or results that has the lowest variance
                auto it = std::min_element(results.begin(), results.end(), internal::comparator::variance);
                if(it == results.end()) return get_initial_guess<Scalar>(initial_mps, {});

                if(it->get_variance() < initial_mps.get_variance()) {
                    tools::log->debug("Previous result is a good initial guess: {} | var {:8.2e}", it->get_name(), it->get_variance());
                    return get_initial_guess<Scalar>(*it, {});
                } else
                    return get_initial_guess<Scalar>(initial_mps, {});
            }
        }

        template<typename Scalar>
        std::vector<opt_mps_init_t<Scalar>> get_initial_guess_mps(const opt_mps &initial_mps, const std::vector<opt_mps> &results, long nev) {
            std::vector<opt_mps_init_t<Scalar>> init;
            if(results.empty()) {
                if constexpr(std::is_same_v<Scalar, real>)
                    init.push_back({initial_mps.get_tensor().real(), 0});
                else
                    init.push_back({initial_mps.get_tensor(), 0});
            } else {
                for(long n = 0; n < nev; n++) {
                    // Take the latest result with idx == n

                    // Start by collecting the results with the correct index
                    std::vector<std::reference_wrapper<const opt_mps>> results_idx_n;
                    for(const auto &r : results) {
                        if(r.get_eigs_idx() == n) results_idx_n.emplace_back(r);
                    }
                    if(not results_idx_n.empty()) {
                        if constexpr(std::is_same_v<Scalar, real>) {
                            init.push_back({results_idx_n.back().get().get_tensor().real(), n});
                        } else {
                            init.push_back({results_idx_n.back().get().get_tensor(), n});
                        }
                    }
                }
            }
            if(init.size() > safe_cast<size_t>(nev)) throw except::logic_error("Found too many initial guesses");
            return init;
        }

        template<typename Scalar>
        double get_largest_eigenvalue_hamiltonian_squared(const TensorsFinite &tensors) {
            const auto &env2                = tensors.get_multisite_env_var_blk();
            auto        hamiltonian_squared = MatVecMPO<Scalar>(env2.L, env2.R, tensors.get_multisite_mpo_squared());
            tools::log->trace("Finding largest-magnitude eigenvalue");
            eig::solver solver; // Define a solver just to find the maximum eigenvalue
            solver.config.tol             = settings::precision::eigs_tol_min;
            solver.config.maxIter         = 200;
            solver.config.maxNev          = 1;
            solver.config.maxNcv          = 16;
            solver.config.compute_eigvecs = eig::Vecs::OFF;
            solver.config.sigma           = std::nullopt;
            solver.config.ritz            = eig::Ritz::LM;
            solver.setLogLevel(2);
            solver.eigs(hamiltonian_squared);
            return eig::view::get_eigvals<double>(solver.result).cwiseAbs().maxCoeff();
        }

    }

    template<typename MatVecType>
    void eigs_folded_spectrum_executor(eig::solver &solver, MatVecType &hamiltonian_squared, const TensorsFinite &tensors, const opt_mps &initial_mps,
                                       std::vector<opt_mps> &results, const OptMeta &meta) {
        using Scalar = typename MatVecType::Scalar;
        if(std::is_same_v<Scalar, cplx> and meta.optType == OptType::REAL)
            throw except::logic_error("eigs_folded_spectrum_executor error: Mixed Scalar:cplx with OptType::REAL");
        if(std::is_same_v<Scalar, real> and meta.optType == OptType::CPLX)
            throw except::logic_error("eigs_folded_spectrum_executor error: Mixed Scalar:real with OptType::CPLX");

        solver.config.primme_effective_ham_sq = &hamiltonian_squared;
        hamiltonian_squared.reset();

        tools::log->trace("eigs_folded_spectrum_executor: Defining the Hamiltonian-squared matrix-vector product");

        if(solver.config.lib == eig::Lib::ARPACK) {
            if(tensors.model->has_compressed_mpo_squared()) throw std::runtime_error("optimize_folded_spectrum_eigs with ARPACK requires non-compressed MPO²");
            if(not solver.config.ritz) solver.config.ritz = eig::Ritz::LM;
            if(not solver.config.sigma)
                solver.config.sigma = folded_spectrum::get_largest_eigenvalue_hamiltonian_squared<Scalar>(tensors) + 1.0; // Add one to shift enough
        } else if(solver.config.lib == eig::Lib::PRIMME) {
            if(not solver.config.ritz) solver.config.ritz = eig::Ritz::SM;
            if(solver.config.sigma and solver.config.sigma.value() != 0.0 and tensors.model->has_compressed_mpo_squared())
                throw except::logic_error("optimize_folded_spectrum_eigs with PRIMME with sigma requires non-compressed MPO²");
        }

        tools::log->debug("eigs_folded_spectrum_executor: Solving [H²x=λx] {} {} | ritz {} | shifts {} | maxIter {} | tol {:.2e} | init on | size {} | mps {}",
                          eig::LibToString(solver.config.lib), eig::RitzToString(solver.config.ritz), solver.config.sigma, solver.config.primme_targetShifts,
                          solver.config.maxIter.value(), solver.config.tol.value(), hamiltonian_squared.rows(), hamiltonian_squared.get_shape_mps());

        auto init =
            folded_spectrum::get_initial_guess_mps<Scalar>(initial_mps, results, solver.config.maxNev.value()); // Init holds the data in memory for this scope
        for(auto &i : init) solver.config.initial_guess.push_back({i.mps.data(), i.idx});
        solver.eigs(hamiltonian_squared);
        internal::extract_results(tensors, initial_mps, meta, solver, results, false);
    }

    template<typename Scalar>
    void eigs_manager_folded_spectrum(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, const OptMeta &meta) {
        eig::solver solver;
        auto       &cfg           = solver.config;
        cfg.loglevel              = 2;
        cfg.compute_eigvecs       = eig::Vecs::ON;
        cfg.tol                   = meta.eigs_tol.value_or(settings::precision::eigs_tol_min); // 1e-12 is good. This Sets "eps" in primme, see link above.
        cfg.maxIter               = meta.eigs_iter_max.value_or(settings::precision::eigs_iter_max);
        cfg.maxNev                = meta.eigs_nev.value_or(1);
        cfg.maxNcv                = meta.eigs_ncv.value_or(settings::precision::eigs_ncv);
        cfg.maxTime               = meta.eigs_time_max.value_or(2 * 60 * 60); // Two hours default
        cfg.primme_minRestartSize = meta.primme_minRestartSize;
        cfg.primme_maxBlockSize   = meta.primme_maxBlockSize;
        cfg.primme_locking        = 0;

        cfg.lib           = eig::Lib::PRIMME;
        cfg.primme_method = eig::stringToMethod(meta.primme_method);
        cfg.tag += meta.label;
        switch(meta.optRitz) {
            case OptRitz::SR: cfg.ritz = eig::Ritz::primme_smallest; break;
            case OptRitz::LR: cfg.ritz = eig::Ritz::primme_largest; break;
            case OptRitz::LM: cfg.ritz = eig::Ritz::primme_largest_abs; break;
            case OptRitz::SM: {
                cfg.ritz                = eig::Ritz::primme_closest_abs; // H² is positive definite!
                cfg.primme_projection   = meta.primme_projection.value_or("primme_proj_refined");
                cfg.primme_targetShifts = {meta.eigv_target.value_or(0.0)};
                break;
            }
            default: throw except::logic_error("undhandled ritz: {}", enum2sv(meta.optRitz));
        }
        cfg.primme_preconditioner = folded_spectrum::preconditioner_jacobi<Scalar>;
        cfg.jcbMaxBlockSize       = meta.eigs_jcbMaxBlockSize;

        const auto &mpos                  = tensors.get_model().get_mpo_active();
        const auto &envv                  = tensors.get_edges().get_var_active();
        auto        hamiltonian_squared   = MatVecMPOS<Scalar>(mpos, envv);
        hamiltonian_squared.factorization = eig::Factorization::LLT; // H² is positive definite so this should work for the preconditioner
        eigs_folded_spectrum_executor(solver, hamiltonian_squared, tensors, initial_mps, results, meta);
    }

    opt_mps internal::optimize_folded_spectrum(const TensorsFinite &tensors, const opt_mps &initial_mps, [[maybe_unused]] const AlgorithmStatus &status,
                                               OptMeta &meta) {
        if(meta.optSolver == OptSolver::EIG) return optimize_folded_spectrum_eig(tensors, initial_mps, status, meta);

        using namespace internal;
        using namespace settings::precision;
        initial_mps.validate_initial_mps();
        reports::eigs_add_entry(initial_mps, spdlog::level::debug);

        auto                 t_var = tid::tic_scope("eigs-xdmrg", tid::level::higher);
        std::vector<opt_mps> results;
        switch(meta.optType) {
            case OptType::REAL: eigs_manager_folded_spectrum<real>(tensors, initial_mps, results, meta); break;
            case OptType::CPLX: eigs_manager_folded_spectrum<cplx>(tensors, initial_mps, results, meta); break;
        }
        auto t_post = tid::tic_scope("post");
        if(results.empty()) {
            meta.optExit = OptExit::FAIL_ERROR;
            return initial_mps; // The solver failed
        }

        // Sort results
        if(results.size() > 1) { std::sort(results.begin(), results.end(), internal::Comparator(meta)); }
        for(const auto &result : results) reports::eigs_add_entry(result, spdlog::level::debug);

        return results.front();
    }
}