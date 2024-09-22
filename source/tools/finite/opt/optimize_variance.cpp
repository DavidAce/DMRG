
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
    void massMatrixMatvec(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr) {
        if(x == nullptr) return;
        if(y == nullptr) return;
        if(primme == nullptr) return;
        const auto H_ptr = static_cast<MatVecMPOS<Scalar> *>(primme->matrix);
        H_ptr->MultBx(x, ldx, y, ldy, blockSize, primme, ierr);
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
                tools::log->info("eval {:38.32f} ({:+8.3e}) | rnorm {:38.32f} ({:+8.3e}) | problemNorm {:.3e}", *eval, diff_evals_sum, *rNorm, diff_rnorms_sum,
                                 problemNorm);
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
    std::vector<opt_bond_init_t<Scalar>> get_initial_guess_bond(const opt_mps &initial_mps, const std::vector<opt_mps> &results, long nev) {
        std::vector<opt_bond_init_t<Scalar>> init;
        if(results.empty()) {
            if constexpr(std::is_same_v<Scalar, real>)
                init.push_back({initial_mps.get_bond().real(), 0});
            else
                init.push_back({initial_mps.get_bond(), 0});
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
                        init.push_back({results_idx_n.back().get().get_bond().real(), n});
                    } else {
                        init.push_back({results_idx_n.back().get().get_bond(), n});
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

    template<typename MatVecType>
    void eigs_variance_executor(eig::solver &solver, MatVecType &hamiltonian_squared, const TensorsFinite &tensors, const opt_mps &initial_mps,
                                std::vector<opt_mps> &results, const OptMeta &meta) {
        using Scalar = typename MatVecType::Scalar;
        if(std::is_same_v<Scalar, cplx> and meta.optType == OptType::REAL)
            throw except::logic_error("eigs_variance_executor error: Mixed Scalar:cplx with OptType::REAL");
        if(std::is_same_v<Scalar, real> and meta.optType == OptType::CPLX)
            throw except::logic_error("eigs_variance_executor error: Mixed Scalar:real with OptType::CPLX");

        // const auto       &mpo = tensors.get_multisite_mpo();
        // const auto       &env = tensors.get_multisite_env_ene_blk();
        // MatVecMPO<Scalar> hamiltonian(env.L, env.R, mpo);
        // solver.config.primme_effective_ham    = &hamiltonian;
        solver.config.primme_effective_ham_sq = &hamiltonian_squared;
        hamiltonian_squared.reset();

        tools::log->trace("eigs_variance_executor: Defining the Hamiltonian-squared matrix-vector product");

        if(solver.config.lib == eig::Lib::ARPACK) {
            if(tensors.model->has_compressed_mpo_squared()) throw std::runtime_error("optimize_variance_eigs with ARPACK requires non-compressed MPO²");
            if(not solver.config.ritz) solver.config.ritz = eig::Ritz::LM;
            if(not solver.config.sigma) solver.config.sigma = get_largest_eigenvalue_hamiltonian_squared<Scalar>(tensors) + 1.0; // Add one to shift enough
        } else if(solver.config.lib == eig::Lib::PRIMME) {
            if(not solver.config.ritz) solver.config.ritz = eig::Ritz::SM;
            if(solver.config.sigma and solver.config.sigma.value() != 0.0 and tensors.model->has_compressed_mpo_squared())
                throw except::logic_error("optimize_variance_eigs with PRIMME with sigma requires non-compressed MPO²");
        }
        tools::log->debug("Finding an eigenstate of [H²{}]{}{} | {} ritz {} | shifts: {} | maxIter {} | tol {:.2e} | init on | size {} | mps {}",
                          solver.config.sigma ? "-σ" : "", solver.config.shift_invert == eig::Shinv::ON ? "⁻¹" : "",
                          solver.config.sigma ? fmt::format(" | σ = {:.16f}", solver.config.sigma->real()) : "", eig::LibToString(solver.config.lib),
                          eig::RitzToString(solver.config.ritz), solver.config.primme_targetShifts, solver.config.maxIter.value(), solver.config.tol.value(),
                          hamiltonian_squared.rows(), hamiltonian_squared.get_shape_mps());

        if constexpr(sfinae::is_any_v<MatVecType, MatVecZero<real>, MatVecZero<cplx>>) {
            auto init = get_initial_guess_bond<Scalar>(initial_mps, results, solver.config.maxNev.value()); // Init holds the data in memory for this scope
            for(auto &i : init) solver.config.initial_guess.push_back({i.bond.data(), i.idx});
            solver.eigs(hamiltonian_squared);
            internal::extract_results(tensors, initial_mps, meta, solver, results, false);
        } else {
            auto init = get_initial_guess_mps<Scalar>(initial_mps, results, solver.config.maxNev.value()); // Init holds the data in memory for this scope
            for(auto &i : init) solver.config.initial_guess.push_back({i.mps.data(), i.idx});
            solver.eigs(hamiltonian_squared);
            internal::extract_results(tensors, initial_mps, meta, solver, results, false);
        }
    }

    template<typename Scalar>
    void eigs_manager(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, const OptMeta &meta) {
        eig::solver solver;
        auto       &cfg     = solver.config;
        cfg.loglevel        = 1;
        cfg.compute_eigvecs = eig::Vecs::ON;
        cfg.tol             = meta.eigs_tol.value_or(settings::precision::eigs_tol_min); // 1e-12 is good. This Sets "eps" in primme, see link above.
        // cfg.tol_other             = std::max(1e-10, cfg.tol.value());                          // We do not need high precision on subsequent solutions
        cfg.maxIter               = meta.eigs_iter_max.value_or(settings::precision::eigs_iter_max);
        cfg.maxNev                = meta.eigs_nev.value_or(1);
        cfg.maxNcv                = meta.eigs_ncv.value_or(settings::precision::eigs_ncv);
        cfg.maxTime               = meta.eigs_time_max.value_or(2 * 60 * 60); // Two hours default
        cfg.primme_minRestartSize = meta.primme_minRestartSize;
        cfg.primme_maxBlockSize   = meta.primme_maxBlockSize;
        cfg.primme_locking        = 0;
        // cfg.primme_convTestFun    = convTestFun<Scalar>;

        cfg.lib           = eig::Lib::PRIMME;
        cfg.primme_method = eig::stringToMethod(meta.primme_method);
        cfg.tag += meta.label;
        switch(meta.optRitz) {
            case OptRitz::SR: cfg.ritz = eig::Ritz::primme_smallest; break;
            case OptRitz::LR: cfg.ritz = eig::Ritz::primme_largest; break;
            case OptRitz::SM: {
                cfg.ritz                = eig::Ritz::primme_closest_abs; // H² is positive definite!
                cfg.primme_targetShifts = {meta.eigv_target.value_or(0.0)};
                cfg.primme_projection   = meta.primme_projection.value_or("primme_proj_refined");
                // if(cfg.primme_projection.value() != "primme_proj_default") {
                //     cfg.ritz                = eig::Ritz::primme_closest_abs;
                //     cfg.primme_targetShifts = {meta.eigv_target.value_or(0.0)};
                // }
                break;
            }
            default: throw except::logic_error("undhandled ritz: {}", enum2sv(meta.optRitz));
        }

        // cfg.primme_preconditioner = preconditioner_jacobi<Scalar>;
        // cfg.jcbMaxBlockSize       = 1;//meta.eigs_jcbMaxBlockSize;

        cfg.primme_massMatrixMatvec = massMatrixMatvec<Scalar>;
        cfg.primme_projection       = meta.primme_projection.value_or("primme_proj_RR");
        cfg.primme_targetShifts.clear();
        cfg.ritz                = eig::Ritz::primme_largest_abs; // H² is positive definite!
        // Overrides from default

        const auto &state = *tensors.state;
        const auto &model = *tensors.model;
        const auto &edges = *tensors.edges;
        const auto &mpss  = state.get_mps_active();
        const auto &mpos  = model.get_mpo_active();
        const auto &envv  = edges.get_var_active();
        const auto &enve  = edges.get_ene_active();
        if(meta.optAlgo == OptAlgo::DIRECTZ) {
            auto hamiltonian_squared = MatVecZero<Scalar>(mpss, mpos, envv);
            // auto matrix = hamiltonian_squared.get_matrix();
            // eig::solver sol;
            // sol.eig<eig::Form::SYMM>(matrix.data(), matrix.rows(), eig::Vecs::ON);
            // auto eval = eig::view::get_eigvals<real>(sol.result);
            // fmt::print("{:^9} {:<20}\n", " ", "H_zero eigv");
            // for(long idx = 0; idx < eval.size(); ++idx) {
            //     // if(std::abs(evals1[idx]) > 1.1) continue;
            //     fmt::print("idx {:2}: {:20.16f}\n", idx, eval[idx]);
            // }

            eigs_variance_executor(solver, hamiltonian_squared, tensors, initial_mps, results, meta);
        } else {
            // auto hamiltonian_squared = MatVecMPOS<Scalar>(mpos, envv, enve);
            auto hamiltonian_squared = MatVecMPOS<Scalar>(mpos, enve, envv);
            if constexpr(hamiltonian_squared.can_shift_invert) {
                if(initial_mps.get_tensor().size() <= settings::precision::eigs_max_size_shift_invert) {
                    cfg.shift_invert                  = eig::Shinv::ON;
                    cfg.sigma                         = cplx(0.0, 0.0);
                    cfg.primme_targetShifts           = {};
                    cfg.ritz                          = eig::Ritz::primme_largest_abs;
                    hamiltonian_squared.factorization = eig::Factorization::LLT;
                }
            }
            eigs_variance_executor(solver, hamiltonian_squared, tensors, initial_mps, results, meta);
        }
    }

    opt_mps internal::optimize_variance_eigs(const TensorsFinite &tensors, const opt_mps &initial_mps, [[maybe_unused]] const AlgorithmStatus &status,
                                             OptMeta &meta) {
        if(meta.optSolver == OptSolver::EIG) return optimize_variance_eig(tensors, initial_mps, status, meta);

        using namespace internal;
        using namespace settings::precision;
        // initial_mps.validate_basis_vector();
        initial_mps.validate_initial_mps();
        // if(not tensors.model->is_shifted()) throw std::runtime_error("optimize_variance_eigs requires energy-shifted MPO²");
        reports::eigs_add_entry(initial_mps, spdlog::level::debug);

        auto                 t_var = tid::tic_scope("eigs-var", tid::level::higher);
        std::vector<opt_mps> results;
        switch(meta.optType) {
            case OptType::REAL: eigs_manager<real>(tensors, initial_mps, results, meta); break;
            case OptType::CPLX: eigs_manager<cplx>(tensors, initial_mps, results, meta); break;
        }
        auto t_post = tid::tic_scope("post");
        if(results.empty()) {
            meta.optExit = OptExit::FAIL_ERROR;
            return initial_mps; // The solver failed
        }

        // auto metax2    = meta;
        // metax2.optCost = OptCost::ENERGY;
        // metax2.optRitz = OptRitz::SM;
        // #pragma message "remove energy opt trial"
        //         std::vector<opt_mps> results2;
        //         for(const auto &result : results) {
        //             if(result.get_rnorm_H() < 1e-5) {
        //                 results2.emplace_back(optimize_energy_eigs(tensors, result, status, metax2));
        //                 results2.back().set_name(fmt::format("{} energy", results2.back().get_name()));
        //                 results2.back().set_eigs_eigval(std::pow(results2.back().get_eigs_eigval(), 2.0)); // set the eigval to E² for easier comparison with
        //                 H². reports::eigs_add_entry(results2.back(), spdlog::level::debug);
        //                 // tools::finite::opt::reports::print_eigs_report();
        //             }
        //         }
        //         results.insert(results.end(), std::make_move_iterator(results2.begin()), std::make_move_iterator(results2.end()));

        // Sort results
        if(results.size() > 1) {
            std::sort(results.begin(), results.end(), internal::Comparator(meta)); // Smallest eigenvalue wins
            // Add some information about the other eigensolutions to the winner
            // We can use that to later calculate the gap and the approximate convergence rate
            for(size_t idx = 1; idx < results.size(); ++idx) { results.front().next_evs.emplace_back(results[idx]); }
        }
        for(const auto &result : results) reports::eigs_add_entry(result, spdlog::level::debug);

        // if(results.size() > 1) {
        //     // Remove solutions whose residual norm is too bad to be taken seriously
        //     std::remove_if(results.rbegin(), results.rend(), [&](opt_mps &r) -> bool {
        //         if(results.size() > 1 and r.get_rnorm() > std::max(1e-6, settings::precision::eigs_tol_max)) { return true; }
        //         return false;
        //     });
        // }

        return results.front();
    }
}