
#include "../opt_meta.h"
#include "../opt_mps.h"
#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "general/iter.h"
#include "math/eig.h"
#include "math/eig/matvec/matvec_mpo.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/state/StateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt/opt-internal.h"
#include "tools/finite/opt/report.h"
#include <primme/primme.h>
// Temporary
//#include <h5pp/h5pp.h>

namespace tools::finite::opt {

    template<typename Scalar>
    struct opt_init_t {
        Eigen::Tensor<Scalar, 3> mps = {};
        long                     idx = 0;
    };

    template<typename Scalar>
    Eigen::Tensor<Scalar, 3> get_initial_guess(const opt_mps &initial_mps, const std::vector<opt_mps> &results) {
        if(results.empty()) {
            if constexpr(std::is_same_v<Scalar, real>)
                return initial_mps.get_tensor().real();
            else
                return initial_mps.get_tensor();
        } else {
            // Return whichever of initial_mps or results that has the lowest variance or gradient
            auto it = std::min_element(results.begin(), results.end(), internal::comparator::gradient);
            if(it == results.end()) return get_initial_guess<Scalar>(initial_mps, {});

            if(it->get_grad_max() < initial_mps.get_grad_max()) {
                tools::log->debug("Previous result is a good initial guess: {} | var {:8.2e}  ∇fᵐᵃˣ {:8.2e}", it->get_name(), it->get_variance(),
                                  it->get_grad_max());
                return get_initial_guess<Scalar>(*it, {});
            } else
                return get_initial_guess<Scalar>(initial_mps, {});
        }
    }

    template<typename Scalar>
    std::vector<opt_init_t<Scalar>> get_initial_guesses(const opt_mps &initial_mps, const std::vector<opt_mps> &results, long nev) {
        std::vector<opt_init_t<Scalar>> init;
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
                //
                //
                //                // Return whichever of initial_mps or results that has the lowest variance or gradient
                //
                //
                //
                //                long min_idx = min_gradient_idx(results, n);
                //                if(min_idx >= 0) {
                //                    auto &res = results.at(static_cast<size_t>(min_idx));
                //                    tools::log->debug("Found good initial guess for nev {}: idx {} {:<34} | lg var {:.8f}  ∇fᵐᵃˣ {:8.2e}", n, min_idx,
                //                    res.get_name(),
                //                                      std::log10(res.get_variance()), res.get_grad_max());
                //
                //                    if constexpr(std::is_same_v<Scalar, real>)
                //                        init.push_back({res.get_tensor().real(), idx});
                //                    else
                //                        init.push_back({res.get_tensor(), idx});
                //                }
            }
        }
        if(init.size() > static_cast<size_t>(nev)) throw except::logic_error("Found too many initial guesses");
        return init;
    }

    template<typename Scalar, typename MatVecMPO<Scalar>::DecompMode mode>
    void preconditioner(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr) {
        if(x == nullptr) return;
        if(y == nullptr) return;
        if(primme == nullptr) return;
        const auto H_ptr = static_cast<MatVecMPO<Scalar> *>(primme->matrix);
        H_ptr->decomp    = mode;
        H_ptr->FactorOP();
        H_ptr->MultOPv(x, ldx, y, ldy, blockSize, primme, ierr);
    }

    void RnormReadConvTest([[maybe_unused]] double *eval, [[maybe_unused]] void *evec, double *rNorm, int *isconv, struct primme_params *primme, int *ierr) {
        if(rNorm == nullptr) return;
        if(primme == nullptr) return;

        double problemNorm;
        if(not primme->massMatrixMatvec) {
            problemNorm = primme->aNorm > 0.0 ? primme->aNorm : primme->stats.estimateLargestSVal;
        } else {
            problemNorm = primme->aNorm > 0.0 && primme->invBNorm > 0.0 ? primme->aNorm * primme->invBNorm : primme->stats.estimateLargestSVal;
        }
        *ierr = 0;
    }

    template<typename MatrixProductType>
    void GradientConvTest(double *eval, void *evec, double *rNorm, int *isconv, struct primme_params *primme, int *ierr) {
        if(rNorm == nullptr) return;
        if(primme == nullptr) return;

        double problemNorm;
        if(not primme->massMatrixMatvec) {
            problemNorm = primme->aNorm > 0.0 ? primme->aNorm : primme->stats.estimateLargestSVal;
        } else {
            problemNorm = primme->aNorm > 0.0 && primme->invBNorm > 0.0 ? primme->aNorm * primme->invBNorm : primme->stats.estimateLargestSVal;
        }
        double prec         = primme->eps * problemNorm;
        int    default_crit = *rNorm < prec;

        *isconv = default_crit;
        *ierr   = 0;

        if(eval == nullptr) return;
        if(evec == nullptr) return;
        if(primme->convtest == nullptr) return;
        auto &solver = *static_cast<eig::solver *>(primme->convtest);
        auto &config = solver.config;
        auto &result = solver.result;

        // Store rNorm
        result.meta.last_res_norm = *rNorm;

        // Terminate if its taking too long
        if(config.maxTime.has_value() and primme->stats.elapsedTime > config.maxTime.value()) {
            solver.log->warn("primme: max time has been exeeded: {:.2f}", config.maxTime.value());
            primme->maxMatvecs = 0;
        }

        if constexpr(MatrixProductType::storage == eig::Storage::MPS) {
            if(config.primme_effective_ham == nullptr) return;
            if(config.primme_effective_ham_sq == nullptr) return;
            long   iter_since_last = primme->stats.numMatvecs - result.meta.last_grad_iter;
            double time_since_last = primme->stats.elapsedTime - result.meta.last_grad_time;

            double grad_tol  = config.primme_grad_tol ? config.primme_grad_tol.value() : 1e-6;
            long   grad_iter = config.primme_grad_iter ? config.primme_grad_iter.value() : 100;
            double grad_time = config.primme_grad_time ? config.primme_grad_time.value() : 5.0;

            // If we do inner iterations we better check gradient convergence every outer iteration instead, so we override the
            // iter/time periods given in config.
            if(primme->correctionParams.maxInnerIterations != 0) {
                grad_iter = 0;
                grad_time = 0;
            }

            if(iter_since_last >= grad_iter or time_since_last >= grad_time) {
                auto t_grad  = tid::tic_token("primme.grad");
                auto H_ptr   = static_cast<MatrixProductType *>(config.primme_effective_ham);
                auto H2_ptr  = static_cast<MatrixProductType *>(config.primme_effective_ham_sq);
                using Scalar = typename MatrixProductType::Scalar;
                using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
                if(H2_ptr->rows() != H_ptr->rows()) throw except::logic_error("H2 and H size mismatch");
                H_ptr->compress();
                long mps_size = static_cast<long>(H2_ptr->rows());
                auto v        = Eigen::Map<Vector>(static_cast<Scalar *>(evec), mps_size);

                Vector Hv(mps_size), H2v(mps_size);

                H_ptr->MultAx(v.data(), Hv.data());
                H2_ptr->MultAx(v.data(), H2v.data());

                if(H_ptr->template get_shift<Scalar>() != 0.0) Hv += v * H_ptr->template get_shift<Scalar>();
                if(H2_ptr->template get_shift<Scalar>() != 0.0) H2v += v * H2_ptr->template get_shift<Scalar>();
                auto vHv      = v.dot(Hv);
                auto vH2v     = v.dot(H2v);
                auto norm_1   = 1.0 / v.norm();
                auto pref     = std::is_same_v<Scalar, cplx> ? 2.0 : 1.0; // Factor 2 for complex
                auto grad     = pref * norm_1 * (H2v - 2.0 * vHv * Hv - (vH2v - 2.0 * vHv * vHv) * v);
                auto grad_max = grad.template lpNorm<Eigen::Infinity>();

                // Store gradient info
                result.meta.last_grad_max  = grad_max;
                result.meta.last_grad_iter = primme->stats.numMatvecs;
                result.meta.last_grad_time = primme->stats.elapsedTime;

                *isconv = std::max<int>(*isconv, grad_max < grad_tol);
                *ierr   = 0;
                //            if(*isconv != 0)
                solver.log->debug(
                    FMT_STRING("ops {:< 5} | mv {:<5} | iter {:<4} | λ {:20.16f} | res {:8.2e} | ∇fᵐᵃˣ {:8.2e} | time {:8.2f} s | dt {:8.2e} s/op"), H2_ptr->mv,
                    primme->stats.numMatvecs, primme->stats.numOuterIterations, primme->stats.estimateMinEVal, *rNorm, grad_max, primme->stats.elapsedTime,
                    primme->stats.timeMatvec / primme->stats.numMatvecs);
            }
        }
    }

    bool try_harder(const std::vector<opt_mps> &results, [[maybe_unused]] const OptMeta &meta, spdlog::level::level_enum level = spdlog::level::debug) {
        if(results.empty()) {
            tools::log->debug("Try harder: true | first run");
            return true;
        }

        //        static constexpr auto full_str     = std::array<std::string_view, 3>{"lapack", "dsyevd", "zheevd"};
        //        bool                  is_full_diag = std::find(full_str.begin(), full_str.end(), results.back().get_name()) != full_str.end();
        const auto &name  = results.back().get_name();
        bool is_full_diag = name.find("lapack") != std::string::npos or name.find("dsyevd") != std::string::npos or name.find("dsyevx") != std::string::npos or
                            name.find("zheevd") != std::string::npos;
        if(is_full_diag) {
            tools::log->log(level, "Try harder: false | full diagonalization is good enough");
            return false; // Full diag. You are finished.
        }
        // It's important that we only consider the "eigenvector 0" results, i.e. idx == 0.
        std::vector<std::reference_wrapper<const opt_mps>> results_eig0;
        for(const auto &r : results) {
            if(r.get_eigs_idx() == 0) results_eig0.emplace_back(r);
        }

        const auto &back0 = results_eig0.back().get();

        // When there is a degeneracy, we may get a situation like this:
        //
        // E/L                    λ                      σ²H      overlap
        // +5.565233108953220     +0.000000000001958     2.74e-12 0.695437332511797
        // +5.565233122209704     +0.000000000003793     2.75e-12 0.718586749494365
        //
        // We are really looking for the state in the second line, but since nev == 1 we don't see the next state.
        // The only way to detect this scenario is to check if the overlap is below 0.7071...
        // If that is true, we should get the next state with nev == 2 and check its overlap too.
        // If their variance is about the same we know that we have met this scenario and should pick
        // the one with highest overlap.
        // TODO: One could perhaps try SVD on each state to check if the entropy changes drastically?
        if(back0.get_overlap() <= std::sqrt(0.5)) {
            tools::log->log(level, "Try harder: true | nearly degenerate");
            return true;
        }

        // TODO: The test below is wrong.  eigs_tol is eps, but the convergence condition is actually res <= eps * aNorm * invBnorm
        //        if(back0.get_eigs_rnorm() <= back0.get_eigs_tol()) { // 1.1 to get a non-exact match as well
        //            tools::log->debug("Try harder: false | residual_norm tolerance reached: {:8.2e}", back0.get_eigs_tol());
        //            return false;
        //        }
        bool tryharder = back0.get_grad_max() >= back0.get_grad_tol();
        tools::log->log(level, "Try harder: {} | grad norm {:8.2e} threshold = {:<8.2e}", tryharder, back0.get_grad_max(), back0.get_grad_tol());
        return tryharder;
    }

    template<typename Scalar>
    double get_largest_eigenvalue_hamiltonian_squared(const TensorsFinite &tensors) {
        const auto &env2                = tensors.get_multisite_env_var_blk();
        auto        hamiltonian_squared = MatVecMPO<Scalar>(env2.L, env2.R, tensors.get_multisite_mpo_squared());
        tools::log->trace("Finding largest-magnitude eigenvalue");
        eig::solver solver; // Define a solver just to find the maximum eigenvalue
        solver.config.tol             = settings::precision::eigs_tolerance;
        solver.config.maxIter         = 200;
        solver.config.maxNev          = 1;
        solver.config.maxNcv          = 16;
        solver.config.compute_eigvecs = eig::Vecs::OFF;
        solver.config.sigma           = std::nullopt;
        solver.config.compress        = settings::precision::use_compressed_mpo_squared_otf;
        solver.config.ritz            = eig::Ritz::LM;
        solver.setLogLevel(2);
        solver.eigs(hamiltonian_squared);
        return eig::view::get_eigvals<double>(solver.result).cwiseAbs().maxCoeff();
    }

    template<typename Scalar>
    void eigs_variance_executor(eig::solver &solver, MatVecMPO<Scalar> &hamiltonian_squared, const TensorsFinite &tensors, const opt_mps &initial_mps,
                                std::vector<opt_mps> &results, const OptMeta &meta) {
        if(std::is_same_v<Scalar, cplx> and meta.optType == OptType::REAL)
            throw except::logic_error("eigs_launcher error: Mixed Scalar:cplx with OptType::REAL");
        if(std::is_same_v<Scalar, real> and meta.optType == OptType::CPLX)
            throw except::logic_error("eigs_launcher error: Mixed Scalar:real with OptType::CPLX");

        const auto       &mpo = tensors.get_multisite_mpo();
        const auto       &env = tensors.get_multisite_env_ene_blk();
        MatVecMPO<Scalar> hamiltonian(env.L, env.R, mpo);
        solver.config.primme_effective_ham    = &hamiltonian;
        solver.config.primme_effective_ham_sq = &hamiltonian_squared;

        hamiltonian_squared.reset();
        auto init = get_initial_guesses<Scalar>(initial_mps, results, solver.config.maxNev.value()); // Init holds the data in memory for this scope
        for(auto &i : init) solver.config.initial_guess.push_back({i.mps.data(), i.idx});
        tools::log->trace("Defining energy-shifted Hamiltonian-squared matrix-vector product");

        if(solver.config.lib == eig::Lib::ARPACK) {
            if(tensors.model->is_compressed_mpo_squared()) throw std::runtime_error("eigs_optimize_variance with ARPACK requires non-compressed MPO²");
            if(not solver.config.ritz) solver.config.ritz = eig::Ritz::LM;
            if(not solver.config.sigma)
                solver.config.sigma = get_largest_eigenvalue_hamiltonian_squared<Scalar>(tensors) + 1.0; // Add one just to make sure we shift enough
            tools::log->debug("Finding excited state as minimum of [(H-E)²-σ] | σ = {:.16f} | arpack {} | init on | dims {} = {}",
                              std::real(solver.config.sigma.value()), eig::RitzToString(solver.config.ritz.value()), hamiltonian_squared.get_shape_mps(),
                              hamiltonian_squared.rows());
            solver.eigs(hamiltonian_squared);
        } else if(solver.config.lib == eig::Lib::PRIMME) {
            if(not solver.config.ritz) solver.config.ritz = eig::Ritz::SA;
            if(solver.config.sigma and tensors.model->is_compressed_mpo_squared())
                throw except::logic_error("eigs_optimize_variance with PRIMME with given sigma requires non-compressed MPO²");
            if(solver.config.sigma)
                tools::log->debug("Finding excited state as minimum of [(H-E)²-σ] | σ = {:.16f} | primme {} | maxIter {} | init on | dims {} = {}",
                                  std::real(solver.config.sigma.value()), eig::RitzToString(solver.config.ritz.value()), solver.config.maxIter.value(),
                                  hamiltonian_squared.get_shape_mps(), hamiltonian_squared.rows());

            else
                tools::log->debug("Finding excited state as minimum of [(H-E)²] | primme {} | maxIter {} | init on | dims {} = {}",
                                  eig::RitzToString(solver.config.ritz.value()), solver.config.maxIter.value(), hamiltonian_squared.get_shape_mps(),
                                  hamiltonian_squared.rows());
            solver.eigs(hamiltonian_squared);
        }

        internal::eigs_extract_results(tensors, initial_mps, meta, solver, results, false);
    }

    template<typename Scalar>
    void eigs_manager(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, const OptMeta &meta) {
        std::vector<eig::settings> configs(1);
        // https://www.cs.wm.edu/~andreas/software/doc/appendix.html#c.primme_params.eps
        configs[0].tol             = 1e-12; // 1e-12 is good. This Sets "eps" in primme, see link above.
        configs[0].maxIter         = 1000;
        configs[0].maxNev          = 1;
        configs[0].maxNcv          = 4;
        configs[0].compress        = false;
        configs[0].maxTime         = 2 * 60 * 60; // Two hours
        configs[0].lib             = eig::Lib::PRIMME;
        configs[0].ritz            = eig::Ritz::SA;
        configs[0].compute_eigvecs = eig::Vecs::ON;
        configs[0].loglevel        = 2;
        configs[0].primme_method   = eig::PrimmeMethod::PRIMME_GD_plusK; // eig::PrimmeMethod::PRIMME_JDQMR;
        //        configs[0].primme_preconditioner = simps_preconditioner<MatVecMPO<Scalar>>;

        // Overrides from default
        if(meta.compress_otf) configs[0].compress = meta.compress_otf;
        if(meta.eigs_max_tol) configs[0].tol = meta.eigs_max_tol;
        if(meta.eigs_max_ncv) configs[0].maxNcv = meta.eigs_max_ncv;
        if(meta.eigs_max_iter) configs[0].maxIter = meta.eigs_max_iter;
        if(meta.eigs_grad_tol) configs[0].primme_grad_tol = meta.eigs_grad_tol;

        const auto                      &env2 = tensors.get_multisite_env_var_blk();
        std::optional<MatVecMPO<Scalar>> hamiltonian_squared;
        for(const auto &config : configs) {
            eig::solver solver;
            solver.config = config;
            if(not hamiltonian_squared) hamiltonian_squared = MatVecMPO<Scalar>(env2.L, env2.R, tensors.get_multisite_mpo_squared());
            eigs_variance_executor<Scalar>(solver, hamiltonian_squared.value(), tensors, initial_mps, results, meta);
            if(&config == &configs.back()) break;
            if(not try_harder(results, meta, spdlog::level::debug)) break;
        }
    }

    opt_mps internal::eigs_optimize_variance(const TensorsFinite &tensors, const opt_mps &initial_mps, [[maybe_unused]] const AlgorithmStatus &status,
                                             OptMeta &meta) {
        if(tensors.active_problem_size() <= settings::precision::max_size_full_diag) return internal::eig_optimize_variance(tensors, initial_mps, status, meta);

        using namespace internal;
        using namespace settings::precision;
        initial_mps.validate_basis_vector();
        if(not tensors.model->is_shifted()) throw std::runtime_error("eigs_optimize_variance requires energy-shifted MPO²");
        reports::eigs_add_entry(initial_mps, spdlog::level::info);

        auto                 t_var = tid::tic_scope("eigs-var");
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

        if(results.size() >= 2) {
            std::sort(results.begin(), results.end(), internal::Comparator(meta)); // Smallest eigenvalue (i.e. variance) wins
        }

        for(const auto &mps : results) reports::eigs_add_entry(mps, spdlog::level::info);

        return results.front();
    }

}