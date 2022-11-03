
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

namespace tools::finite::opt {

    template<typename Scalar>
    struct opt_mps_init_t {
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
    std::vector<opt_mps_init_t<Scalar>> get_initial_guesses(const opt_mps &initial_mps, const std::vector<opt_mps> &results, long nev) {
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

    template<typename Scalar, typename eig::Factorization factorization>
    void preconditioner(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr) {
        if(x == nullptr) return;
        if(y == nullptr) return;
        if(primme == nullptr) return;
        const auto H_ptr     = static_cast<MatVecMPO<Scalar> *>(primme->matrix);
        H_ptr->factorization = factorization;
        H_ptr->FactorOP();
        H_ptr->MultOPv(x, ldx, y, ldy, blockSize, primme, ierr);
    }

    template<typename Scalar>
    double get_largest_eigenvalue_hamiltonian_squared(const TensorsFinite &tensors) {
        const auto &env2                = tensors.get_multisite_env_var_blk();
        auto        hamiltonian_squared = MatVecMPO<Scalar>(env2.L, env2.R, tensors.get_multisite_mpo_squared());
        tools::log->trace("Finding largest-magnitude eigenvalue");
        eig::solver solver; // Define a solver just to find the maximum eigenvalue
        solver.config.tol             = settings::solver::eigs_tol_min;
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
            if(not solver.config.sigma) solver.config.sigma = get_largest_eigenvalue_hamiltonian_squared<Scalar>(tensors) + 1.0; // Add one to shift enough
        } else if(solver.config.lib == eig::Lib::PRIMME) {
            if(not solver.config.ritz) solver.config.ritz = eig::Ritz::SA;
            if(solver.config.sigma and solver.config.sigma.value() != 0.0 and tensors.model->is_compressed_mpo_squared())
                throw except::logic_error("eigs_optimize_variance with PRIMME with sigma requires non-compressed MPO²");
        }
        tools::log->debug("Finding excited state of operator [(H-E)²{}]{}{} | {} {} | maxIter {} | tol {:.2e} | init on | size {} | mps {} | mpo {}",
                          solver.config.sigma ? "-σ" : "", solver.config.shift_invert == eig::Shinv::ON ? "⁻¹" : "",
                          solver.config.sigma ? fmt::format(" | σ = {:.16f}", solver.config.sigma->real()) : "", eig::LibToString(solver.config.lib),
                          eig::RitzToString(solver.config.ritz), solver.config.maxIter.value(), solver.config.tol.value(), hamiltonian_squared.rows(),
                          hamiltonian_squared.get_shape_mps(), hamiltonian_squared.get_shape_mpo());

        solver.eigs(hamiltonian_squared);
        internal::eigs_extract_results(tensors, initial_mps, meta, solver, results, false);
    }

    template<typename Scalar>
    void eigs_manager(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, const OptMeta &meta) {
        eig::solver solver;
        auto       &cfg = solver.config;
        // https://www.cs.wm.edu/~andreas/software/doc/appendix.html#c.primme_params.eps
        cfg.tol             = 1e-12; // 1e-12 is good. This Sets "eps" in primme, see link above.
        cfg.maxIter         = 1000;
        cfg.maxNev          = 1;
        cfg.maxNcv          = 4;
        cfg.compress        = false;
        cfg.maxTime         = 2 * 60 * 60; // Two hours
        cfg.lib             = eig::Lib::PRIMME;
        cfg.ritz            = eig::Ritz::primme_smallest; // eig::Ritz::SA;
        cfg.compute_eigvecs = eig::Vecs::ON;
        cfg.loglevel        = 2;
        cfg.primme_method   = eig::PrimmeMethod::PRIMME_DYNAMIC; // eig::PrimmeMethod::PRIMME_JDQMR;

        //         Apply preconditioner if applicable, usually faster on small matrices
        //                if(initial_mps.get_tensor().size() > settings::solver::max_size_full_eigs and initial_mps.get_tensor().size() <= 8000)
        //                    configs[0].primme_preconditioner = preconditioner<Scalar, MatVecMPO<Scalar>::DecompMode::LLT>;
        //        if(initial_mps.get_tensor().size() > 8192 and initial_mps.get_tensor().size() <= 20000)
        //        configs[0].primme_preconditioner = preconditioner<Scalar, MatVecMPO<Scalar>::DecompMode::MATRIXFREE>;

        // Overrides from default
        if(meta.compress_otf) cfg.compress = meta.compress_otf;
        if(meta.eigs_tol) cfg.tol = meta.eigs_tol;
        if(meta.eigs_ncv) cfg.maxNcv = meta.eigs_ncv;
        if(meta.eigs_iter_max) cfg.maxIter = meta.eigs_iter_max;

        const auto &env2                = tensors.get_multisite_env_var_blk();
        auto        hamiltonian_squared = MatVecMPO<Scalar>(env2.L, env2.R, tensors.get_multisite_mpo_squared());
        if(initial_mps.get_tensor().size() <= settings::solver::max_size_shift_invert) {
            cfg.shift_invert                  = eig::Shinv::ON;
            cfg.sigma                         = 0.0;
            cfg.primme_target_shifts          = {};
            cfg.ritz                          = eig::Ritz::primme_largest_abs;
            hamiltonian_squared.factorization = eig::Factorization::LLT;
            hamiltonian_squared.set_readyCompress(tensors.model->is_compressed_mpo_squared());
        }
        eigs_variance_executor<Scalar>(solver, hamiltonian_squared, tensors, initial_mps, results, meta);
    }

    opt_mps internal::eigs_optimize_variance(const TensorsFinite &tensors, const opt_mps &initial_mps, [[maybe_unused]] const AlgorithmStatus &status,
                                             OptMeta &meta) {
        if(tensors.active_problem_size() <= settings::solver::max_size_full_eigs) return internal::eig_optimize_variance(tensors, initial_mps, status, meta);

        using namespace internal;
        using namespace settings::precision;
        initial_mps.validate_basis_vector();
        if(not tensors.model->is_shifted()) throw std::runtime_error("eigs_optimize_variance requires energy-shifted MPO²");
        reports::eigs_add_entry(initial_mps, spdlog::level::debug);

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

        for(const auto &mps : results) reports::eigs_add_entry(mps, spdlog::level::debug);

        return results.front();
    }

}