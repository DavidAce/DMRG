
#include "../opt_meta.h"
#include "../opt_mps.h"
#include "opt-internal.h"
#include "report.h"
#include <algorithms/AlgorithmStatus.h>
#include <config/settings.h>
#include <general/iter.h>
#include <math/eig.h>
#include <math/eig/matvec/matvec_mpo_eigen.h>
#include <math/tenx.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/state/StateFinite.h>
#include <tensors/TensorsFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>

void tools::finite::opt::internal::krylov_extract_solutions(const TensorsFinite &tensors, const opt_mps &initial_mps, const eig::solver &solver,
                                                            std::vector<opt_mps> &results, const OptMeta &meta, bool converged_only) {
    auto dims_mps = initial_mps.get_tensor().dimensions();
    if(solver.result.meta.eigvals_found and solver.result.meta.eigvecsR_found) {
        auto eigvecs = eig::view::get_eigvecs<cplx>(solver.result, eig::Side::R, converged_only);
        auto eigvals = eig::view::get_eigvals<real>(solver.result, converged_only);
        if(eigvecs.cols() == eigvals.size()) {
            for(long idx = 0; idx < eigvals.size(); idx++) {
                // It's important to normalize the eigenvectors - they are not always well normalized when we get them from the eig::solver
                auto eigvec_i = tenx::TensorCast(eigvecs.col(idx).normalized(), dims_mps);
                auto overlap  = std::abs(initial_mps.get_vector().dot(eigvecs.col(idx)));
                auto energy   = tools::finite::measure::energy(eigvec_i, tensors);
                auto eigval   = energy - initial_mps.get_energy_reduced();
                auto variance = tools::finite::measure::energy_variance(eigvec_i, tensors);
                results.emplace_back(fmt::format("{:<8} eigenvector {}", solver.config.tag, idx), eigvec_i, tensors.active_sites, eigval,
                                     initial_mps.get_energy_reduced(), variance, overlap, tensors.get_length());
                auto &mps = results.back();
                mps.set_time(solver.result.meta.time_total);
                mps.set_counter(static_cast<size_t>(solver.result.meta.counter));
                mps.set_iter(static_cast<size_t>(solver.result.meta.iter));
                mps.set_grad_norm(tools::finite::measure::grad_max_norm(eigvec_i, tensors));
                mps.is_basis_vector = true;
                mps.validate_candidate();
                mps.set_krylov_nev(solver.result.meta.nev_converged);
                mps.set_krylov_ncv(solver.result.meta.ncv);
                mps.set_krylov_tol(solver.result.meta.tol);
                mps.set_krylov_eigval(eigvals(idx));
                mps.set_krylov_ritz(solver.result.meta.ritz);
                mps.set_krylov_shift(solver.result.meta.sigma);
                mps.set_optmode(meta.optMode);
                mps.set_optspace(meta.optSpace);
                if(not solver.result.meta.residual_norms.empty()) mps.set_krylov_resid(solver.result.meta.residual_norms.at(idx));
            }
        }
    }
}

namespace tools::finite::opt::internal {
    // Make a handy variance comparator
    auto comp_variance = [](const opt_mps &lhs, const opt_mps &rhs) {
        return lhs.get_variance() < rhs.get_variance();
    };
    auto comp_overlap = [](const opt_mps &lhs, const opt_mps &rhs) {
        return lhs.get_overlap() > rhs.get_overlap();
    };
    auto comp_gradient = [](const opt_mps &lhs, const opt_mps &rhs) {
        return lhs.get_grad_norm() < rhs.get_grad_norm();
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
            auto it = std::min_element(results.begin(), results.end(), comp_gradient);
            if (it == results.end()) return get_initial_guess<Scalar>(initial_mps, {});

            if(it->get_grad_norm() < initial_mps.get_grad_norm()) {
                tools::log->debug("Previous result is a good initial guess: {} | log₁₀ var {:.8f}  |∇f|∞ {:8.2e}", it->get_name(),
                                  std::log10(it->get_variance()), it->get_grad_norm());
                return get_initial_guess<Scalar>(*it, {});
            } else
                return get_initial_guess<Scalar>(initial_mps, {});
        }
    }

    template<typename Scalar>
    MatVecMPOEigen<Scalar> get_matvec_hamiltonian_squared(const TensorsFinite &tensors) {
        const auto &env2 = tensors.get_multisite_var_blk();
        return MatVecMPOEigen<Scalar>(env2.L, env2.R, tensors.get_multisite_mpo_squared());
    }

    template<typename Scalar>
    double get_largest_eigenvalue_hamiltonian_squared(const TensorsFinite &tensors) {
        auto hamiltonian_squared = get_matvec_hamiltonian_squared<Scalar>(tensors);
        tools::log->trace("Finding largest-magnitude eigenvalue");
        eig::solver solver; // Define a solver just to find the maximum eigenvalue
        solver.config.tol             = settings::precision::eig_tolerance;
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
    void krylov_launcher(const TensorsFinite &tensors, const opt_mps &initial_mps, eig::solver &solver, std::vector<opt_mps> &results, const OptMeta &meta) {
        if(std::is_same_v<Scalar, cplx> and meta.optType == OptType::REAL)
            throw std::logic_error("krylov_launcher error: Mixed Scalar:cplx with OptType::REAL");
        if(std::is_same_v<Scalar, real> and meta.optType == OptType::CPLX)
            throw std::logic_error("krylov_launcher error: Mixed Scalar:real with OptType::CPLX");
        auto size = tensors.active_problem_size();

        if(size <= 1024) {
            tools::log->trace("Full diagonalization of (H-E)²");
            auto matrix       = tools::finite::opt::internal::get_multisite_hamiltonian_squared_matrix<double>(*tensors.model, *tensors.edges);
            solver.config.tag = "lapack";
            solver.eig(matrix.data(), matrix.rows());
        } else {
            auto init            = get_initial_guess<Scalar>(initial_mps, results);
            solver.config.residp = init.data();
            tools::log->trace("Defining reduced Hamiltonian-squared matrix-vector product");

            if(tensors.model->is_compressed_mpo_squared()) {
                tools::log->warn("Finding excited state using H² with ritz SM because all MPO²'s are compressed!");
                auto hamiltonian_squared = get_matvec_hamiltonian_squared<Scalar>(tensors);
                solver.config.ritz       = eig::Ritz::SM;
                solver.eigs(hamiltonian_squared);
            } else {
                if(solver.config.lib == eig::Lib::ARPACK) {
                    double largest_eigval      = get_largest_eigenvalue_hamiltonian_squared<Scalar>(tensors) + 1.0; // Add one just to make sure we shift enough
                    auto   hamiltonian_squared = get_matvec_hamiltonian_squared<Scalar>(tensors);
                    tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | arpack LM | init on ...",
                                      largest_eigval);
                    solver.config.sigma = largest_eigval;
                    solver.eigs(hamiltonian_squared);
                    tools::log->trace(
                        "Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | init on ... OK: iters {} | counter {} | "
                        "time {:.3f} | matprod {:.3f}s",
                        largest_eigval, solver.result.meta.iter, solver.result.meta.counter, solver.result.meta.time_total, solver.result.meta.time_matprod);
                } else if(solver.config.lib == eig::Lib::PRIMME) {
                    // Reset the matrix
                    tools::log->debug("Finding excited state using shifted operator [(H-E)²-λ], with λ = 1.0 | primme SA | init on ...");
                    auto hamiltonian_squared = get_matvec_hamiltonian_squared<Scalar>(tensors);
                    solver.config.sigma      = 1.0;
                    solver.eigs(hamiltonian_squared);
                }
            }
        }
        krylov_extract_solutions(tensors, initial_mps, solver, results, meta, false);
    }

    bool try_harder(const std::vector<opt_mps> & results){
        if(results.empty()) return true;
        if(results.back().get_name().find("lapack") != std::string::npos) return false; // Full diag. You are finished.

        auto it = std::min_element(results.begin(), results.end(), comp_gradient);
        if(it == results.end()) return true;
        if(it->get_grad_norm() > 1e+2) return true;  // Very large gradient, just try harder
        if(it->get_grad_norm() < 1e-0) return false; // We rarely improve variance once grad max norm is smaller than one
        if(results.size() >= 2){
            // In this case we can compare the gradient improvement between the last two runs.
            // If there hasn't been an improvement, we can assume it won't help to try anymore
            const auto & res_back0 = results.back();
            const auto & res_back1 = results.rbegin()[1];
            auto grad_diff = std::abs(res_back0.get_grad_norm()/res_back1.get_grad_norm());
            return grad_diff < 1e-1; // Try harder if the last run made a 90% reduction on the gradient
        }

        return it->get_grad_norm() > 1e-0; // We rarely improve variance once grad max norm is smaller than one
    }

}

tools::finite::opt::opt_mps tools::finite::opt::internal::krylov_optimization(const TensorsFinite &tensors, const opt_mps &initial_mps,
                                                                              [[maybe_unused]] const AlgorithmStatus &status, OptMeta &meta) {
    using namespace internal;
    using namespace settings::precision;
    initial_mps.validate_candidate();
    if(not tensors.model->is_reduced()) throw std::runtime_error("krylov_optimization requires energy-reduced MPO²");
    if(tensors.model->is_compressed_mpo_squared()) throw std::runtime_error("krylov_optimization requires non-compressed MPO²");

    auto t_var    = tid::tic_scope("krylov");
    auto dims_mps = initial_mps.get_tensor().dimensions();
    auto size     = initial_mps.get_tensor().size();

    tools::log->debug("Finding excited state: minimum eigenstate of (H-E)² | dims {} = {}", dims_mps, size);

    eig::solver solver_arpack;
    solver_arpack.config.tol             = settings::precision::eig_tolerance;
    solver_arpack.config.maxIter         = 400;
    solver_arpack.config.maxTime         = 30 * 60;
    solver_arpack.config.maxNev          = 1;
    solver_arpack.config.maxNcv          = 64;
    solver_arpack.config.compress        = settings::precision::use_compressed_mpo_squared_otf;
    solver_arpack.config.lib             = eig::Lib::ARPACK;
    solver_arpack.config.ritz            = eig::Ritz::LM;
    solver_arpack.config.compute_eigvecs = eig::Vecs::ON;
    solver_arpack.config.tag             = "arpack tol * 1e+0";
    solver_arpack.setLogLevel(2);

    eig::solver solver_arpack_tol1e6;
    solver_arpack_tol1e6.config = solver_arpack.config;
    solver_arpack_tol1e6.config.tol             = 1e6*settings::precision::eig_tolerance;
    solver_arpack_tol1e6.config.tag             = "arpack tol * 1e+6";



    eig::solver solver_primme;
    solver_primme.config.tol             = settings::precision::eig_tolerance;
    solver_primme.config.maxIter         = 80000;
    solver_primme.config.maxTime         = 60 * 60;
    solver_primme.config.maxNev          = 1;
    solver_primme.config.maxNcv          = 32; // 128
    solver_primme.config.compress        = settings::precision::use_compressed_mpo_squared_otf;
    solver_primme.config.lib             = eig::Lib::PRIMME;
    solver_primme.config.ritz            = eig::Ritz::SA;
    solver_primme.config.compute_eigvecs = eig::Vecs::ON;
    solver_primme.config.tag             = "primme tol * 1e+0";
    solver_primme.config.primme_method   = "PRIMME_GD_plusK";
    solver_primme.setLogLevel(2);

    eig::solver solver_primme_tol1e4;
    solver_primme_tol1e4.config = solver_primme.config;
    solver_primme_tol1e4.config.tol             = 1e4*settings::precision::eig_tolerance;
    solver_primme_tol1e4.config.tag             = "primme tol * 1e+4";

    eig::solver solver_primme_tol1e2;
    solver_primme_tol1e2.config = solver_primme.config;
    solver_primme_tol1e2.config.tol             = 1e2*settings::precision::eig_tolerance;
    solver_primme_tol1e2.config.tag             = "primme tol * 1e+2";






    std::vector<opt_mps> results;

#pragma message "testing primme"
    if(try_harder(results)){
        switch(meta.optType) {
            case OptType::REAL: krylov_launcher<real>(tensors, initial_mps, solver_arpack_tol1e6, results, meta); break;
            case OptType::CPLX: krylov_launcher<cplx>(tensors, initial_mps, solver_arpack_tol1e6, results, meta); break;
        }
    }


    if(try_harder(results)){
        switch(meta.optType) {
            case OptType::REAL: krylov_launcher<real>(tensors, initial_mps, solver_primme_tol1e4, results, meta); break;
            case OptType::CPLX: krylov_launcher<cplx>(tensors, initial_mps, solver_primme_tol1e4, results, meta); break;
        }
    }

    if(try_harder(results)){
        switch(meta.optType) {
            case OptType::REAL: krylov_launcher<real>(tensors, initial_mps, solver_primme_tol1e2, results, meta); break;
            case OptType::CPLX: krylov_launcher<cplx>(tensors, initial_mps, solver_primme_tol1e2, results, meta); break;
        }
    }


    if(try_harder(results)){
        switch(meta.optType) {
            case OptType::REAL: krylov_launcher<real>(tensors, initial_mps, solver_primme, results, meta); break;
            case OptType::CPLX: krylov_launcher<cplx>(tensors, initial_mps, solver_primme, results, meta); break;
        }
    }


    t_var.toc();

    if(results.empty()) {
        meta.optExit = OptExit::FAIL_ERROR;
        return initial_mps; // The solver failed
    }
    constexpr size_t max_print = settings::debug ? 32 : 8;
    for(const auto &[num, mps] : iter::enumerate(results)) {
        if(num >= max_print) break;
        reports::krylov_add_entry(mps);
    }

    if(results.size() >= 2) {
#pragma message "Sorting according to gradient"
        if(results.back().get_name().find("lapack") != std::string::npos)
            std::sort(results.begin(), results.end(), comp_variance);
        else if(meta.optMode == OptMode::VARIANCE) std::sort(results.begin(), results.end(), comp_gradient);
        else if(meta.optMode == OptMode::OVERLAP) std::sort(results.begin(), results.end(), comp_overlap);
    }
    return results.front();
}
