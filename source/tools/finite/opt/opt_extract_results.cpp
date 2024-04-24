#include "../opt_meta.h"
#include "../opt_mps.h"
#include "config/settings.h"
#include "math/eig.h"
#include "math/num.h"
#include "measure/MeasurementsTensorsFinite.h"
#include "tensors/state/StateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt/opt-internal.h"

void tools::finite::opt::internal::extract_results(const TensorsFinite &tensors, const opt_mps &initial_mps, const OptMeta &meta, const eig::solver &solver,
                                                   std::vector<opt_mps> &results, bool converged_only, std::optional<std::vector<long>> indices) {
    auto t_ext    = tid::tic_scope("extract");
    auto dims_mps = initial_mps.get_tensor().dimensions();
    if(solver.result.meta.eigvals_found and solver.result.meta.eigvecsR_found) {
        auto eigvecs  = eig::view::get_eigvecs<cplx>(solver.result, eig::Side::R, converged_only);
        auto eigvals  = eig::view::get_eigvals<real>(solver.result, converged_only);

        if(eigvecs.cols() == eigvals.size()) /* Checks if eigenvectors converged for each eigenvalue */ {
            double overlap_sq_sum = 0;
            size_t num_solutions  = 0;
            // Eigenvalues are normally sorted small to large, so we reverse when looking for large.
            for(const auto &idx : indices.value_or(num::range<long>(0, eigvals.size()))) {
                if(idx >= eigvals.size()) throw except::logic_error("idx ({}) >= eigvals.size() ({})", idx, eigvals.size());
                results.emplace_back(opt_mps());
                auto &mps           = results.back();
                mps.is_basis_vector = true;
                if constexpr(settings::debug) tools::log->trace("Extracting result: idx {} | eigval {:.16f}", idx, eigvals(idx));
                mps.set_name(fmt::format("eigenvector {} [{:^8}]", idx, solver.config.tag));
                mps.set_tensor(eigvecs.col(idx).normalized(), dims_mps); // eigvecs are not always well normalized when we get them from eig::solver
                mps.set_sites(initial_mps.get_sites());
                mps.set_energy_shift(initial_mps.get_energy_shift()); // Will set energy if also given the eigval
                mps.set_overlap(std::abs(initial_mps.get_vector().dot(mps.get_vector())));
                mps.set_length(initial_mps.get_length());
                mps.set_time(solver.result.meta.time_total);
                mps.set_op(safe_cast<size_t>(solver.result.meta.num_op));
                mps.set_mv(safe_cast<size_t>(solver.result.meta.num_mv));
                mps.set_pc(safe_cast<size_t>(solver.result.meta.num_pc));
                mps.set_iter(safe_cast<size_t>(solver.result.meta.iter));
                mps.set_eigs_idx(idx);
                mps.set_eigs_nev(solver.result.meta.nev_converged);
                mps.set_eigs_ncv(solver.result.meta.ncv);
                mps.set_eigs_tol(solver.result.meta.tol);
                mps.set_eigs_eigval(eigvals[idx]);
                mps.set_eigs_ritz(solver.result.meta.ritz);
                mps.set_eigs_shift(solver.result.meta.sigma);
                mps.set_optcost(meta.optCost);
                mps.set_optalgo(meta.optAlgo);
                mps.set_optsolver(meta.optSolver);
                if(solver.result.meta.residual_norms.empty())
                    mps.set_rnorm(tools::finite::measure::residual_norm(mps.get_tensor(), tensors.get_multisite_mpo_squared(),
                                                                        tensors.get_multisite_env_var_blk().L, tensors.get_multisite_env_var_blk().R));
                else
                    mps.set_rnorm(solver.result.meta.residual_norms.at(safe_cast<size_t>(idx))); // primme convergence precision
                auto   measurements = MeasurementsTensorsFinite();
                double energy       = tools::finite::measure::energy(mps.get_tensor(), tensors, meta.svd_cfg, &measurements);
                double eigval       = energy - initial_mps.get_energy_shift();
                double variance     = tools::finite::measure::energy_variance(mps.get_tensor(), tensors, meta.svd_cfg, &measurements);

                mps.set_energy(energy);
                mps.set_eshift_eigval(eigval);
                mps.set_variance(variance);
                // tools::log->info("extract_results: set variance: {:.16f}", variance);

                mps.validate_basis_vector();

                // Sum up the contributions. Since full diag gives an orthonormal basis, this adds up to one. Normally only
                // a few eigenvectors contribute to most of the sum.
                overlap_sq_sum += mps.get_overlap() * mps.get_overlap();
                num_solutions++; // Count the number of solutions added
            }
        }
    }
}

void tools::finite::opt::internal::extract_results_subspace(const TensorsFinite &tensors, const opt_mps &initial_mps, const OptMeta &meta,
                                                            const eig::solver &solver, const std::vector<opt_mps> &subspace_mps,
                                                            std::vector<opt_mps> &results) {
    auto t_ext    = tid::tic_scope("extract");
    auto dims_mps = initial_mps.get_tensor().dimensions();
    if(solver.result.meta.eigvals_found and solver.result.meta.eigvecsR_found) {
        auto eigvecs = eig::view::get_eigvecs<cplx>(solver.result, eig::Side::R);
        auto eigvals = eig::view::get_eigvals<real>(solver.result);
        if(eigvecs.cols() == eigvals.size()) /* Checks if eigenvectors converged for each eigenvalue */ {
            auto indices = num::range<long>(0, eigvals.size());
            // Eigenvalues are normally sorted small to large, so we reverse when looking for large.
            if(meta.optRitz == OptRitz::LR) std::reverse(indices.begin(), indices.end());
            for(const auto &idx : indices) {
                results.emplace_back(opt_mps());
                auto &mps           = results.back();
                mps.is_basis_vector = false;
                mps.set_name(fmt::format("eigenvector {} [{:^8}]", idx, solver.config.tag));
                // eigvecs are not always well normalized when we get them from eig::solver
                mps.set_tensor(subspace::get_vector_in_fullspace(subspace_mps, eigvecs.col(idx).normalized()), dims_mps);
                mps.set_sites(initial_mps.get_sites());
                mps.set_energy_shift(initial_mps.get_energy_shift()); // Will set energy if also given the eigval
                mps.set_overlap(std::abs(initial_mps.get_vector().dot(mps.get_vector())));
                mps.set_length(initial_mps.get_length());
                mps.set_time(solver.result.meta.time_total);
                mps.set_op(safe_cast<size_t>(solver.result.meta.num_op));
                mps.set_mv(safe_cast<size_t>(solver.result.meta.num_mv));
                mps.set_pc(safe_cast<size_t>(solver.result.meta.num_pc));
                mps.set_iter(safe_cast<size_t>(solver.result.meta.iter));
                mps.set_eigs_idx(idx);
                mps.set_eigs_nev(solver.result.meta.nev_converged);
                mps.set_eigs_ncv(solver.result.meta.ncv);
                mps.set_eigs_tol(solver.result.meta.tol);
                mps.set_eigs_eigval(eigvals[idx]);
                mps.set_eigs_ritz(solver.result.meta.ritz);
                mps.set_eigs_shift(solver.result.meta.sigma);
                mps.set_optcost(meta.optCost);
                mps.set_optalgo(meta.optAlgo);
                mps.set_optsolver(meta.optSolver);
                if(solver.result.meta.residual_norms.empty())
                    mps.set_rnorm(tools::finite::measure::residual_norm(mps.get_tensor(), tensors.get_multisite_mpo_squared(),
                                                                        tensors.get_multisite_env_var_blk().L, tensors.get_multisite_env_var_blk().R));
                else
                    mps.set_rnorm(solver.result.meta.residual_norms.at(safe_cast<size_t>(idx))); // primme convergence precision
                auto   measurements = MeasurementsTensorsFinite();
                double energy       = tools::finite::measure::energy(mps.get_tensor(), tensors, meta.svd_cfg, &measurements);
                double eigval       = energy - initial_mps.get_energy_shift();
                double variance     = tools::finite::measure::energy_variance(mps.get_tensor(), tensors, meta.svd_cfg, &measurements);

                mps.set_energy(energy);
                mps.set_eshift_eigval(eigval);
                mps.set_variance(variance);
                // tools::log->info("extract_results_subspace: set variance: {:.16f}", variance);
                //                mps.set_grad_max(grad_max);
                mps.validate_basis_vector();
            }
        }
    }
}
