#include "math/tenx.h"
// -- (textra first)
#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "math/eig.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/state/StateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt/bfgs_callback.h"
#include "tools/finite/opt/bfgs_subspace_functor.h"
#include "tools/finite/opt/opt-internal.h"
#include "tools/finite/opt/report.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"
#include <ceres/gradient_problem.h>

using namespace tools::finite::opt;
using namespace tools::finite::opt::internal;

namespace settings {
    static constexpr bool debug_subspace = false;
}

opt_mps tools::finite::opt::internal::eigs_optimize_subspace(const TensorsFinite &tensors, const opt_mps &initial_mps,
                                                             [[maybe_unused]] const AlgorithmStatus &status, OptMeta &meta) {
    tools::log->trace("Optimizing subspace");
    auto t_sub = tid::tic_scope("subspace");
    initial_mps.validate_initial_mps();
    if(meta.optMode != OptMode::SUBSPACE)
        throw except::runtime_error("Wrong optimization mode [{}]. Expected [{}]", enum2sv(meta.optMode), enum2sv(OptMode::SUBSPACE));

    /*
     * Subspace optimization
     *
     *
     * In subspace optimization we consider a local set of sites l of the L-site system,
     * usually this corresponds to l = 1 or up to l = 8 adjacent sites.
     *
     * A subspace in this context is a truncated basis, i.e. a small subset of eigenvectors
     * of the local "effective" Hamiltonian. The subspace is a set of k eigenvectors [x]
     * which have significant overlap with the current state |y⟩, i.e. we define k such that
     *
     *          ε > 1 - Σ_i^k |⟨x_i|y⟩|²,
     *
     * where ε is a small number that controls error of the truncation. A value ε ~1e-10 is
     * reasonable. Note that the truncation implies that k is smaller than the local
     * Hilbert space dimension.
     *
     * After having found a good subspace, the objective is to find a linear combination
     * of eigenvectors which minimizes the energy variance.
     *
     * It is worth noting some observations. Let {x} be the set of all eigenvectors
     * to the local effective Hamiltonian H_local.
     * Then, when the DMRG process is fully converged:
     *      - only one x_i has overlap ⟨x_i|y⟩ = 1
     *        Since the sum of all overlaps must add to 1, the rest have <x_j|y> = 0 when i != j.
     *      - This x is also the one that minimizes the energy variance.
     *
     * However, before the DMRG process has converged this is not true. Instead:
     *      - we have ⟨x_i|y⟩ > 0 for several i.
     *      - a linear combination of several x can have lower variance than any
     *        single x.
     *
     * Fully diagonalizing H_local yields all K eigenvectors {x}, but if H_local is too big this operation
     * becomes prohibitively expensive. Instead we resort to finding a subset with k << K eigenvectors [x],
     * whose eigenvalues are the k energies closest to the current energy. Usually the eigenvectors
     * which have some overlap ⟨x_i|y⟩ > 0 are found in the subset [x] if k is large enough.
     *
     * Subspace optimization steps
     *
     * Step 1)  Find a subspace [x], i.e. take a set of k eigenvectors of the local effective Hamiltonian.
     *          Empirically, eigenvectors whose eigenvalues (energy) are closest to the current energy,
     *          tend to have nonzero overlap with the current vector |y⟩.
     *          On an iterative solver we keep increasing "nev" (number of requested eigenvectors) until
     *          the subspace error ε is small enough.
     *          If any eigenvectors have to removed, (e.g. due to memory/performance constraints),
     *          then sort the eigenvectors in order of decreasing overlap ⟨x_i|y⟩, and start deleting
     *          from the end.
     *
     * Step 2)  Project the squared effective K*K Hamiltonian, down to the k*k subspace, H².
     *          Using BFGS, find the linear combination |w⟩ of eigenvectors that minimizes the variance.
     *
     *              min_w Var H = ⟨H²⟩ - ⟨H⟩² = ⟨w|H²|w⟩ - ⟨E⟩²
     *
     *          where E are the energy eigenvalues from step 1.
     *
     */

    // Handy references
    const auto &model = *tensors.model;
    const auto &edges = *tensors.edges;

    /*
     *  Step 1) Find the subspace.
     *  The subspace is a set of eigenstates obtained from full or partial diagonalization
     */

    std::vector<opt_mps> subspace;
    switch(meta.optType) {
        case OptType::CPLX: subspace = internal::subspace::find_subspace<cplx>(tensors, settings::precision::target_subspace_error, meta); break;
        case OptType::REAL: subspace = internal::subspace::find_subspace<real>(tensors, settings::precision::target_subspace_error, meta); break;
    }

    tools::log->trace("Subspace found with {} eigenvectors", subspace.size());

    /*
     * Filter the eigenvectors
     *
     */

    internal::subspace::filter_subspace(subspace, settings::precision::max_subspace_size);

    /*
     *
     * Step 2) Optimize variance in the subspace of k eigenvectors
     *
     */

    // We need the eigenvalues in a convenient format as well
    auto eigvals = internal::subspace::get_eigvals(subspace);

    Eigen::VectorXcd subspace_vector = internal::subspace::get_vector_in_subspace(subspace, initial_mps.get_vector());
    tools::log->trace("Starting EIGS optimization of subspace");

    /*
     *
     *  Start the EIGS optimization process for the subspace
     *
     */
    std::vector<opt_mps> results;
    eig::solver          solver;

    auto options = internal::bfgs_default_options;
    auto summary = ceres::GradientProblemSolver::Summary();
    auto t_eigs  = tid::tic_scope("eigs");
    switch(meta.optType) {
        case OptType::CPLX: {
            tools::log->trace("Running EIGS subspace cplx");
            auto H2_subspace = subspace::get_hamiltonian_squared_in_subspace<cplx>(model, edges, subspace);
            solver.eig<eig::Form::SYMM>(H2_subspace.data(), H2_subspace.rows(), 'I', 1, 1, 0.0, 1.0, 1);
            break;
        }
        case OptType::REAL: {
            tools::log->trace("Running EIGS subspace real");
            auto H2_subspace = subspace::get_hamiltonian_squared_in_subspace<real>(model, edges, subspace);
            solver.eig<eig::Form::SYMM>(H2_subspace.data(), H2_subspace.rows(), 'I', 1, 1, 0.0, 1.0, 1);
            break;
        }
    }
    eigs_extract_results_subspace(tensors, initial_mps, meta, solver, subspace, results);

    if(results.empty()) {
        meta.optExit = OptExit::FAIL_ERROR;
        return initial_mps; // The solver failed
    }

    if constexpr(settings::debug or settings::debug_subspace) {
        for(const auto &optimized_mps : results) {
            auto t_dbg = tid::tic_scope("debug");
            // Check that Ceres results are correct
            double energy_check   = tools::finite::measure::energy(optimized_mps.get_tensor(), tensors);
            double variance_check = tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors);
            if(std::abs(1.0 - std::abs(optimized_mps.get_energy() / energy_check)) > 1e-3)
                tools::log->warn("Energy mismatch: Ceres: {:.16f} | DMRG {:.16f}", optimized_mps.get_energy(), energy_check);
            if(std::abs(1.0 - std::abs(optimized_mps.get_variance() / variance_check)) > 1e-3)
                tools::log->warn("Variance mismatch: Ceres: {:8.2e} | DMRG {:8.2e}", optimized_mps.get_variance(), variance_check);
        }
    }

    if(results.size() >= 2) {
        std::sort(results.begin(), results.end(), Comparator(meta)); // Smallest eigenvalue (i.e. variance) wins
    }

    for(const auto &mps : results) reports::eigs_add_entry(mps, spdlog::level::debug);
    return results.front();
}
