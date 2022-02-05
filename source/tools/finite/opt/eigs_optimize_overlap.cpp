#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/state/StateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt/opt-internal.h"
#include "tools/finite/opt/report.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"
#include <ceres/gradient_problem.h>

using namespace tools::finite::opt;
using namespace tools::finite::opt::internal;

opt_mps tools::finite::opt::internal::eigs_optimize_overlap(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status,
                                                            OptMeta &meta) {
    tools::log->trace("Optimizing in OVERLAP mode");
    auto t_sub = tid::tic_scope("eigs-ovl");
    initial_mps.validate_basis_vector();

    /*
     * Overlap optimization
     *
     * Introduction
     * In subspace optimization we consider a "local" subsection of the "global" L-site system,
     * usually this corresponds to n=2 local sites but in multisite dmrg we may consider
     * any n=[2,8] adjacent sites.
     *
     * The point here is to minimize the global energy variance Var H_global, by tuning only the parameters of
     * the local wavefunction |psi>_local corresponding to the local Hamiltonian H_local for n sites.
     * In other words, we seek the local vector which minimizes the global energy variance.
     *
     * It is worth noting some properties which hold when the DMRG process is fully converged:
     *      - If {x} are all the eigenvectors of H_local, only one of them has overlap <x_i|y> = 1
     *        where y is the current local vector for n sites.
     *        Since the sum of all overlaps must add to 1, the rest have <x_j|y> = 0 when i != j.
     *      - This particular eigenvector is also the one that minimizes the Var H_global,
     *        and no further optimization should be needed.
     *
     * However, before the DMRG process has converged this is not true. Instead:
     *      - If {x} are all the eigenvectors of H_local we have <x_i|y> > 0 for several i.
     *      - None {x} minimizes Var H_global, but a linear combination of several x could.
     *
     * Fully diagonalizing H_local yields all K eigenvectors {x}, but if H_local is too big this operation
     * becomes prohibitively expensive. Instead we resort to finding a subset with k << K eigenvectors [x],
     * whose eigenvalues are the k energies closest to the current energy. Usually the eigenvectors
     * which have some overlap <x_i|y> > 0 are found in the subset [x] if k is large enough.
     *
     * In Subspace Optimization we find a linear combination eigenvectors in the subspace [x] which
     * minimizes Var H_global, hence the name.
     *
     *
     *
     * Subspace optimization steps
     *
     * Step 0)  Find a subspace, i.e. a set of k eigenvectors [x] to the local Hamiltonian H_local with
     *          energy eigenvalue closest to the current energy. Note that H_local is a K * K matrix,
     *          and k << K. The set [x] is sorted in order of descending overlap <x_i|y>,
     *          where y is the current vector.
     *
     *
     * Step 1)
     *      - (O) If  OptMode::OVERLAP:
     *            Find the index for the best eigvec inside of the energy window:
     *              - OA) idx >= 0: eigvec found in window. Return that eigvec
     *              - OB) idx = -1: No eigvec found in window. Return old tensor
     *       -(V) If OptMode::SUBSPACE
     *          Make sure k is as small as possible. I.e. filter out eigenvectors from [x] down
     *          to an even smaller set of "relevant" eigvecs for doing subspace optimization.
     *          Allowing a maximum of k == 64 eigvecs keeps ram below 2GB when the problem
     *          size == 4096 (the linear size of H_local and the mps multisite_tensor).
     *          This means that we filter out
     *              * eigvecs outside of the energy window (if one is enabled)
     *              * eigvecs with little or no overlap to the current state.
     *          The filtering process collects the eigvec state with best overlap until
     *          either of these two criteria is met:
     *              * More than N=max_accept eigvecs have been collected
     *              * The subspace error "eps" is low enough¹
     *          The filtering returns the eigvecs sorted in decreasing overlap.
     *          ¹ We define eps = 1 - Σ_i |<x_i|y>|². A value eps ~1e-10 is reasonable.
     *
     * Step 2)
     *          Find the index for the best eigvec inside of the energy window:
     *              - VA) idx = -1: No eigvec found in window. Return old tensor
     *              - VB) We have many eigvecs, including the current tensor.
     *                    They should be viewed as good starting guesses for LBFGS optimization.
     *                    Optimize them one by one, and keep the
     *
     * Step 2)  Find the best overlapping state among the relevant eigvecs.
     * Step 3)  We can now make different decisions based on the overlap.
     *          A)  If best_overlap_idx == -1
     *              No state is in energy window -> discard! Return old multisite_mps_tensor.
     *          B)  If overlap_high <= best_overlap.
     *              This can happen if the environments have been modified just slightly since the last time considered
     *              these sites, but the signal is still clear -- we are still targeting the same state.
     *              However we can't be sure that the contributions from nearby states is just noise. Instead of just
     *              keeping the state we should optimize its variance. This is important in the later stages when variance
     *              is low and we don't want to ruin those last decimals.
     *              We just need to decide which initial guess to use.
     *                  B1) If best_overlap_variance <= theta_variance: set theta_initial = best_overlap_theta.
     *                  B2) Else, set theta_initial = multisite_mps_tensor.
     *          C)  If overlap_cat <= best_overlap and best_overlap < overlap_high
     *              This can happen for one reasons:
     *                  1) There are a few eigvec states with significant overlap (superposition)
     *              It's clear that we need to optimize, but we have to think carefully about the initial guess.
     *              Right now it makes sense to always choose best overlap theta, since that forces the algorithm to
     *              choose a particular state and not get stuck in superposition. Choosing the old theta may just entrench
     *              the algorithm into a local minima.
     *          D)  If 0 <= best_overlap and best_overlap < overlap_cat
     *              This can happen for three reasons, most often early in the simulation.
     *                  1) There are several eigvec states with significant overlap (superposition)
     *                  2) The highest overlapping states were outside of the energy window, leaving just these eigvecs.
     *                  3) The energy targeting of states has failed for some reason, perhaps the spectrum is particularly dense_lu.
     *              In any case, it is clear we are lost Hilbert space.
     *              Also, the subspace_error is no longer a good measure of how useful the subspace is to us, since it's only
     *              measuring how well the old state can be described, but the old state is likely very different from what
     *              we're looking for.
     *              So to address all three cases, do LBFGS optimization with best_overlap_theta as initial guess.
     *
     * In particular, notice that we never use the eigvec that happens to have the best variance.
     */
    // Handy references
    const auto &state = *tensors.state;
    const auto &model = *tensors.model;
    const auto &edges = *tensors.edges;

    /*
     *  Step 0) Find the subspace.
     *  The subspace set of eigenstates obtained from full or partial diagonalization
     */

    std::vector<opt_mps> subspace;
    switch(meta.optType) {
        case OptType::CPLX: subspace = internal::subspace::find_subspace<cplx>(tensors, settings::precision::target_subspace_error, meta); break;
        case OptType::REAL: subspace = internal::subspace::find_subspace<double>(tensors, settings::precision::target_subspace_error, meta); break;
    }

    tools::log->trace("eigs_optimize_overlap: subspace found with {} eigenvectors", subspace.size());

    auto max_overlap_idx = internal::subspace::get_idx_to_eigvec_with_highest_overlap(subspace, status.energy_llim_per_site, status.energy_ulim_per_site);
    if(max_overlap_idx) {
        auto &eigvec_max_overlap = *std::next(subspace.begin(), static_cast<long>(max_overlap_idx.value()));
        eigvec_max_overlap.set_variance(tools::finite::measure::energy_variance(eigvec_max_overlap.get_tensor(), tensors));
        if(tools::log->level() == spdlog::level::trace) {
            tools::log->trace("eigs_optimize_overlap: eigvec {:<2} has highest overlap {:.16f} | energy {:>20.16f} | variance {:>8.2e}",
                              max_overlap_idx.value(), eigvec_max_overlap.get_overlap(), eigvec_max_overlap.get_energy_per_site(),
                              eigvec_max_overlap.get_variance());
        }
        if(eigvec_max_overlap.get_overlap() < 0.1)
            tools::log->debug("eigs_optimize_overlap: Overlap fell below < 0.1: {:20.16f}", eigvec_max_overlap.get_overlap());
        return eigvec_max_overlap;
    } else {
        tools::log->warn("eigs_optimize_overlap: No overlapping states in energy range. Returning the state unchanged.");
        return initial_mps;
    }
}